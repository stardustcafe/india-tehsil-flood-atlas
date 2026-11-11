"""Cyclone wind: the 1-in-100-year tropical cyclone wind per sub-district, today and in a warmer climate.

Source: STORM synthetic tropical cyclones, 10 km wind-speed return periods for the North Indian basin
(Bloemendaal et al. 2020, 2022; 4TU.ResearchData, CC0). Present: 10,000 years of synthetic cyclones resampled
from IBTrACS 1980-2017. Future: the same with the 2015-2050 climate signal of four high-resolution HighResMIP
models (CMCC-CM2-VHR4, CNRM-CM6-1-HR, EC-Earth3P-HR, HadGEM3-GC31-HM) under high emissions (SSP5-8.5); the
median of the four, and each model for the range. Only the NI-basin GeoTIFFs are read out of the archives.

Winds are 10-minute sustained at 10 m height. IMD categories use 3-minute winds, about 5-10% higher, so the
classes here (17 / 24.5 / 33 / 46 m/s, IMD's cyclonic, severe, very severe and extremely severe thresholds)
are slightly conservative.

Per sub-district: area-weighted mean and maximum 100-year wind, and 2020 people (GPW / GHS-POP) where it reaches
33 m/s, very severe cyclonic storm strength.

Outputs:
  data/raw/storm/                       NI-basin GeoTIFFs
  data/processed/india_adm3_cyclone.csv one row per sub-district reached by cyclone winds
  viewer/data/cyclone.json              those sub-districts, TopoJSON index order
  viewer/data/wind_today/, wind_future/ tiles z4-8: 1-in-100-year wind classes
Usage: python scripts/14_cyclone_wind.py
"""
import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import requests
from rasterio.transform import from_origin

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import remote_zip  # noqa: E402

spec = importlib.util.spec_from_file_location("aq", ROOT / "scripts" / "10_aqueduct_floods.py")
aq = importlib.util.module_from_spec(spec)
spec.loader.exec_module(aq)

PRESENT = "https://data.4tu.nl/file/0ea98bdd-5772-4da8-ae97-99735e891aff/0b98c6f5-c7af-45bb-bb6e-dac53e1b8d55"
FUTURE = "https://data.4tu.nl/file/504c838e-2bd8-4d61-85a1-d495bdc560c3/"
ZIPS = {"present": (PRESENT, "STORM_FIXED_RETURN_PERIODS_NI_{rp}_YR_RP.tif"),
        "median": (FUTURE + "9740b235-c2fe-4677-a8b9-f90eb13e4938", "STORM_FIXED_RETURN_PERIODS_MEDIAN_NI_{rp}_YR_RP.tif"),
        "cmcc": (FUTURE + "2a3759e8-f4ad-4190-abb3-434ad5c4679b", "CMCC-CM2-VHR4_NI_{rp}_YR_RP.tif"),
        "cnrm": (FUTURE + "c460a0c8-f918-4ede-a734-75e77e99b102", "CNRM-CM6-1-HR_NI_{rp}_YR_RP.tif"),
        "ecearth": (FUTURE + "304d1441-bd71-47c7-8231-b20253c1cc2a", "EC-Earth3P-HR_NI_{rp}_YR_RP.tif"),
        "hadgem": (FUTURE + "856f9530-56d7-489e-8005-18ae36db4804", "HadGEM3-GC31-HM_NI_{rp}_YR_RP.tif")}
GCMS = ["cmcc", "cnrm", "ecearth", "hadgem"]
RPS = [10, 100]
VSCS = 33.0                          # m/s, very severe cyclonic storm
BREAKS = [17, 24.5, 33, 46]          # IMD cyclonic, severe, very severe, extremely severe (m/s)
COLORS = ["d5efe6", "8fd1bd", "3fa38a", "1d7a6b", "0b3f3a"]
RAW = ROOT / "data" / "raw" / "storm"
PROC = ROOT / "data" / "processed"
VIEW = ROOT / "viewer" / "data"
WORK = ROOT / "data" / "interim" / "wind_tiles"
W, N, RES = aq.W, aq.N, aq.RES
NROW, NCOL = round((aq.N - aq.S) / RES), round((aq.E - aq.W) / RES)


def fetch(key, rp):
    out = RAW / f"{key}_NI_{rp}.tif"
    if out.exists():
        return out
    url, pat = ZIPS[key]
    url = requests.get(url, headers={"Range": "bytes=0-0"}, timeout=60).url
    idx = remote_zip.index(url)
    name = next(k for k in idx if k.upper().endswith(pat.format(rp=rp).upper()) and "_NI_" in k)
    return remote_zip.extract(url, idx[name], out)


def on_grid(path):
    """STORM 0.1 deg wind -> the 30" analysis grid (nearest), 0 outside cyclone reach."""
    with rasterio.open(path) as d:
        a, t = d.read(1), d.transform
    lat = N - (np.arange(NROW) + 0.5) * RES
    lon = W + (np.arange(NCOL) + 0.5) * RES
    r = np.clip(((t.f - lat) / -t.e).astype(int), 0, a.shape[0] - 1)
    c = np.clip(((lon - t.c) / t.a).astype(int), 0, a.shape[1] - 1)
    g = a[r][:, c].astype(np.float32)
    g[~np.isfinite(g) | (g < 0)] = 0
    return g


def build_tiles(name, wind, india):
    env, res = aq.tiles.gdal_env(), 40075016.685578488 / (256 * 2 ** aq.tiles.MAXZ)
    WORK.mkdir(parents=True, exist_ok=True)
    cls = np.digitize(wind, BREAKS).astype(np.uint8)  # 0 below 17 m/s
    cls[~india] = 0
    src = WORK / f"{name}_cls.tif"
    prof = {"driver": "GTiff", "dtype": "uint8", "count": 1, "width": NCOL, "height": NROW, "nodata": 0,
            "crs": "EPSG:4326", "transform": from_origin(W, N, RES, RES), "compress": "deflate"}
    with rasterio.open(src, "w", **prof) as d:
        d.write(cls, 1)
    colors = WORK / "wind_colors.txt"
    colors.write_text("nv 0 0 0 0\n" + "".join(
        f"{k} {int(c[:2], 16)} {int(c[2:4], 16)} {int(c[4:], 16)} 255\n" for k, c in enumerate(COLORS[1:], 1)))
    merc, rgba, xyz = WORK / f"{name}_3857.tif", WORK / f"{name}_rgba.vrt", WORK / f"{name}_xyz"
    subprocess.run(["gdalwarp", "-q", "-overwrite", "-t_srs", "EPSG:3857", "-tr", str(res), str(res), "-tap",
                    "-r", "near", "-co", "COMPRESS=DEFLATE", str(src), str(merc)], check=True, env=env)
    subprocess.run(["gdaldem", "color-relief", "-q", "-alpha", "-exact_color_entry", "-of", "VRT",
                    str(merc), str(colors), str(rgba)], check=True, env=env)
    shutil.rmtree(xyz, ignore_errors=True)
    subprocess.run([str(aq.tiles.CONDA / "python.exe"), str(aq.tiles.CONDA / "Scripts" / "gdal2tiles.py"), "--xyz", "-z",
                    f"4-{aq.tiles.MAXZ}", "-r", "near", "-w", "none", "-x", "-q",
                    f"--processes={max(1, (os.cpu_count() or 2) // 2)}", str(rgba), str(xyz)], check=True, env=env)
    aq.tiles.bundler.bundle(xyz, VIEW / f"wind_{name}", aq.tiles.bundler.parse_palette(",".join(COLORS[1:])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-tiles", action="store_true")
    args = ap.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)
    grids = {(k, rp): on_grid(fetch(k, rp)) for k in ZIPS for rp in RPS if rp == 100 or k in ("present", "median")}

    adm3 = aq.adm3mod.load_adm3()
    n = len(adm3) + 1
    uid, cell, frac = aq.weights(adm3, NROW, NCOL)
    cover = np.bincount(cell, weights=frac, minlength=NROW * NCOL)
    share = frac / np.maximum(cover[cell], 1e-9)
    order_uid = adm3["uid"].values
    india = np.zeros(NROW * NCOL, bool); india[cell] = True
    lat = N - (np.arange(NROW) + 0.5) * RES
    area = np.repeat(np.cos(np.radians(lat)), NCOL)[cell] * frac
    tot_area = np.bincount(uid, weights=area, minlength=n)
    pops = {}
    for s in aq.SOURCES:
        with rasterio.open(PROC / f"{s}2020_india_pop.tif") as d:
            p = d.read(1).astype(np.float64).ravel()
        p[p < 0] = 0
        pops[s] = p[cell] * share

    def mean(g):
        return (np.bincount(uid, weights=g.ravel()[cell] * area, minlength=n) / np.maximum(tot_area, 1e-9))[order_uid]

    def peak(g):
        out = np.zeros(n)
        np.maximum.at(out, uid, g.ravel()[cell])
        return out[order_uid]

    def people(g, s):
        return np.bincount(uid, weights=pops[s] * (g.ravel()[cell] >= VSCS), minlength=n)[order_uid]

    now, fut = grids[("present", 100)], grids[("median", 100)]
    gcm_mean = np.stack([mean(grids[(k, 100)]) for k in GCMS])
    df = pd.DataFrame({"adm3_id": adm3["adm3_id"].values,
                       "wind100_today_mean": mean(now), "wind100_today_max": peak(now),
                       "wind100_future_mean": mean(fut), "wind100_future_max": peak(fut),
                       "wind100_future_lo": gcm_mean.min(0), "wind100_future_hi": gcm_mean.max(0),
                       "wind10_today_mean": mean(grids[("present", 10)]), "wind10_future_mean": mean(grids[("median", 10)])})
    for s in aq.SOURCES:
        df[f"{s}_vscs_today"] = people(now, s).round()
        df[f"{s}_vscs_future"] = people(fut, s).round()
    df = df[(df["wind100_today_max"] >= BREAKS[0]) | (df["wind100_future_max"] >= BREAKS[0])]
    df.round(2).to_csv(PROC / "india_adm3_cyclone.csv", index=False)

    g = pd.read_csv(PROC / "india_adm3_gfsm.csv")
    order = list(g.loc[g["qa_flag"] != "sliver", "adm3_id"])
    pos = {a: i for i, a in enumerate(order)}
    df = df[df["adm3_id"].isin(pos)].assign(i=lambda d: d["adm3_id"].map(pos)).sort_values("i")
    x10 = lambda c: [int(round(v * 10)) for v in df[c]]
    payload = {"note": "wind m/s x10, 10-minute sustained; people 2020", "vscs": VSCS, "idx": df["i"].tolist(),
               "m0": x10("wind100_today_mean"), "x0": x10("wind100_today_max"), "m1": x10("wind100_future_mean"),
               "x1": x10("wind100_future_max"), "lo": x10("wind100_future_lo"), "hi": x10("wind100_future_hi"),
               "t0": x10("wind10_today_mean"), "t1": x10("wind10_future_mean"),
               **{f"{s}{k}": [int(v) for v in df[f"{s}_vscs_{w}"]] for s in aq.SOURCES for k, w in ((0, "today"), (1, "future"))}}
    (VIEW / "cyclone.json").write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")

    print(f"{len(df)} sub-districts reached by cyclone-strength 100-year winds (>= {BREAKS[0]} m/s)")
    for s in aq.SOURCES:
        print(f"  {s}: people where the 100-year wind reaches {VSCS} m/s: today {df[f'{s}_vscs_today'].sum() / 1e7:.2f} crore,"
              f" future {df[f'{s}_vscs_future'].sum() / 1e7:.2f} crore")
    d = df["wind100_future_mean"] - df["wind100_today_mean"]
    print(f"  change in mean 100-year wind: median {d.median():+.2f} m/s, {(d > 0).mean():.0%} of sub-districts higher")
    top = df.nlargest(5, "wind100_today_max")[["adm3_id", "wind100_today_max", "wind100_future_max"]]
    print(top.to_string(index=False))
    print(f"{VIEW / 'cyclone.json'} {(VIEW / 'cyclone.json').stat().st_size / 1e6:.2f} MB")
    if not args.no_tiles:
        india = india.reshape(NROW, NCOL)
        for name, gr in (("today", now), ("future", fut)):
            build_tiles(name, gr, india)


if __name__ == "__main__":
    main()
