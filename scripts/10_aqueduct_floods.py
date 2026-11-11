"""People and land in the 1-in-100-year river and coastal flood zone per sub-district, today and under
RCP4.5 / RCP8.5 in 2030, 2050 and 2080, from WRI Aqueduct Floods v2 (Ward et al. 2020; CC BY 4.0).

Source maps (https://aqueduct.wridata.org/AqueductFloods20/): global 30" inundation depth (m), no flood defences.
  river  inunriver_historical_000000000WATCH_1980_rp00100           baseline (1960-1999 climate)
         inunriver_{rcp}_{gcm}_{2030|2050|2080}_rp00100             5 CMIP5 models (GLOFRIS)
  coast  inuncoast_historical_nosub_hist_rp0100_0                   baseline
         inuncoast_{rcp}_wtsub_{year}_rp0100_0_perc_50              median sea-level rise, land subsidence included
Future periods are centred on the year: 2030 = 2010-2049, 2050 = 2030-2069, 2080 = 2060-2099.

Method: a 30" cell is in the flood zone when river or coastal depth > 0. Each cell's 2020 population
(data/processed/{gpw,ghs}2020_india_pop.tif, same grid) is shared among the sub-districts it overlaps, in proportion
to overlap area (10" rasterisation). For river flooding each model is overlaid separately; the sub-district value is
the median of the 5 models, with their minimum and maximum as the range. Population is held at 2020 levels.

Outputs:
  data/raw/aqueduct/india/*.tif           India window of each map, depth in cm (int16)
  data/processed/india_adm3_floodzone.csv one row per sub-district x state (today / rcp x year)
  viewer/data/flood.json                  compact arrays in TopoJSON order
  viewer/data/flood_{state}/              depth tiles z4-8 (median model river depth, max with coast)
Usage: python scripts/10_aqueduct_floods.py [--workers 6]
"""
import argparse
import importlib.util
import json
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import requests
from rasterio.features import rasterize
from rasterio.transform import from_origin
from rasterio.windows import from_bounds

ROOT = Path(__file__).resolve().parents[1]


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


adm3mod = load("adm3", "02_adm3_flood_susceptibility.py")
tiles = load("poptiles", "07_build_pop_tiles.py")  # gdal_env, CONDA, bundler

BASE_URL = "https://aqueduct.wridata.org/AqueductFloods20"
RAW = ROOT / "data" / "raw" / "aqueduct"
CROP = RAW / "india"
PROC = ROOT / "data" / "processed"
VIEW = ROOT / "viewer" / "data"
WORK = ROOT / "data" / "interim" / "flood_tiles"
W, S, E, N, RES, SUB = 68.0, 6.0, 98.0, 38.0, 1 / 120, 3  # same window and grid as the population rasters
GCMS = ["0000GFDL-ESM2M", "0000HadGEM2-ES", "00IPSL-CM5A-LR", "MIROC-ESM-CHEM", "00000NorESM1-M"]
RCPS = {"rcp4p5": "RCP4.5", "rcp8p5": "RCP8.5"}
YEARS = ["2030", "2050", "2080"]
SOURCES = ["gpw", "ghs"]
# depth classes (m): (0, 0.5], (0.5, 1], (1, 2], (2, 4], > 4
DEPTH_BREAKS = [0.5, 1, 2, 4]
COLORS = ["bfe6f0", "7cc4dc", "3a9ac2", "1c6aa0", "0b3d6e"]


def river_name(rcp, gcm, year):
    if rcp == "historical":
        return "inunriver_historical_000000000WATCH_1980_rp00100"
    return f"inunriver_{rcp}_{gcm}_{year}_rp00100"


def coast_name(rcp, year):
    if rcp == "historical":
        return "inuncoast_historical_nosub_hist_rp0100_0"
    return f"inuncoast_{rcp}_wtsub_{year}_rp0100_0_perc_50"


def all_names():
    names = [river_name("historical", None, None), coast_name("historical", None)]
    for rcp in RCPS:
        for year in YEARS:
            names += [river_name(rcp, g, year) for g in GCMS] + [coast_name(rcp, year)]
    return names


def fetch_crop(name):
    """Download one global map, keep the India window as int16 cm, drop the global file."""
    out = CROP / f"{name}.tif"
    if out.exists():
        return out
    tmp = RAW / f"{name}.part"
    for attempt in range(5):
        try:
            with requests.get(f"{BASE_URL}/{name}.tif", stream=True, timeout=600) as r:
                r.raise_for_status()
                with open(tmp, "wb") as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk)
            break
        except requests.RequestException as e:
            if attempt == 4:
                raise SystemExit(f"{name}: {e}")
    with rasterio.open(tmp) as d:
        win = from_bounds(W, S, E, N, transform=d.transform).round_offsets().round_lengths()
        a = d.read(1, window=win)
        nd = d.nodata
    a = np.where((a == nd) | ~np.isfinite(a) | (a < 0), 0, np.round(a * 100)).clip(0, 32767).astype(np.int16)
    prof = {"driver": "GTiff", "dtype": "int16", "count": 1, "width": a.shape[1], "height": a.shape[0], "nodata": None,
            "crs": "EPSG:4326", "transform": from_origin(W, N, RES, RES), "compress": "deflate", "predictor": 2, "tiled": True}
    with rasterio.open(out.with_suffix(".tmp.tif"), "w", **prof) as dst:
        dst.write(a, 1)
    out.with_suffix(".tmp.tif").replace(out)
    tmp.unlink()
    print(f"  {name}: {(a > 0).mean():.1%} of window flooded")
    return out


def depth(name):
    with rasterio.open(CROP / f"{name}.tif") as d:
        return d.read(1)


def weights(adm3, nrow, ncol):
    """(uid, cell, overlap fraction) for every sub-district / 30" cell pair, from a 10" rasterisation in bands."""
    keys, counts = [], []
    band = 384
    for r0 in range(0, nrow, band):
        nr = min(band, nrow - r0)
        t = from_origin(W, N - r0 * RES, RES / SUB, RES / SUB)
        ids = rasterize(zip(adm3.geometry, adm3["uid"]), out_shape=(nr * SUB, ncol * SUB), transform=t, fill=0, dtype="int32")
        cell = ((np.arange(nr * SUB)[:, None] // SUB + r0) * ncol + np.arange(ncol * SUB)[None, :] // SUB)
        m = ids > 0
        k, c = np.unique(ids[m].astype(np.int64) * (nrow * ncol) + cell[m], return_counts=True)
        keys.append(k); counts.append(c)
    k, c = np.concatenate(keys), np.concatenate(counts)
    uid, cell, frac = k // (nrow * ncol), k % (nrow * ncol), c / SUB ** 2
    missing = set(adm3["uid"]) - set(np.unique(uid))
    if missing:  # smaller than a 10" cell: representative point
        rp = adm3[adm3["uid"].isin(missing)]
        pts = rp.geometry.representative_point()
        r = ((N - pts.y) / RES).astype(int).clip(0, nrow - 1)
        cc = ((pts.x - W) / RES).astype(int).clip(0, ncol - 1)
        uid = np.concatenate([uid, rp["uid"].values]); cell = np.concatenate([cell, (r * ncol + cc).values])
        frac = np.concatenate([frac, np.full(len(rp), 0.01)])
        print(f"{len(missing)} sub-districts smaller than 10\" use their representative point")
    return uid, cell, frac


def build_tiles(name, dep, india):
    """Depth classes -> EPSG:3857 -> colour relief -> gdal2tiles z4-8 -> palette bundles."""
    env, res = tiles.gdal_env(), 40075016.685578488 / (256 * 2 ** tiles.MAXZ)
    WORK.mkdir(parents=True, exist_ok=True)
    cls = np.digitize(dep / 100, DEPTH_BREAKS, right=True).astype(np.uint8) + 1
    cls[(dep <= 0) | ~india] = 0
    src = WORK / f"{name}_cls.tif"
    prof = {"driver": "GTiff", "dtype": "uint8", "count": 1, "width": dep.shape[1], "height": dep.shape[0], "nodata": 0,
            "crs": "EPSG:4326", "transform": from_origin(W, N, RES, RES), "compress": "deflate"}
    with rasterio.open(src, "w", **prof) as d:
        d.write(cls, 1)
    colors = WORK / "flood_colors.txt"
    colors.write_text("nv 0 0 0 0\n" + "".join(
        f"{k} {int(c[:2], 16)} {int(c[2:4], 16)} {int(c[4:], 16)} 255\n" for k, c in enumerate(COLORS, 1)))
    merc, rgba, xyz = WORK / f"{name}_3857.tif", WORK / f"{name}_rgba.vrt", WORK / f"{name}_xyz"
    subprocess.run(["gdalwarp", "-q", "-overwrite", "-t_srs", "EPSG:3857", "-tr", str(res), str(res), "-tap",
                    "-r", "near", "-co", "COMPRESS=DEFLATE", str(src), str(merc)], check=True, env=env)
    subprocess.run(["gdaldem", "color-relief", "-q", "-alpha", "-exact_color_entry", "-of", "VRT",
                    str(merc), str(colors), str(rgba)], check=True, env=env)
    shutil.rmtree(xyz, ignore_errors=True)
    subprocess.run([str(tiles.CONDA / "python.exe"), str(tiles.CONDA / "Scripts" / "gdal2tiles.py"), "--xyz", "-z",
                    f"4-{tiles.MAXZ}", "-r", "near", "-w", "none", "-x", "-q",
                    f"--processes={max(1, (os.cpu_count() or 2) // 2)}", str(rgba), str(xyz)], check=True, env=env)
    tiles.bundler.bundle(xyz, VIEW / f"flood_{name}", tiles.bundler.parse_palette(",".join(COLORS)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--no-tiles", action="store_true")
    args = ap.parse_args()
    CROP.mkdir(parents=True, exist_ok=True)
    names = all_names()
    print(f"{sum(not (CROP / f'{n}.tif').exists() for n in names)} of {len(names)} maps to download")
    with ThreadPoolExecutor(args.workers) as ex:
        list(ex.map(fetch_crop, names))

    adm3 = adm3mod.load_adm3()
    nrow, ncol = round((N - S) / RES), round((E - W) / RES)
    uid, cell, frac = weights(adm3, nrow, ncol)
    n = len(adm3) + 1
    india = np.zeros(nrow * ncol, bool); india[cell] = True
    india = india.reshape(nrow, ncol)
    # share of each cell's people and area going to each pair (cells cut by the border or coast: India's part only)
    cover = np.bincount(cell, weights=frac, minlength=nrow * ncol)
    lat = N - (np.arange(nrow) + 0.5) * RES
    cell_km2 = (6371.0072 ** 2 * np.radians(RES) * np.abs(np.sin(np.radians(lat + RES / 2)) - np.sin(np.radians(lat - RES / 2))))
    area = np.repeat(cell_km2, ncol)[cell] * frac
    tot_area = np.bincount(uid, weights=area, minlength=n)
    pops = {}
    for s in SOURCES:
        with rasterio.open(PROC / f"{s}2020_india_pop.tif") as d:
            assert d.width == ncol and d.height == nrow and abs(d.transform.c - W) < 1e-9 and abs(d.transform.f - N) < 1e-9
            p = d.read(1).astype(np.float64).ravel()
        p[p < 0] = 0
        pops[s] = p[cell] * frac / np.maximum(cover[cell], 1e-9)
    order_uid = adm3["uid"].values

    def zone_stats(mask):
        f = mask.ravel()[cell]
        out = {s: np.bincount(uid, weights=pops[s] * f, minlength=n)[order_uid] for s in SOURCES}
        out["area_pct"] = (100 * np.bincount(uid, weights=area * f, minlength=n) / np.maximum(tot_area, 1e-9))[order_uid]
        return out

    states = [("today", None, None)] + [(f"{r}_{y}", r, y) for r in RCPS for y in YEARS]
    rows, view = [], {}
    for key, rcp, year in states:
        print(key)
        coast = depth(coast_name(rcp or "historical", year)) > 0
        rivers = [depth(river_name("historical", None, None))] if rcp is None else [depth(river_name(rcp, g, year)) for g in GCMS]
        per = []
        for rd in rivers:  # each model separately, then median / min / max across models
            r = zone_stats(rd > 0)
            t = zone_stats((rd > 0) | coast)
            per.append((r, t))
        c = zone_stats(coast)
        stack = lambda part, f: np.stack([p[part][f] for p in per])
        rec = {"adm3_id": adm3["adm3_id"].values, "state": key, "scenario": RCPS.get(rcp, "today"), "year": year or "today",
               "area_pct": np.median(stack(1, "area_pct"), 0).round(2)}
        for s in SOURCES:
            tt = stack(1, s)
            rec[f"{s}_total"] = np.median(tt, 0).round()
            rec[f"{s}_total_min"] = tt.min(0).round()
            rec[f"{s}_total_max"] = tt.max(0).round()
            rec[f"{s}_river"] = np.median(stack(0, s), 0).round()
            rec[f"{s}_coast"] = c[s].round()
        rows.append(pd.DataFrame(rec))
        view[key] = rec
        if not args.no_tiles:
            md = np.median(np.stack(rivers), 0)
            build_tiles(key, np.maximum(md, depth(coast_name(rcp or "historical", year))).astype(np.int16), india)

    out = pd.concat(rows, ignore_index=True)
    out.to_csv(PROC / "india_adm3_floodzone.csv", index=False)

    g = pd.read_csv(PROC / "india_adm3_gfsm.csv")
    order = g.loc[g["qa_flag"] != "sliver", "adm3_id"]
    def arr(v, scale=1):
        x = pd.Series(np.asarray(v, float), index=adm3["adm3_id"].values).reindex(order).values
        return [0 if not np.isfinite(y) else int(round(y * scale)) for y in x]
    payload = {"rp": 100, "years": YEARS, "rcps": RCPS, "sources": SOURCES, "depth_breaks": DEPTH_BREAKS, "colors": COLORS,
               "note": "people in 2020; area_pct x10", "data": {}}
    for key, rec in view.items():
        d = {"a": arr(rec["area_pct"], 10)}
        for s in SOURCES:
            d[s] = arr(rec[f"{s}_total"]); d[s + "R"] = arr(rec[f"{s}_river"]); d[s + "C"] = arr(rec[f"{s}_coast"])
            if key != "today":
                d[s + "Lo"] = arr(rec[f"{s}_total_min"]); d[s + "Hi"] = arr(rec[f"{s}_total_max"])
        payload["data"][key] = d
    (VIEW / "flood.json").write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")

    summ = out.groupby("state")[[f"{s}_total" for s in SOURCES] + [f"{s}_total_min" for s in SOURCES] + [f"{s}_total_max" for s in SOURCES]].sum()
    print((summ / 1e7).round(2).to_string(), "\n(crore people)")
    print(f"{VIEW / 'flood.json'} {(VIEW / 'flood.json').stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
