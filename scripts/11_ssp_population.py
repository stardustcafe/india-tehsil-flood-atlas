"""Projected population per sub-district under the SSPs, and people in the 1-in-100-year flood zone with
that population, 2030-2090.

Source: Wang, Meng & Long (2022), 1 km grid population 2020-2100 under SSP1-5, Scientific Data 9:563,
https://doi.org/10.6084/m9.figshare.19608594 (CC BY 4.0). Random-forest downscaling of the SSP national
populations from a WorldPop 2020 base. Only the needed years are read out of the 2.5 GB zips with HTTP range
requests, cropped to India and kept in data/raw/ssp_pop/india/.

Method: our 2020 figures (GPW / GHS-POP, scripts/06) stay the anchor; the SSP grids supply growth.
  growth   g = SSP people in the sub-district in year y / in 2020
  flood    f = SSP people inside the future flood zone in year y / inside the same zone in 2020
               (zone = cells where >= 3 of 5 river models or the coastal map flood, from scripts/10), so growth
               that concentrates on floodplains counts; falls back to g where the zone held < 50 SSP people
  people in the flood zone with projected population = people in the zone with 2020 population (scripts/10) x f
SSP and flood maps are paired by nearest RCP: SSP1-2.6 and SSP2-4.5 with RCP4.5, SSP3-7.0 and SSP5-8.5 with RCP8.5.

Outputs:
  data/processed/india_adm3_ssp_pop.csv  one row per sub-district x SSP x year
  viewer/data/popssp.json                growth (x1000) for 2030/2050/2070/2090 (climate period midpoints) and
                                         flood-zone people for 2030/2050/2080 (flood-map years), TopoJSON order
Usage: python scripts/11_ssp_population.py
"""
import importlib.util
import json
import struct
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import requests
from rasterio.transform import from_origin

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("aq", ROOT / "scripts" / "10_aqueduct_floods.py")
aq = importlib.util.module_from_spec(spec)
spec.loader.exec_module(aq)

FIGSHARE = {1: 34829160, 2: 34829370, 3: 45894312, 5: 34829391}  # file ids of SPP1/SSP2/SSP3/SSP5 zips
SSPS = {"ssp126": 1, "ssp245": 2, "ssp370": 3, "ssp585": 5}
RCP_OF = {"ssp126": "rcp4p5", "ssp245": "rcp4p5", "ssp370": "rcp8p5", "ssp585": "rcp8p5"}
G_YEARS = [2030, 2050, 2070, 2090]   # midpoints of the climate periods 2020-39 ... 2080-99
F_YEARS = [2030, 2050, 2080]         # Aqueduct flood-map years
YEARS = sorted({2020, *G_YEARS, *F_YEARS})
RAW = ROOT / "data" / "raw" / "ssp_pop" / "india"
PROC = ROOT / "data" / "processed"
VIEW = ROOT / "viewer" / "data"
W, N, RES = aq.W, aq.N, aq.RES
NROW, NCOL = round((aq.N - aq.S) / RES), round((aq.E - aq.W) / RES)


def ranged(fid, a, b, **kw):
    """GET bytes a..b of a figshare file; its signed storage URL expires after 10 s, so sign every request."""
    url = requests.head(f"https://ndownloader.figshare.com/files/{fid}", allow_redirects=False, timeout=60).headers["location"]
    r = requests.get(url, headers={"Range": f"bytes={a}-{b}"}, timeout=900, **kw)
    r.raise_for_status()
    return r


def zip_index(fid):
    """name -> (local header offset, compressed size, method) from the zip's central directory."""
    size = int(ranged(fid, 0, 0).headers["Content-Range"].split("/")[1])
    tail = ranged(fid, size - 65536, size - 1).content
    i = tail.rfind(b"PK\x05\x06")
    cdsize, cdoff = struct.unpack("<II", tail[i + 12:i + 20])
    cd = ranged(fid, cdoff, cdoff + cdsize - 1).content
    out, p = {}, 0
    while p < len(cd):
        f = struct.unpack("<IHHHHHHIIIHHHHHII", cd[p:p + 46])
        name = cd[p + 46:p + 46 + f[10]].decode()
        out[name] = (f[16], f[8], f[4])
        p += 46 + f[10] + f[11] + f[12]
    return out


def fetch(ssp, year):
    out = RAW / f"SSP{ssp}_{year}.tif"
    if out.exists():
        return out
    fid = FIGSHARE[ssp]
    idx = zip_index(fid)
    name = next(k for k in idx if k.endswith(f"SSP{ssp}_{year}.tif"))
    lho, csize, method = idx[name]
    h = ranged(fid, lho, lho + 29).content
    start = lho + 30 + sum(struct.unpack("<HH", h[26:30]))
    tmp = RAW.parent / f"SSP{ssp}_{year}.global.tif"
    dec = zlib.decompressobj(-15) if method == 8 else None
    with ranged(fid, start, start + csize - 1, stream=True) as r, open(tmp, "wb") as f:
        for ch in r.iter_content(1 << 20):
            f.write(dec.decompress(ch) if dec else ch)
        if dec:
            f.write(dec.flush())
    with rasterio.open(tmp) as d:
        # grid is offset by -0.05 cell from the 30" lattice: round to the nearest cell
        r0, c0 = round((d.transform.f - N) / RES), round((W - d.transform.c) / RES)
        a = d.read(1, window=((r0, r0 + NROW), (c0, c0 + NCOL))).astype(np.float32)
        a[(a == d.nodata) | ~np.isfinite(a) | (a < 0)] = 0
    prof = {"driver": "GTiff", "dtype": "float32", "count": 1, "width": NCOL, "height": NROW, "crs": "EPSG:4326",
            "transform": from_origin(W, N, RES, RES), "compress": "deflate", "predictor": 3, "tiled": True}
    with rasterio.open(out, "w", **prof) as dst:
        dst.write(a, 1)
    tmp.unlink()
    print(f"  SSP{ssp} {year}: {a.sum() / 1e7:.1f} crore people in the India window")
    return out


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    jobs = [(s, y) for s in SSPS.values() for y in YEARS]
    with ThreadPoolExecutor(4) as ex:  # ~120 MB compressed per year
        list(ex.map(lambda j: fetch(*j), jobs))

    adm3 = aq.adm3mod.load_adm3()
    uid, cell, frac = aq.weights(adm3, NROW, NCOL)
    n = len(adm3) + 1
    cover = np.bincount(cell, weights=frac, minlength=NROW * NCOL)
    share = frac / np.maximum(cover[cell], 1e-9)
    order_uid = adm3["uid"].values

    def tehsil_sum(grid, mask=None):
        v = grid.ravel()[cell] * share
        if mask is not None:
            v = v * mask.ravel()[cell]
        return np.bincount(uid, weights=v, minlength=n)[order_uid]

    def zone(rcp, year):
        wet = sum((aq.depth(aq.river_name(rcp, g, str(year))) > 0).astype(np.int8) for g in aq.GCMS)
        return (wet >= 3) | (aq.depth(aq.coast_name(rcp, str(year))) > 0)

    zones = {(r, y): zone(r, y) for r in set(RCP_OF.values()) for y in F_YEARS}
    fz = pd.read_csv(PROC / "india_adm3_floodzone.csv").set_index(["state", "adm3_id"])
    rows, view = [], {}
    for key, s in SSPS.items():
        grids = {y: rasterio.open(RAW / f"SSP{s}_{y}.tif").read(1) for y in YEARS}
        base = tehsil_sum(grids[2020])
        rec = {}
        for y in YEARS[1:]:
            tot = tehsil_sum(grids[y])
            g = np.where(base > 0, tot / np.maximum(base, 1e-9), 1.0)
            row = {"adm3_id": adm3["adm3_id"].values, "ssp": key, "year": y, "growth": g.round(4)}
            if y in F_YEARS:
                z = zones[(RCP_OF[key], y)]
                z0, z1 = tehsil_sum(grids[2020], z), tehsil_sum(grids[y], z)
                f = np.where(z0 >= 50, z1 / np.maximum(z0, 1e-9), g)
                row["flood_growth"] = f.round(4)
                aq_rows = fz.loc[f"{RCP_OF[key]}_{y}"].reindex(adm3["adm3_id"].values)
                for src in aq.SOURCES:
                    row[f"{src}_floodzone"] = (aq_rows[f"{src}_total"].values * f).round()
                    rec[f"{src}_f{y}"] = row[f"{src}_floodzone"]
            if y in G_YEARS:
                rec[f"g{y}"] = g
            rows.append(pd.DataFrame(row))
        view[key] = rec
        print(key, "India growth to 2050/2090: "
              f"{tehsil_sum(grids[2050]).sum() / base.sum():.3f} / {tehsil_sum(grids[2090]).sum() / base.sum():.3f}")
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(PROC / "india_adm3_ssp_pop.csv", index=False)

    g = pd.read_csv(PROC / "india_adm3_gfsm.csv")
    order = g.loc[g["qa_flag"] != "sliver", "adm3_id"]
    def arr(v, scale=1):
        x = pd.Series(np.asarray(v, float), index=adm3["adm3_id"].values).reindex(order).values
        return [0 if not np.isfinite(y) else int(round(y * scale)) for y in x]
    payload = {"g_years": G_YEARS, "f_years": F_YEARS, "rcp_of": RCP_OF, "note": "g = growth factor x1000",
               "data": {k: {f: arr(v, 1000 if f[1:].isdigit() else 1) for f, v in rec.items()} for k, rec in view.items()}}
    (VIEW / "popssp.json").write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")

    summ = out[out["year"].isin(F_YEARS)].groupby(["ssp", "year"])[[f"{s}_floodzone" for s in aq.SOURCES]].sum() / 1e7
    print(summ.round(2).to_string(), "\n(crore people in the flood zone with projected population)")
    print(f"{VIEW / 'popssp.json'} {(VIEW / 'popssp.json').stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
