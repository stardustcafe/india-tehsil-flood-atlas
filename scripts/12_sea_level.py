"""Sea-level rise for coastal sub-districts, and land and people below the future mean sea level.

Sea level: IPCC AR6 sea-level projections (Fox-Kemper et al. 2021; Garner et al. 2021, Zenodo 6382554,
CC BY 4.0), relative sea-level change at the 1 deg grid points, medium confidence, against 1995-2014, for
SSP1-2.6 / 2-4.5 / 3-7.0 / 5-8.5. Median and likely range (17th-83rd percentile). Tide-gauge locations are
left out: their local land-motion terms are inconsistent within a few km (e.g. 4 cm vs 2 m by 2100 around
Kolkata). Each coastal sub-district takes the grid point nearest its representative point.

Elevation: DeltaDTM v1.1 (Pronk et al. 2024, CC BY 4.0, built from Copernicus GLO-30 (c) DLR e.V.), a 30 m
coastal terrain model (buildings and trees removed) relative to EGM2008, roughly mean sea level. Only the
India tiles are read out of the 17 GB Asia archive.

Method: every 30 m pixel is binned by elevation in 5 cm steps from -2 m to 3 m and counted per 30" cell. A
sub-district's land below a level T is the area of pixels below T; its people below T share each cell's 2020
population (GPW / GHS-POP) in proportion to the cell's land below T. T = 0 gives "below sea level today",
T = projected rise gives "below the future mean sea level". A bathtub estimate: no check that the land is
connected to the sea, no tides or surges (those are in the Aqueduct flood zone), and much of this land is
already protected by embankments, so it shows land that depends on protection, not land that will flood.
Coastal sub-districts: at least 1 km2 of land below 3 m.

Outputs:
  data/raw/sealevel/                       AR6 NetCDF files and DeltaDTM tiles (India)
  data/processed/india_adm3_sealevel.csv   one row per coastal sub-district x SSP x year
  viewer/data/slr.json                     coastal sub-districts only, TopoJSON index order
  viewer/data/lowland/                     tiles z4-10: land height above today's sea level, below 3 m
Usage: python scripts/12_sea_level.py [--workers 8]
"""
import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import geopandas as gpd
import h5py
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import remote_zip as rz  # noqa: E402

spec = importlib.util.spec_from_file_location("aq", ROOT / "scripts" / "10_aqueduct_floods.py")
aq = importlib.util.module_from_spec(spec)
spec.loader.exec_module(aq)

SLR_URL = "https://zenodo.org/api/records/6382554/files/ar6-regional-confidence.zip/content"
SLR_MEMBER = "ar6-regional-confidence/regional/confidence_output_files/medium_confidence/{s}/total_{s}_medium_confidence_values.nc"
DTM_URL = "https://data.4tu.nl/file/1da2e70f-6c4d-4b03-86bd-b53e789cc629/672eba4c-1334-44c6-8119-8879ded25912"  # Asia.zip
DTM_TILES = "https://data.4tu.nl/file/1da2e70f-6c4d-4b03-86bd-b53e789cc629/60a69899-2e67-4f9f-8761-3b57094acd12"
SSPS = ["ssp126", "ssp245", "ssp370", "ssp585"]
YEARS = [2030, 2050, 2070, 2090, 2100]
QUANT = {"lo": 0.17, "med": 0.5, "hi": 0.83}
RAW = ROOT / "data" / "raw" / "sealevel"
PROC = ROOT / "data" / "processed"
VIEW = ROOT / "viewer" / "data"
W, N, RES = aq.W, aq.N, aq.RES
NROW, NCOL = round((aq.N - aq.S) / RES), round((aq.E - aq.W) / RES)
LOW, STEP, HIGH = -2.0, 0.05, 3.0
NB = round((HIGH - LOW) / STEP)  # elevation bins; edges LOW + j*STEP, j = 0..NB
LOWLAND_BREAKS = [0, 0.5, 1, 2, 3]   # tile classes: < 0, 0-0.5, 0.5-1, 1-2, 2-3 m
LOWLAND_COLORS = ["6a1b6a", "a8327f", "d6619a", "eea2bf", "f8d8e4"]


def fetch_all(workers):
    RAW.mkdir(parents=True, exist_ok=True)
    (RAW / "deltadtm").mkdir(exist_ok=True)
    need = [s for s in SSPS if not (RAW / f"total_{s}.nc").exists()]
    if need:
        idx = rz.index(SLR_URL)
        for s in need:
            rz.extract(SLR_URL, idx[SLR_MEMBER.format(s=s)], RAW / f"total_{s}.nc")
            print(f"  AR6 {s}")
    gp = RAW / "deltadtm_tiles.gpkg"
    if not gp.exists():
        import requests
        gp.write_bytes(requests.get(DTM_TILES, timeout=300).content)
    india = aq.adm3mod.load_adm3().to_crs(4326).union_all()
    tiles = gpd.read_file(gp).cx[W:aq.E, aq.S:N]
    tiles = tiles[tiles.intersects(india)]["tile"].tolist()
    idx_path = RAW / "asia_index.json"
    idx = json.loads(idx_path.read_text()) if idx_path.exists() else rz.index(DTM_URL)
    idx_path.write_text(json.dumps(idx))
    todo = [t for t in tiles if t in idx and not (RAW / "deltadtm" / t).exists()]
    print(f"{len(tiles)} DeltaDTM tiles touch India, {len([t for t in tiles if t in idx])} in the archive, {len(todo)} to fetch")
    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(lambda t: rz.extract(DTM_URL, idx[t], RAW / "deltadtm" / t), todo))
    return sorted((RAW / "deltadtm").glob("*.tif"))


def cell_hist(paths):
    """{30" cell index: pixel counts per elevation bin (NB bins, below LOW folded into the first)}."""
    cells, counts = [], []
    for p in paths:
        with rasterio.open(p) as d:
            a = d.read(1)
            lat0, lon0 = round(d.transform.f), round(d.transform.c)
            nd = d.nodata
        r0, c0 = round((N - lat0) / RES), round((lon0 - W) / RES)
        if not (0 <= r0 < NROW and 0 <= c0 < NCOL):
            continue
        ok = (a != nd) & np.isfinite(a) & (a < HIGH)
        k = np.clip(((a - LOW) / STEP).astype(np.int32), 0, NB - 1)
        yy, xx = np.nonzero(ok)
        cid = (yy // 30) * 120 + xx // 30
        h = np.bincount(cid * NB + k[ok], minlength=14400 * NB).reshape(14400, NB)
        keep = np.nonzero(h.sum(1))[0]
        gr, gc = r0 + keep // 120, c0 + keep % 120
        inside = (gr < NROW) & (gc < NCOL)
        cells.append((gr * NCOL + gc)[inside]); counts.append(h[keep][inside])
    return np.concatenate(cells), np.concatenate(counts)


def slr_table(adm3):
    """Per sub-district: nearest AR6 grid point's sea-level rise (m) by SSP, year and quantile."""
    pts = adm3.to_crs(4326).geometry.representative_point()
    out = {}
    for s in SSPS:
        with h5py.File(RAW / f"total_{s}.nc", "r") as f:
            lat, lon, ids = f["lat"][:], f["lon"][:], f["locations"][:]
            years, q = list(f["years"][:]), f["quantiles"][:]
            grid = np.nonzero((ids > 1_000_000) & (lat > aq.S - 3) & (lat < N + 3) & (lon > W - 3) & (lon < aq.E + 3))[0]
            slc = f["sea_level_change"][:, :, grid.min():grid.max() + 1][:, :, grid - grid.min()]
        valid = (slc[:, years.index(2100), :] > -32768).all(0)
        g, glat, glon = grid[valid], lat[grid][valid], lon[grid][valid]
        coslat = np.cos(np.radians(pts.y.values))[:, None]
        near = np.argmin((pts.y.values[:, None] - glat[None]) ** 2 + ((pts.x.values[:, None] - glon[None]) * coslat) ** 2, 1)
        sl = slc[:, :, valid]
        for y in YEARS:
            for name, qq in QUANT.items():
                out[(s, y, name)] = sl[int(np.argmin(abs(q - qq))), years.index(y), near] / 1000.0
    return out


def build_tiles(paths):
    """Height above today's sea level for land below 3 m, India only, as palette tiles z4-10."""
    env, res = aq.tiles.gdal_env(), 40075016.685578488 / (256 * 2 ** 10)
    work = ROOT / "data" / "interim" / "lowland_tiles"
    work.mkdir(parents=True, exist_ok=True)
    vrt, merc, cls, rgba, xyz = (work / n for n in ("dtm.vrt", "dtm_3857.tif", "cls.tif", "rgba.vrt", "xyz"))
    lst = work / "tiles.txt"
    lst.write_text("\n".join(str(p) for p in paths))
    india = work / "india.gpkg"
    gpd.GeoDataFrame(geometry=[aq.adm3mod.load_adm3().to_crs(4326).union_all()], crs=4326).to_file(india)
    subprocess.run(["gdalbuildvrt", "-q", "-overwrite", "-input_file_list", str(lst), str(vrt)], check=True, env=env)
    subprocess.run(["gdalwarp", "-q", "-overwrite", "-t_srs", "EPSG:3857", "-tr", str(res), str(res), "-tap", "-r", "average",
                    "-cutline", str(india), "-dstnodata", "-9999", "-wm", "512", "-co", "COMPRESS=DEFLATE", "-co", "TILED=YES",
                    str(vrt), str(merc)], check=True, env=env)
    with rasterio.open(merc) as src:
        prof = src.profile
        prof.update(dtype="uint8", nodata=0, compress="deflate")
        with rasterio.open(cls, "w", **prof) as dst:
            for _, win in src.block_windows(1):
                a = src.read(1, window=win)
                c = (np.digitize(a, LOWLAND_BREAKS) + 1).astype(np.uint8)  # < 0 -> 1, 0-0.5 -> 2, ..., 2-3 -> 5
                c[(a == -9999) | (a >= LOWLAND_BREAKS[-1]) | ~np.isfinite(a)] = 0
                dst.write(c, 1, window=win)
    colors = work / "colors.txt"
    colors.write_text("nv 0 0 0 0\n" + "".join(
        f"{k} {int(c[:2], 16)} {int(c[2:4], 16)} {int(c[4:], 16)} 255\n" for k, c in enumerate(LOWLAND_COLORS, 1)))
    subprocess.run(["gdaldem", "color-relief", "-q", "-alpha", "-exact_color_entry", "-of", "VRT", str(cls), str(colors), str(rgba)],
                   check=True, env=env)
    shutil.rmtree(xyz, ignore_errors=True)
    subprocess.run([str(aq.tiles.CONDA / "python.exe"), str(aq.tiles.CONDA / "Scripts" / "gdal2tiles.py"), "--xyz", "-z", "4-10",
                    "-r", "near", "-w", "none", "-x", "-q", f"--processes={max(1, (os.cpu_count() or 2) // 2)}", str(rgba), str(xyz)],
                   check=True, env=env)
    aq.tiles.bundler.bundle(xyz, VIEW / "lowland", aq.tiles.bundler.parse_palette(",".join(LOWLAND_COLORS)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--no-tiles", action="store_true")
    args = ap.parse_args()
    paths = fetch_all(args.workers)

    adm3 = aq.adm3mod.load_adm3()
    n = len(adm3) + 1
    uid, cell, frac = aq.weights(adm3, NROW, NCOL)
    cover = np.bincount(cell, weights=frac, minlength=NROW * NCOL)
    lat = N - (np.arange(NROW) + 0.5) * RES
    cell_km2 = np.repeat(6371.0072 ** 2 * np.radians(RES) * np.abs(np.sin(np.radians(lat + RES / 2)) - np.sin(np.radians(lat - RES / 2))), NCOL)

    hc, hcount = cell_hist(paths)
    print(f"{len(hc):,} cells with land below {HIGH} m")
    # cumulative share of each cell's area below each edge (900 pixels of 1" per 30" cell)
    pos = np.full(NROW * NCOL, -1, np.int64); pos[hc] = np.arange(len(hc))
    cumshare = np.concatenate([np.zeros((len(hc), 1)), np.cumsum(hcount, 1) / 900.0], 1)  # edges 0..NB
    sel = pos[cell] >= 0
    u, c, fr = uid[sel], cell[sel], frac[sel]
    cs = cumshare[pos[c]]  # pairs x edges
    land = np.maximum(cover[c], 1e-9)
    pops = {}
    for s in aq.SOURCES:
        with rasterio.open(PROC / f"{s}2020_india_pop.tif") as d:
            p = d.read(1).astype(np.float64).ravel()
        p[p < 0] = 0
        pops[s] = p[c] * fr / land
    order_uid = adm3["uid"].values
    below = {"area": np.zeros((len(adm3), NB + 1)), **{s: np.zeros((len(adm3), NB + 1)) for s in aq.SOURCES}}
    for j in range(NB + 1):
        share = np.minimum(1.0, cs[:, j] / land)  # share of the cell's (Indian) land below this edge
        below["area"][:, j] = np.bincount(u, weights=cell_km2[c] * cs[:, j] * fr / land, minlength=n)[order_uid]
        for s in aq.SOURCES:
            below[s][:, j] = np.bincount(u, weights=pops[s] * share, minlength=n)[order_uid]

    def at(level, key):
        """Linear interpolation of the cumulative curves at a per-sub-district level (m)."""
        x = np.clip((np.asarray(level) - LOW) / STEP, 0, NB)
        j0 = np.floor(x).astype(int).clip(0, NB - 1)
        t = x - j0
        r = np.arange(len(adm3))
        return below[key][r, j0] * (1 - t) + below[key][r, j0 + 1] * t

    coastal = below["area"][:, -1] >= 1.0
    sl = slr_table(adm3)
    rows = []
    for s in SSPS:
        for y in YEARS:
            rec = {"adm3_id": adm3["adm3_id"].values, "ssp": s, "year": y,
                   **{f"slr_{k}_m": sl[(s, y, k)].round(3) for k in QUANT},
                   "area_below_today_km2": at(np.zeros(len(adm3)), "area").round(2),
                   "area_below_future_km2": at(sl[(s, y, "med")], "area").round(2)}
            for src in aq.SOURCES:
                rec[f"{src}_below_today"] = at(np.zeros(len(adm3)), src).round()
                for k in QUANT:
                    rec[f"{src}_below_{k}"] = at(sl[(s, y, k)], src).round()
            rows.append(pd.DataFrame(rec)[coastal])
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(PROC / "india_adm3_sealevel.csv", index=False)

    g = pd.read_csv(PROC / "india_adm3_gfsm.csv")
    order = list(g.loc[g["qa_flag"] != "sliver", "adm3_id"])
    pos_of = {a: i for i, a in enumerate(order)}
    cids = [a for a in adm3["adm3_id"].values[coastal] if a in pos_of]
    r_of = {a: i for i, a in enumerate(adm3["adm3_id"].values)}
    rr = np.array([r_of[a] for a in cids])
    payload = {"idx": [pos_of[a] for a in cids], "years": YEARS, "ssps": SSPS,
               "note": "slr in cm; area km2 x10; people 2020", "a0": [int(round(v * 10)) for v in at(np.zeros(len(adm3)), "area")[rr]],
               **{f"{src}0": [int(round(v)) for v in at(np.zeros(len(adm3)), src)[rr]] for src in aq.SOURCES}, "data": {}}
    for s in SSPS:
        d = {}
        for y in YEARS:
            for k in QUANT:
                d[f"{k}{y}"] = [int(round(v * 100)) for v in sl[(s, y, k)][rr]]
            d[f"a{y}"] = [int(round(v * 10)) for v in at(sl[(s, y, "med")], "area")[rr]]
            for src in aq.SOURCES:
                d[f"{src}{y}"] = [int(round(v)) for v in at(sl[(s, y, "med")], src)[rr]]
        payload["data"][s] = d
    (VIEW / "slr.json").write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")

    c = out[out["year"] == 2100].groupby("ssp")[["slr_med_m"]].median()
    t = out[out["year"].isin([2050, 2100])].groupby(["ssp", "year"])[[f"{s}_below_med" for s in aq.SOURCES] + ["area_below_future_km2"]].sum()
    print(f"{coastal.sum()} coastal sub-districts; today below sea level: "
          f"{out.loc[(out.ssp == 'ssp245') & (out.year == 2030), 'area_below_today_km2'].sum():,.0f} km2, "
          + ", ".join(f"{s} {out.loc[(out.ssp == 'ssp245') & (out.year == 2030), f'{s}_below_today'].sum() / 1e5:.1f} lakh" for s in aq.SOURCES))
    print("median coastal sub-district rise by 2100 (m):\n" + c.round(2).to_string())
    print((t / [1e5, 1e5, 1]).round(1).to_string(), "\n(people in lakh, area km2, below the median future sea level)")
    print(f"{VIEW / 'slr.json'} {(VIEW / 'slr.json').stat().st_size / 1e6:.2f} MB")
    if not args.no_tiles:
        build_tiles(paths)


if __name__ == "__main__":
    main()
