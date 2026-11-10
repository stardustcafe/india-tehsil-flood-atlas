"""Population exposed to each GFSM susceptibility class, per sub-district and on a 30" grid.

Population sources (30 arc-second, ~1 km, 2020), chosen with --source:
  gpw  GPW v4.11 population count adjusted to 2015 UN WPP country totals (NASA SEDAC, CC BY 4.0);
       spreads census counts uniformly within each census unit.
  ghs  GHS-POP R2023A epoch 2020, EPSG:4326 30ss (EC JRC, CC BY 4.0); distributes census counts
       by built-up volume, so it places people within sub-districts.
Any other 30" population GeoTIFF on the global grid works via --pop with a new --source name.

Method (dasymetric split of each population cell by GFSM land area):
  * every 1x1 deg GFSM tile covers exactly 120x120 population cells, so each 30 m pixel is
    assigned to (sub-district, population cell, GFSM class) and its true ellipsoidal area summed;
  * a cell's population is shared out in proportion to the area of each (sub-district, class)
    inside it, over the cell's total GFSM land area (classes 1-5, including the part outside India,
    so border cells are not over-allocated); GFSM NoData pixels (water / masked) get no people;
  * population in cells with no GFSM land is reported as unallocated.

Outputs (data/processed/):
  india_adm3_pop_{source}.csv              adm3_id + population by class
  india_adm3_flood_exposure.geojson/.csv   GFSM statistics joined with every available source,
                                           columns prefixed {source}_ (e.g. ghs_pop_high_or_above)
  {source}2020_india_pop.tif, {source}2020_india_pop_high.tif
                                           30" grids: population, population in High+Very High

Usage: python scripts/06_population_exposure.py --source {gpw,ghs} [--pop PATH] [--workers N]
"""
import argparse
import importlib.util
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import rasterize
from rasterio.transform import from_origin
from rasterio.windows import Window, from_bounds
from pyproj import Transformer

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("adm3", ROOT / "scripts" / "02_adm3_flood_susceptibility.py")
adm3mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adm3mod)

POP_ROOT = ROOT / "data" / "raw" / "population"
OUT_DIR = ROOT / "data" / "processed"
MERC, E2, NCLS = adm3mod.MERC, adm3mod.E2, adm3mod.NCLS
CLASSES = adm3mod.CLASSES
CPD = 120  # 30 arc-second cells per degree
# India window on the global 30" grid (lon 68..98, lat 6..38)
LON0, LAT1, NCOL, NROW = 68, 38, 30 * CPD, 32 * CPD


def process_tile(args):
    """Return (keys, area_m2) with key = (uid * NROW*NCOL + cell) * NCLS + class, uid 0 = outside India."""
    tif, footprint_merc, shapes = args
    with rasterio.open(tif) as src:
        win = from_bounds(*footprint_merc.bounds, transform=src.transform)
        win = win.round_offsets().round_lengths().intersection(Window(0, 0, src.width, src.height))
        data = src.read(1, window=win)
        transform = src.window_transform(win)
    ids = rasterize(shapes, out_shape=data.shape, transform=transform, fill=0, dtype="int32")

    to_geo = Transformer.from_crs(MERC, "EPSG:4326", always_xy=True)
    ys = transform.f + transform.e * (np.arange(data.shape[0]) + 0.5)
    xs = transform.c + transform.a * (np.arange(data.shape[1]) + 0.5)
    lon, _ = to_geo.transform(xs, np.zeros_like(xs))
    _, lat = to_geo.transform(np.zeros_like(ys), ys)
    row = np.clip(np.floor((LAT1 - lat) * CPD).astype(np.int64), 0, NROW - 1)
    col = np.clip(np.floor((np.asarray(lon) - LON0) * CPD).astype(np.int64), 0, NCOL - 1)
    latr = np.radians(lat)
    row_area = abs(transform.a * transform.e) * np.cos(latr) ** 2 / (1 - E2 * np.sin(latr) ** 2)

    # strips of rows keep per-worker memory small (a full tile is ~14 M pixels)
    keys, sums = [], []
    for r0 in range(0, data.shape[0], 256):
        sl = slice(r0, r0 + 256)
        cell = row[sl, None] * NCOL + col[None, :]
        key = (ids[sl].astype(np.int64) * (NROW * NCOL) + cell) * NCLS + np.clip(data[sl], 0, NCLS - 1)
        uniq, inv = np.unique(key.ravel(), return_inverse=True)
        w = np.repeat(row_area[sl], data.shape[1])
        keys.append(uniq); sums.append(np.bincount(inv.ravel(), weights=w))
    uniq, inv = np.unique(np.concatenate(keys), return_inverse=True)
    return uniq, np.bincount(inv, weights=np.concatenate(sums))


def find_pop(path, source):
    if path:
        return Path(path)
    d = POP_ROOT / source
    hits = sorted(d.rglob("*.tif"))
    if not hits:
        raise SystemExit(f"No population GeoTIFF under {d}. Download the 2020 30 arc-second "
                         "population count (see README) and unzip it there.")
    return hits[0]


def read_pop_window(path):
    with rasterio.open(path) as src:
        res = src.res[0]
        if abs(res - 1 / CPD) > 1e-9:
            raise SystemExit(f"{path.name}: expected 30 arc-second cells, got {res}")
        c0f, r0f = (LON0 - src.bounds.left) * CPD, (src.bounds.top - LAT1) * CPD
        c0, r0 = round(c0f), round(r0f)
        if max(abs(c0f - c0), abs(r0f - r0)) > 0.01:
            print(f"  warning: grid offset by ({c0f - c0:.3f}, {r0f - r0:.3f}) cells from the 30\" lattice")
        pop = src.read(1, window=Window(c0, r0, NCOL, NROW), boundless=True, fill_value=0).astype(np.float64)
        if src.nodata is not None:
            pop[pop == src.nodata] = 0
    pop[~np.isfinite(pop) | (pop < 0)] = 0
    return pop.ravel()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default="gpw", help="short name: gpw or ghs (sets input folder and output names)")
    ap.add_argument("--pop", help="30 arc-second population GeoTIFF (default: first .tif under data/raw/population/<source>)")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    args = ap.parse_args()

    pop_path = find_pop(args.pop, args.source)
    pop = read_pop_window(pop_path)
    print(f"Population raster: {pop_path.name}, {pop.sum():,.0f} people in the India window")

    adm3 = adm3mod.load_adm3()
    adm3_merc = adm3.to_crs(MERC)
    tiles = gpd.read_file(adm3mod.TILE_SEL)
    jobs = []
    for tid, fp in zip(tiles["tile_id"], tiles.to_crs(MERC).geometry):
        tif = adm3mod.TILE_DIR / f"{tid}_fsm_ei5.tif"
        if not tif.exists():
            continue
        hits = adm3_merc.iloc[adm3_merc.sindex.query(fp, predicate="intersects")]
        shapes = [(g.intersection(fp), uid) for g, uid in zip(hits.geometry, hits["uid"])]
        shapes = [(g, u) for g, u in shapes if not g.is_empty]
        if shapes:  # population cells never cross tiles, so tiles without India need no processing
            jobs.append((str(tif), fp, shapes))
    print(f"Processing {len(jobs)} tiles with {args.workers} workers")

    keys, areas = [], []
    with ProcessPoolExecutor(args.workers) as ex:
        for i, (k, a) in enumerate(ex.map(process_tile, jobs, chunksize=1), 1):
            keys.append(k); areas.append(a)
            if i % 50 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)} tiles")
    keys, areas = np.concatenate(keys), np.concatenate(areas)
    cls = keys % NCLS
    cell = (keys // NCLS) % (NROW * NCOL)
    uid = keys // NCLS // (NROW * NCOL)

    land = cls > 0
    denom = np.bincount(cell[land], weights=areas[land], minlength=NROW * NCOL)  # GFSM land per cell, any side of border
    share = np.zeros_like(areas)
    ok = land & (denom[cell] > 0)
    share[ok] = areas[ok] / denom[cell[ok]]
    people = pop[cell] * share

    n = len(adm3) + 1
    by_cls = np.zeros((n, NCLS))
    np.add.at(by_cls, (uid, cls), people)
    allocated = people[uid > 0].sum()
    india_cells = np.unique(cell[uid > 0])
    unalloc = pop[india_cells][denom[india_cells] == 0].sum()
    print(f"Allocated {allocated:,.0f} people to sub-districts; {unalloc:,.0f} in cells with no GFSM land")

    # sub-district table
    P = by_cls[adm3["uid"].values]
    tot = P[:, 1:].sum(axis=1)
    t = pd.DataFrame({"adm3_id": adm3["adm3_id"].values, "pop_total": tot.round()})
    for c, name in CLASSES.items():
        t[f"pop_{name}"] = P[:, c].round()
    t["pop_high_or_above"] = P[:, 4:6].sum(axis=1).round()
    with np.errstate(invalid="ignore", divide="ignore"):
        t["pct_pop_high_or_above"] = (100 * P[:, 4:6].sum(axis=1) / tot).round(2)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t.to_csv(OUT_DIR / f"india_adm3_pop_{args.source}.csv", index=False)

    joined = gpd.read_file(OUT_DIR / "india_adm3_gfsm.geojson")
    for f in sorted(OUT_DIR.glob("india_adm3_pop_*.csv")):
        src = f.stem.removeprefix("india_adm3_pop_")
        s = pd.read_csv(f).set_index("adm3_id").add_prefix(f"{src}_")
        joined = joined.merge(s, left_on="adm3_id", right_index=True, how="left")
    joined.to_file(OUT_DIR / "india_adm3_flood_exposure.geojson", driver="GeoJSON", COORDINATE_PRECISION=6)
    joined.drop(columns="geometry").to_csv(OUT_DIR / "india_adm3_flood_exposure.csv", index=False)

    # 30" grids for mapping: population and population in High + Very High, India cells only
    in_india = np.zeros(NROW * NCOL, bool)
    in_india[india_cells] = True
    high = np.bincount(cell[(uid > 0) & (cls >= 4)], weights=people[(uid > 0) & (cls >= 4)], minlength=NROW * NCOL)
    total = np.bincount(cell[uid > 0], weights=people[uid > 0], minlength=NROW * NCOL)
    profile = dict(driver="GTiff", width=NCOL, height=NROW, count=1, dtype="float32", crs="EPSG:4326",
                   transform=from_origin(LON0, LAT1, 1 / CPD, 1 / CPD), nodata=-1, compress="deflate", tiled=True)
    for name, grid in [(f"{args.source}2020_india_pop.tif", total), (f"{args.source}2020_india_pop_high.tif", high)]:
        g = np.where(in_india, grid, -1).astype(np.float32).reshape(NROW, NCOL)
        with rasterio.open(OUT_DIR / name, "w", **profile) as dst:
            dst.write(g, 1)

    print(f"India: {tot.sum():,.0f} people, {t['pop_high_or_above'].sum():,.0f} "
          f"({100 * t['pop_high_or_above'].sum() / tot.sum():.1f}%) in High or Very High areas")


if __name__ == "__main__":
    main()
