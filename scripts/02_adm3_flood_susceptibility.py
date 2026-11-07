"""Summarise GFSM v1 flood susceptibility classes per sub-district (tehsil/taluka/mandal).

Boundaries: geoBoundaries gbOpen IND ADM3 (LGD-derived, 2018, ODbL), with district/state
names attached from geoBoundaries ADM2/ADM1 via each polygon's representative point.

Method (per 1x1 deg GFSM tile, processed in parallel):
  * read the tile once, restricted to its nominal 1 deg cell so overlapping tile edges
    are not double counted;
  * rasterize the ADM3 polygons (pixel-centre rule) into an id raster in EPSG:3395;
  * accumulate per-polygon, per-class *true* area: World Mercator pixels are 30 m x 30 m
    on the map but cover 900 * cos^2(phi) / (1 - e^2 sin^2(phi)) m^2 on the ellipsoid.

Output: data/processed/india_adm3_gfsm.geojson (+ CSV, + simplified GeoJSON for web maps).

Usage: python scripts/02_adm3_flood_susceptibility.py [--workers N]
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import rasterize
from rasterio.windows import from_bounds
from pyproj import Transformer

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
GB = RAW / "boundaries" / "geoboundaries"
TILE_DIR = RAW / "gfsm" / "tiles_india"
TILE_SEL = RAW / "gfsm" / "india_tile_selection.gpkg"
OUT_DIR = ROOT / "data" / "processed"

MERC = "EPSG:3395"
E2 = 0.00669437999014  # WGS84 first eccentricity squared
CLASSES = {1: "very_low", 2: "low", 3: "moderate", 4: "high", 5: "very_high"}
NCLS = 6  # 0 = NoData, 1..5


def load_adm3():
    adm3 = gpd.read_file(GB / "geoBoundaries-IND-ADM3.geojson")
    adm3 = adm3.rename(columns={"shapeName": "adm3_name", "shapeID": "adm3_id"})
    pts = adm3[["adm3_id", "geometry"]].copy()
    pts["geometry"] = adm3.representative_point()
    for level, col in [("ADM2", "district"), ("ADM1", "state")]:
        parent = gpd.read_file(GB / f"geoBoundaries-IND-{level}.geojson")[["shapeName", "geometry"]]
        joined = gpd.sjoin(pts, parent, how="left", predicate="within").drop_duplicates("adm3_id")
        names = joined.set_index("adm3_id")["shapeName"]
        # coastlines differ between levels: fall back to the nearest parent (equal-area CRS)
        miss = names[names.isna()].index
        if len(miss):
            near = gpd.sjoin_nearest(pts[pts["adm3_id"].isin(miss)].to_crs("ESRI:102028"),
                                     parent.to_crs("ESRI:102028"), how="left")
            names.update(near.drop_duplicates("adm3_id").set_index("adm3_id")["shapeName"])
        adm3[col] = names.loc[adm3["adm3_id"]].values
    adm3["uid"] = np.arange(1, len(adm3) + 1, dtype=np.int32)  # 0 reserved for background
    return adm3[["uid", "adm3_id", "adm3_name", "district", "state", "geometry"]]


def process_tile(args):
    """Return {uid: [area_m2 for class 0..5]} for one tile."""
    tif, footprint_merc, shapes = args
    with rasterio.open(tif) as src:
        win = from_bounds(*footprint_merc.bounds, transform=src.transform)
        win = win.round_offsets().round_lengths().intersection(
            rasterio.windows.Window(0, 0, src.width, src.height))
        data = src.read(1, window=win)
        transform = src.window_transform(win)
    ids = rasterize(shapes, out_shape=data.shape, transform=transform, fill=0, dtype="int32")

    # true area per pixel row (all pixels in a row share a latitude)
    ys = transform.f + transform.e * (np.arange(data.shape[0]) + 0.5)
    lat = np.radians(Transformer.from_crs(MERC, "EPSG:4326", always_xy=True)
                     .transform(np.zeros_like(ys), ys)[1])
    row_area = abs(transform.a * transform.e) * np.cos(lat) ** 2 / (1 - E2 * np.sin(lat) ** 2)
    weights = np.broadcast_to(row_area[:, None], data.shape)

    m = ids > 0
    key = ids[m].astype(np.int64) * NCLS + np.clip(data[m], 0, NCLS - 1)
    uniq, inv = np.unique(key, return_inverse=True)
    sums = np.bincount(inv, weights=weights[m])
    out = {}
    for k, s in zip(uniq, sums):
        out.setdefault(int(k // NCLS), np.zeros(NCLS))[int(k % NCLS)] += s
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=max(1, (__import__("os").cpu_count() or 2) - 1))
    args = ap.parse_args()

    adm3 = load_adm3()
    adm3_merc = adm3.to_crs(MERC)
    tiles = gpd.read_file(TILE_SEL)
    tiles_merc = tiles.to_crs(MERC)

    jobs, missing = [], []
    for tid, fp in zip(tiles["tile_id"], tiles_merc.geometry):
        tif = TILE_DIR / f"{tid}_fsm_ei5.tif"
        if not tif.exists():
            missing.append(tid)
            continue
        hits = adm3_merc.iloc[adm3_merc.sindex.query(fp, predicate="intersects")]
        shapes = [(g.intersection(fp), uid) for g, uid in zip(hits.geometry, hits["uid"])]
        shapes = [(g, u) for g, u in shapes if not g.is_empty]
        if shapes:
            jobs.append((str(tif), fp, shapes))
    if missing:
        print(f"WARNING: {len(missing)} selected tiles not on disk: {missing[:10]}...")
    print(f"Processing {len(jobs)} tiles with {args.workers} workers")

    totals = np.zeros((len(adm3) + 1, NCLS))
    with ProcessPoolExecutor(args.workers) as ex:
        for i, res in enumerate(ex.map(process_tile, jobs, chunksize=1), 1):
            for uid, arr in res.items():
                totals[uid] += arr
            if i % 25 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)} tiles")

    km2 = totals[adm3["uid"].values] / 1e6
    valid = km2[:, 1:].sum(axis=1)
    poly_km2 = adm3.to_crs("ESRI:102028").area.values / 1e6  # Asia South Albers equal-area
    stats = pd.DataFrame({"area_km2": poly_km2.round(3),
                          "gfsm_valid_km2": valid.round(3),
                          "gfsm_nodata_km2": km2[:, 0].round(3)})
    with np.errstate(invalid="ignore", divide="ignore"):
        for c, name in CLASSES.items():
            stats[f"km2_{name}"] = km2[:, c].round(3)
        for c, name in CLASSES.items():
            stats[f"pct_{name}"] = (100 * km2[:, c] / valid).round(2)
        stats["pct_high_or_above"] = (100 * km2[:, 4:6].sum(axis=1) / valid).round(2)
        stats["mean_class"] = ((km2[:, 1:] * np.arange(1, 6)).sum(axis=1) / valid).round(3)
        stats["coverage"] = ((valid + km2[:, 0]) / poly_km2).round(3)
    # ok | sliver (degenerate polygon in source) | no_gfsm (outside GFSM coverage, e.g. Lakshadweep)
    stats["qa_flag"] = np.where(poly_km2 < 0.01, "sliver",
                                np.where(valid == 0, "no_gfsm", "ok"))

    out = pd.concat([adm3.drop(columns="uid").reset_index(drop=True), stats], axis=1)
    out = gpd.GeoDataFrame(out, geometry="geometry", crs=adm3.crs)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out.to_file(OUT_DIR / "india_adm3_gfsm.geojson", driver="GeoJSON", COORDINATE_PRECISION=6)
    out.drop(columns="geometry").to_csv(OUT_DIR / "india_adm3_gfsm.csv", index=False)
    web = out.copy()
    web["geometry"] = web.geometry.simplify(0.002, preserve_topology=True)
    web.to_file(OUT_DIR / "india_adm3_gfsm_web.geojson", driver="GeoJSON", COORDINATE_PRECISION=5)

    print(f"Wrote {len(out)} sub-districts to {OUT_DIR}")
    print(f"  qa_flag: {stats['qa_flag'].value_counts().to_dict()}")
    print(f"  missing district: {out['district'].isna().sum()}, missing state: {out['state'].isna().sum()}")


if __name__ == "__main__":
    main()
