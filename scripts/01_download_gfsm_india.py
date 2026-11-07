"""Download GFSM v1 (Zenodo record 20568218) tiles covering India.

Steps:
  1. Select 1x1 degree GFSM tiles intersecting the India boundary
     (Natural Earth 10m, India point-of-view) buffered by ~0.1 deg for coasts.
  2. Download the regional 20x20 degree ZIPs containing those tiles (md5-verified).
  3. Extract only the selected tiles to data/raw/gfsm/tiles_india/.

Usage: python scripts/01_download_gfsm_india.py [--keep-zips]
"""
import argparse
import hashlib
import time
import zipfile
from pathlib import Path

import geopandas as gpd
import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
TILE_INDEX = RAW / "gfsm" / "tile_index" / "GFSM_Tile_Index_v1.gpkg"
BOUNDARY = RAW / "boundaries" / "ne_india_pov" / "ne_10m_admin_0_countries_ind.shp"
ZIP_DIR = RAW / "gfsm" / "zips"
TILE_DIR = RAW / "gfsm" / "tiles_india"
RECORD_API = "https://zenodo.org/api/records/20568218"


def select_tiles():
    tiles = gpd.read_file(TILE_INDEX)
    countries = gpd.read_file(BOUNDARY)
    india = countries[countries["ADM0_A3"] == "IND"].to_crs(tiles.crs)
    aoi = india.union_all().buffer(0.1)
    sel = tiles[tiles.intersects(aoi)].copy()
    out = RAW / "gfsm" / "india_tile_selection.gpkg"
    sel.to_file(out, driver="GPKG")
    india.to_file(RAW / "boundaries" / "india_boundary.gpkg", driver="GPKG")
    print(f"{len(sel)} tiles selected across {sel['zip_filename'].nunique()} archives -> {out}")
    return sel


def md5sum(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while b := f.read(chunk):
            h.update(b)
    return h.hexdigest()


def download(fname, meta):
    ZIP_DIR.mkdir(parents=True, exist_ok=True)
    dest = ZIP_DIR / fname
    expected = meta["checksum"].split(":", 1)[1]
    if dest.exists() and dest.stat().st_size == meta["size"] and md5sum(dest) == expected:
        print(f"[skip] {fname} already verified")
        return dest
    part = dest.with_suffix(".zip.part")
    for attempt in range(1, 11):
        pos = part.stat().st_size if part.exists() else 0
        if pos == meta["size"]:
            break
        headers = {"Range": f"bytes={pos}-"} if pos else {}
        try:
            with requests.get(meta["links"]["self"], headers=headers, stream=True, timeout=120) as r:
                r.raise_for_status()
                mode = "ab" if r.status_code == 206 else "wb"
                with open(part, mode) as f:
                    for block in r.iter_content(1 << 20):
                        f.write(block)
            break
        except (requests.ConnectionError, requests.exceptions.ChunkedEncodingError,
                requests.Timeout) as e:
            print(f"[retry {attempt}] {fname} interrupted at {pos / 1e6:.0f} MB: {type(e).__name__}",
                  flush=True)
            time.sleep(min(60, 5 * attempt))
    else:
        raise RuntimeError(f"Giving up on {fname} after repeated connection failures")
    if md5sum(part) != expected:
        raise RuntimeError(f"Checksum mismatch for {fname}")
    part.replace(dest)
    print(f"[ok] {fname} ({meta['size'] / 1e6:.0f} MB) verified")
    return dest


def extract(zip_path, tile_ids):
    TILE_DIR.mkdir(parents=True, exist_ok=True)
    n = 0
    with zipfile.ZipFile(zip_path) as z:
        for member in z.namelist():
            name = Path(member).name
            if name.lower().endswith(".tif") and any(tid in name for tid in tile_ids):
                (TILE_DIR / name).write_bytes(z.read(member))
                n += 1
    print(f"[extract] {zip_path.name}: {n} tiles")
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep-zips", action="store_true")
    args = ap.parse_args()

    sel = select_tiles()
    files = {f["key"]: f for f in requests.get(RECORD_API, timeout=60).json()["files"]}
    total = 0
    for fname, grp in sel.groupby("zip_filename"):
        zp = download(fname, files[fname])
        total += extract(zp, set(grp["tile_id"]))
        if not args.keep_zips:
            zp.unlink()
    print(f"Done: {total} tiles in {TILE_DIR} (expected {len(sel)})")


if __name__ == "__main__":
    main()
