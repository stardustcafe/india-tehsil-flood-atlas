"""Build compact TopoJSON for the web viewer from data/processed/india_adm3_gfsm.geojson.

Keeps short property keys, then uses mapshaper (via npx) for topology-preserving
simplification (no gaps/overlaps between neighbours) and adds a dissolved state layer.

Output: viewer/data/india_adm3.topo.json  with objects `adm3` and `states`.

Usage: python scripts/03_build_viewer_data.py [--keep 6%]
"""
import argparse
import json
import subprocess
from pathlib import Path

import geopandas as gpd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "processed" / "india_adm3_gfsm.geojson"
OUT = ROOT / "viewer" / "data" / "india_adm3.topo.json"
TMP = ROOT / "data" / "interim" / "viewer_adm3.geojson"

KEYS = {"adm3_name": "n", "district": "d", "state": "s", "area_km2": "a",
        "pct_very_low": "p1", "pct_low": "p2", "pct_moderate": "p3", "pct_high": "p4",
        "pct_very_high": "p5", "pct_high_or_above": "h", "mean_class": "m"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep", default="6%", help="mapshaper -simplify retention")
    args = ap.parse_args()

    g = gpd.read_file(SRC)
    g = g[g["qa_flag"] != "sliver"].reset_index(drop=True)
    out = g[list(KEYS)].rename(columns=KEYS)
    for k in ["a", "p1", "p2", "p3", "p4", "p5", "h"]:
        out[k] = out[k].round(1)
    out["m"] = out["m"].round(2)
    out["hk"] = (g["km2_high"] + g["km2_very_high"]).round(1)  # high+very high area, km2
    out["va"] = g["gfsm_valid_km2"].round(1)  # area with GFSM data, km2
    out["q"] = (g["qa_flag"] == "no_gfsm").astype(int)
    out["i"] = range(len(out))
    out = gpd.GeoDataFrame(out, geometry=g.geometry, crs=g.crs)

    TMP.parent.mkdir(parents=True, exist_ok=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_file(TMP, driver="GeoJSON")

    cmd = ["npx", "-y", "mapshaper", "-i", str(TMP), "name=adm3",
           "-simplify", args.keep, "keep-shapes", "planar",
           "-clean",
           "-dissolve", "s", "+", "name=states",
           "-o", str(OUT), "format=topojson", "quantization=100000", "target=adm3,states"]
    subprocess.run(cmd, check=True, shell=True)

    topo = json.loads(OUT.read_text(encoding="utf-8"))
    n = len(topo["objects"]["adm3"]["geometries"])
    print(f"{OUT} : {OUT.stat().st_size / 1e6:.2f} MB, {n} sub-districts, "
          f"{len(topo['objects']['states']['geometries'])} states")


if __name__ == "__main__":
    main()
