"""Urban growth on flood-prone land: built-up area inside the 1-in-100-year flood zone, 2000, 2020 and 2030.

Source: GHS-BUILT-S R2023A (Pesaresi & Politis 2023, European Commission JRC, CC BY 4.0), built-up surface at
30 arc-seconds from Landsat and Sentinel-2, for the 2000 and 2020 epochs and the 2030 projection (extrapolated
from 1975-2020 trends). The 4326 / 30" grid matches the analysis grid, so cells are cropped, not resampled.

Flood zone: today's 1-in-100-year river or coastal flood zone from WRI Aqueduct (scripts/10), the same cells as
"people in the flood zone". Built-up share of land is the urban flash-flood proxy: sealed surfaces send rain
straight to drains, which the river and coastal maps do not cover.

Outputs:
  data/raw/ghsl/india/                  India crops of the three epochs
  data/processed/india_adm3_urban.csv   one row per sub-district
  viewer/data/urban.json                TopoJSON order
Usage: python scripts/17_urban_landuse.py
"""
import importlib.util
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import requests
from rasterio.transform import from_origin
from rasterio.windows import Window

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("aq", ROOT / "scripts" / "10_aqueduct_floods.py")
aq = importlib.util.module_from_spec(spec)
spec.loader.exec_module(aq)

BASE = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_S_GLOBE_R2023A"
EPOCHS = [2000, 2020, 2030]
RAW = ROOT / "data" / "raw" / "ghsl"
CROP = RAW / "india"
PROC = ROOT / "data" / "processed"
VIEW = ROOT / "viewer" / "data"
W, N, RES = aq.W, aq.N, aq.RES
NROW, NCOL = round((aq.N - aq.S) / RES), round((aq.E - aq.W) / RES)


def crop(epoch):
    out = CROP / f"built_{epoch}.tif"
    if out.exists():
        return out
    name = f"GHS_BUILT_S_E{epoch}_GLOBE_R2023A_4326_30ss_V1_0"
    z = RAW / f"{name}.zip"
    if not z.exists():
        with requests.get(f"{BASE}/{name[:-5]}/V1-0/{name}.zip", stream=True, timeout=900) as r, open(z, "wb") as f:
            r.raise_for_status()
            for ch in r.iter_content(1 << 20):
                f.write(ch)
    tif = RAW / f"{name}.tif"
    if not tif.exists():
        with zipfile.ZipFile(z) as zf:
            zf.extract(f"{name}.tif", RAW)
    with rasterio.open(tif) as d:
        # GHSL's 30" lattice is offset from ours by a fraction of a cell: take the nearest cell
        r0, c0 = round((d.transform.f - N) / RES), round((W - d.transform.c) / RES)
        a = d.read(1, window=Window(c0, r0, NCOL, NROW)).astype(np.float32)
        a[(a == d.nodata) | ~np.isfinite(a) | (a < 0)] = 0
    CROP.mkdir(parents=True, exist_ok=True)
    prof = {"driver": "GTiff", "dtype": "float32", "count": 1, "width": NCOL, "height": NROW, "crs": "EPSG:4326",
            "transform": from_origin(W, N, RES, RES), "compress": "deflate", "predictor": 3, "tiled": True}
    with rasterio.open(out, "w", **prof) as dst:
        dst.write(a, 1)
    tif.unlink()
    print(f"  {epoch}: {a.sum() / 1e6:,.0f} km2 built up in the India window")
    return out


def main():
    built = {e: rasterio.open(crop(e)).read(1).astype(np.float64).ravel() / 1e6 for e in EPOCHS}  # km2 per cell
    zone = ((aq.depth(aq.river_name("historical", None, None)) > 0)
            | (aq.depth(aq.coast_name("historical", None)) > 0)).ravel()

    adm3 = aq.adm3mod.load_adm3()
    n = len(adm3) + 1
    uid, cell, frac = aq.weights(adm3, NROW, NCOL)
    cover = np.bincount(cell, weights=frac, minlength=NROW * NCOL)
    share = frac / np.maximum(cover[cell], 1e-9)
    lat = N - (np.arange(NROW) + 0.5) * RES
    cell_km2 = np.repeat(6371.0072 ** 2 * np.radians(RES) * np.abs(np.sin(np.radians(lat + RES / 2)) - np.sin(np.radians(lat - RES / 2))), NCOL)
    order_uid = adm3["uid"].values
    tsum = lambda v: np.bincount(uid, weights=v, minlength=n)[order_uid]

    land = tsum(cell_km2[cell] * frac)
    df = pd.DataFrame({"adm3_id": adm3["adm3_id"].values, "land_km2": land.round(2)})
    for e in EPOCHS:
        b = built[e][cell] * share
        df[f"built_{e}_km2"] = tsum(b).round(3)
        df[f"built_{e}_floodzone_km2"] = tsum(b * zone[cell]).round(3)
    new = (df["built_2020_km2"] - df["built_2000_km2"]).clip(lower=0)
    new_z = (df["built_2020_floodzone_km2"] - df["built_2000_floodzone_km2"]).clip(lower=0)
    df["new_2000_2020_km2"], df["new_2000_2020_floodzone_km2"] = new.round(3), new_z.round(3)
    df["new_in_floodzone_pct"] = np.where(new > 0.05, 100 * new_z / new.where(new > 0, 1), np.nan).round(1)
    df["built_share_2020_pct"] = (100 * df["built_2020_km2"] / df["land_km2"].where(df["land_km2"] > 0)).round(2)
    df.to_csv(PROC / "india_adm3_urban.csv", index=False)

    g = pd.read_csv(PROC / "india_adm3_gfsm.csv")
    d = df.set_index("adm3_id").loc[g.loc[g["qa_flag"] != "sliver", "adm3_id"]]
    x100 = lambda c: [int(round(v * 100)) for v in d[c].fillna(0)]
    payload = {"note": "areas in km2 x100; s: built-up share of land x100 (%)", "epochs": EPOCHS,
               **{f"b{e}": x100(f"built_{e}_km2") for e in EPOCHS}, **{f"z{e}": x100(f"built_{e}_floodzone_km2") for e in EPOCHS},
               "s": x100("built_share_2020_pct")}
    (VIEW / "urban.json").write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")

    tot = df.sum(numeric_only=True)
    print(f"India built-up: 2000 {tot.built_2000_km2:,.0f} km2, 2020 {tot.built_2020_km2:,.0f}, 2030 {tot.built_2030_km2:,.0f}")
    print(f"  in today's flood zone: 2000 {tot.built_2000_floodzone_km2:,.0f}, 2020 {tot.built_2020_floodzone_km2:,.0f}, "
          f"2030 {tot.built_2030_floodzone_km2:,.0f} km2")
    print(f"  of land built 2000-2020, {100 * tot.new_2000_2020_floodzone_km2 / tot.new_2000_2020_km2:.1f}% is in the flood zone "
          f"(flood zone covers {100 * tsum(cell_km2[cell] * frac * zone[cell]).sum() / land.sum():.1f}% of land)")
    print(f"{VIEW / 'urban.json'} {(VIEW / 'urban.json').stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
