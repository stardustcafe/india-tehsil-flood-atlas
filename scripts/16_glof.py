"""Glacial lake outburst floods (GLOF): dangerous glacial lakes upstream of each sub-district, today and as glaciers
retreat.

Source: Zheng, G. et al. (2021), Increasing risk of glacial lake outburst floods from future Third Pole deglaciation,
Nature Climate Change 11, 411-417; data https://doi.org/10.5281/zenodo.4477945 (CC BY 4.0). About 7,000 lakes
mapped in 2015 with a hazard score (0-1) from lake volume, dam type, slope and the chance of ice or rock falling
in, and lakes expected to form in glacier-bed hollows by 2050 and 2100 under RCP2.6 / 4.5 / 8.5 as glaciers shrink.
High and Very High hazard use the paper's present-day class edges (0.571 and 0.721). A future scenario counts
today's lakes with today's class plus the new lakes reaching those edges (the projected files re-score existing
lakes about 0.05 lower, relative to the larger set, so they are not mixed with today's scores).

Downstream reach: from each lake's 0.5 deg cell along the DDM30 river network (scripts/13) for up to 250 km of
flow path; outburst floods have caused damage this far downstream (Teesta, October 2023). A sub-district is
reached if one of those cells covers at least 10% of its area. A screening link at 50 km cells, not a flood model:
it says which sub-districts lie down-valley from dangerous lakes, not how deep a flood would be.

Outputs:
  data/processed/india_adm3_glof.csv   one row per sub-district reached by a High+ lake in any scenario
  viewer/data/glof.json                those sub-districts, TopoJSON index order
  viewer/data/glof_lakes.geojson       High+ lakes that drain to India, with their class per scenario
Usage: python scripts/16_glof.py
"""
import importlib.util
import json
import zipfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("rp", ROOT / "scripts" / "13_return_period.py")
rp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rp)
aq = rp.aq

URL = "https://zenodo.org/api/records/4477945/files/2.%20hazard%20and%20risk%20assessments.zip/content"
RAW = ROOT / "data" / "raw" / "glof"
PROC, VIEW = rp.PROC, rp.VIEW
H, VH = 0.57125, 0.72075   # present-day class edges, High and Very High
REACH_KM = 250
SCEN = {"now": "present_hazard_and_risk_2015"} | {f"{r}_{y}": f"projected_hazard_and_risk_{y}_{r}"
                                                    for r in ("rcp26", "rcp45", "rcp85") for y in (2050, 2100)}


def lakes():
    d = RAW / "hazard"
    if not (d / f"{SCEN['now']}.shp").exists():
        RAW.mkdir(parents=True, exist_ok=True)
        z = RAW / "hazard.zip"
        z.write_bytes(requests.get(URL, timeout=600).content)
        with zipfile.ZipFile(z) as zf:
            zf.extractall(d)
    out = {}
    for k, name in SCEN.items():
        g = gpd.read_file(d / f"{name}.shp")
        if k != "now":  # projected files re-score existing lakes on a shifted scale: keep today's, add the new ones
            g = g[g["TIME"] == "F"]
        df = pd.DataFrame({"id": g["GLID"], "lat": g.geometry.y, "lon": g.geometry.x,
                           "hazard": g["HAZARD" if k == "now" else "HRGI"].astype(float),
                           "volume": g["VOLUME"].astype(float), "new": k != "now"})
        out[k] = df if k == "now" else pd.concat([out["now"], df], ignore_index=True)
    return out


def downstream_cells(t, lat, lon, la, lo):
    """(India-window cell, km of flow path from the lake) along the path below a point, up to REACH_KM."""
    nr, nc = len(lat), len(lon)
    r, c = int(np.argmin(np.abs(lat - la))), int(np.argmin(np.abs(lon - lo)))
    i, dist, cells = r * nc + c, 0.0, []
    while i >= 0 and dist <= REACH_KM:
        rr, cc = divmod(i, nc)
        wr = round((rp.BBOX[3] - rp.RES / 2 - lat[rr]) / rp.RES)
        wc = round((lon[cc] - rp.BBOX[0] - rp.RES / 2) / rp.RES)
        if 0 <= wr < rp.NR and 0 <= wc < rp.NC:
            cells.append((wr * rp.NC + wc, dist))
        j = t[i]
        if j < 0:
            break
        jr, jc = divmod(j, nc)
        dy = (lat[jr] - lat[rr]) * 111.2
        dx = (lon[jc] - lon[cc]) * 111.2 * np.cos(np.radians(lat[rr]))
        dist += float(np.hypot(dx, dy))
        i = j
    return cells


def main():
    L = lakes()
    t, lat, lon, _ = rp.ddm30()

    adm3 = aq.adm3mod.load_adm3()
    nrow, ncol = round((aq.N - aq.S) / aq.RES), round((aq.E - aq.W) / aq.RES)
    uid, cell, frac = aq.weights(adm3, nrow, ncol)
    k = round(rp.RES / aq.RES)
    big = (cell // ncol // k) * rp.NC + (cell % ncol) // k
    pair = pd.DataFrame({"uid": uid, "big": big, "a": frac}).groupby(["uid", "big"], as_index=False)["a"].sum()
    pair["share"] = pair["a"] / pair.groupby("uid")["a"].transform("sum")
    pair = pair[pair["share"] >= 0.1]
    india_cells = set(pair["big"])
    by_cell = pair.groupby("big")["uid"].apply(list).to_dict()

    reach, keep = {}, {}
    for key, df in L.items():
        df = df[df["hazard"] >= H]
        rows = []
        for r in df.itertuples():
            if r.id not in reach:
                reach[r.id] = [(c, d) for c, d in downstream_cells(t, lat, lon, r.lat, r.lon) if c in india_cells]
            if reach[r.id]:
                near = {}
                for c, d in reach[r.id]:
                    for u in by_cell[c]:
                        near[u] = min(near.get(u, 1e9), d)
                for u, d in near.items():
                    rows.append((u, r.id, r.hazard >= VH, r.volume, r.new, d))
                keep.setdefault(r.id, {"lat": r.lat, "lon": r.lon})[key] = 2 if r.hazard >= VH else 1
        keep_rows = pd.DataFrame(rows, columns=["uid", "id", "vh", "volume", "new", "km"])
        L[key] = keep_rows
        print(f"{key}: {keep_rows['id'].nunique()} High+ lakes drain to India within {REACH_KM} km, "
              f"{keep_rows['uid'].nunique()} sub-districts reached")

    out = pd.DataFrame({"uid": adm3["uid"].values, "adm3_id": adm3["adm3_id"].values})
    for key, df in L.items():
        g = df.groupby("uid")
        out[f"{key}_lakes"] = out["uid"].map(g["id"].nunique()).fillna(0).astype(int)
        out[f"{key}_vh"] = out["uid"].map(df[df["vh"]].groupby("uid")["id"].nunique()).fillna(0).astype(int)
        out[f"{key}_new"] = out["uid"].map(df[df["new"]].groupby("uid")["id"].nunique()).fillna(0).astype(int)
        out[f"{key}_volume_mcm"] = out["uid"].map(g["volume"].sum() / 1e6).fillna(0).round(1)
        out[f"{key}_nearest_km"] = out["uid"].map(g["km"].min()).fillna(-1).round()
    out = out[out.filter(like="_lakes").sum(axis=1) > 0].drop(columns="uid")
    out.to_csv(PROC / "india_adm3_glof.csv", index=False)

    g = pd.read_csv(PROC / "india_adm3_gfsm.csv")
    order = list(g.loc[g["qa_flag"] != "sliver", "adm3_id"])
    pos = {a: i for i, a in enumerate(order)}
    out = out[out["adm3_id"].isin(pos)].assign(i=lambda d: d["adm3_id"].map(pos)).sort_values("i")
    payload = {"note": "High+ hazard glacial lakes within 250 km upstream; v: volume in million m3; d: km to the nearest",
               "scen": list(SCEN), "lakes": {k: int(L[k]["id"].nunique()) for k in SCEN},
               "vhLakes": {k: int(L[k].loc[L[k]["vh"], "id"].nunique()) for k in SCEN},
               "idx": out["i"].tolist(),
               **{f"{k}{f}": [int(round(v)) for v in out[f"{k}_{c}"]] for k in SCEN
                  for f, c in (("n", "lakes"), ("vh", "vh"), ("new", "new"), ("v", "volume_mcm"), ("d", "nearest_km"))}}
    (VIEW / "glof.json").write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
    feats = [{"type": "Feature", "geometry": {"type": "Point", "coordinates": [round(v["lon"], 4), round(v["lat"], 4)]},
              "properties": {s: v.get(s, 0) for s in SCEN}} for v in keep.values()]
    (VIEW / "glof_lakes.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats},
                                                        separators=(",", ":")), encoding="utf-8")
    print(f"{len(out)} sub-districts; {len(feats)} lakes on the map")
    print(out[[c for c in out.columns if c.endswith("_lakes")]].astype(bool).sum().to_string())
    d = out["now_nearest_km"]
    print(f"today, nearest High+ lake: within 50 km {(d.between(0, 50)).sum()}, 50-150 km {(d.between(50.01, 150)).sum()}, "
          f"150-250 km {(d > 150).sum()} sub-districts")
    for f in ("glof.json", "glof_lakes.geojson"):
        print(f"{VIEW / f} {(VIEW / f).stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
