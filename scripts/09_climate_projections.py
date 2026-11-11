"""Projected change in extreme rainfall per sub-district under CMIP6 scenarios, and a flood-pressure index.

Source: World Bank Climate Change Knowledge Portal (CCKP), CMIP6 0.25 deg collection (CC BY 4.0),
multi-model ensemble climatologies (median, 10th and 90th percentile across models), reference 1995-2014.
  rx5day  annual maximum 5-day precipitation (mm)   - river and widespread flooding
  rx1day  annual maximum 1-day precipitation (mm)   - flash and urban flooding
  pr      mean precipitation                         - wetness overall
Scenarios SSP1-2.6, SSP2-4.5, SSP3-7.0, SSP5-8.5 (CMIP6 successors of RCP2.6/4.5/7.0/8.5);
periods 2020-2039, 2040-2059, 2060-2079, 2080-2099.

Robustness (IPCC AR6 Atlas method, shipped with the anomaly files): robust = >=66% of models show a
significant change and >=80% agree on its sign; otherwise "no change" or "conflicting signal". Over India
this strict test passes mostly late in the century under high emissions, so the viewer uses a simpler
agreement measure: at least 90% of models project an increase (the ensemble 10th percentile is above zero).

Per sub-district (area-weighted over 0.25 deg cells, from a 0.025 deg rasterisation; sub-districts smaller
than one fine cell use the cell at their representative point):
  base_rx5day_mm, d_rx5day_pct (median), d_rx5day_p10_pct, d_rx5day_p90_pct, robust_rx5day_pct (% of area
  with a robust change), d_rx1day_pct, d_pr_pct
Flood-pressure index = pct_high_or_above x (1 + d_rx5day_pct / 100): today's share of High + Very High land,
scaled by how much heavier the 5-day rainfall extreme becomes. A relative screening indicator, not a flood
probability: GFSM is static and runoff, river flow, dams and land-use change are not modelled.

Outputs:
  data/processed/india_adm3_climate.csv   long table, one row per sub-district x scenario x period
  viewer/data/climate.json               compact arrays + 0.25 deg India grids for the map layer
Usage: python scripts/09_climate_projections.py
"""
import base64
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import requests
from rasterio.features import rasterize
from rasterio.transform import from_origin
from rasterio.windows import from_bounds

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("adm3", ROOT / "scripts" / "02_adm3_flood_susceptibility.py")
adm3mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adm3mod)

BASE_URL = "https://wbg-cckp.s3.amazonaws.com/data/cmip6-x0.25"
RAW = ROOT / "data" / "raw" / "climate" / "cckp"
OUT_CSV = ROOT / "data" / "processed" / "india_adm3_climate.csv"
OUT_VIEW = ROOT / "viewer" / "data" / "climate.json"
SCENARIOS = {"ssp126": "SSP1-2.6", "ssp245": "SSP2-4.5", "ssp370": "SSP3-7.0", "ssp585": "SSP5-8.5"}
PERIODS = ["2020-2039", "2040-2059", "2060-2079", "2080-2099"]
REF = "1995-2014"
# India window on the 0.25 deg grid (cell centres at multiples of 0.25)
W, S, E, N = 67.875, 5.875, 98.125, 38.125
RES, SUB = 0.25, 10  # fine grid 0.025 deg for area weights


def fname(var, kind, ens, stat, period):
    return f"{kind}-{var}-annual-mean_cmip6-x0.25_ensemble-all-{ens}_climatology_{stat}_{period}.nc"


def fetch(var, kind, ens, stat, period):
    name = fname(var, kind, ens, stat, period)
    path = RAW / name
    if not path.exists():
        url = f"{BASE_URL}/{var}/ensemble-all-{ens}/{name}"
        for attempt in range(5):
            try:
                r = requests.get(url, timeout=300)
                r.raise_for_status()
                path.write_bytes(r.content)
                break
            except requests.RequestException as e:
                if attempt == 4:
                    raise SystemExit(f"{url}: {e}")
    return path


def read(path, layer):
    """India window of one NetCDF variable, north-up, NaN for missing."""
    with rasterio.open(f'netcdf:"{path}":{layer}') as d:
        win = from_bounds(W, S, E, N, transform=d.transform).round_offsets().round_lengths()
        a = d.read(1, window=win).astype(np.float64)
        t = d.window_transform(win)
        if t.e > 0:  # stored south-up
            a = a[::-1]
        a[(a > 1e19) | (a < -1e19)] = np.nan
    return a


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    nrow, ncol = round((N - S) / RES), round((E - W) / RES)
    transform = from_origin(W, N, RES, RES)

    adm3 = adm3mod.load_adm3()
    g = pd.read_csv(ROOT / "data" / "processed" / "india_adm3_gfsm.csv")  # row order = viewer TopoJSON order
    # fine rasterisation -> per sub-district weights on each 0.25 deg cell (weight ~ ground area)
    fine_t = from_origin(W, N, RES / SUB, RES / SUB)
    ids = rasterize(zip(adm3.geometry, adm3["uid"]), out_shape=(nrow * SUB, ncol * SUB), transform=fine_t,
                    fill=0, dtype="int32")
    lat_f = N - (np.arange(nrow * SUB) + 0.5) * RES / SUB
    wrow = np.cos(np.radians(lat_f))
    cell = (np.arange(nrow * SUB)[:, None] // SUB) * ncol + (np.arange(ncol * SUB)[None, :] // SUB)
    m = ids > 0
    key = ids[m].astype(np.int64) * (nrow * ncol) + cell[m]
    uk, inv = np.unique(key, return_inverse=True)
    wts = np.bincount(inv, weights=np.broadcast_to(wrow[:, None], ids.shape)[m])
    k_uid, k_cell = uk // (nrow * ncol), uk % (nrow * ncol)
    india_cell = np.zeros(nrow * ncol, bool)
    india_cell[k_cell] = True
    # sub-districts too small for the fine grid: representative point
    missing = set(adm3["uid"]) - set(np.unique(k_uid))
    if missing:
        rp = adm3[adm3["uid"].isin(missing)]
        pts = rp.geometry.representative_point()
        r = np.clip(((N - pts.y) / RES).astype(int), 0, nrow - 1)
        c = np.clip(((pts.x - W) / RES).astype(int), 0, ncol - 1)
        k_uid = np.concatenate([k_uid, rp["uid"].values])
        k_cell = np.concatenate([k_cell, (r * ncol + c).values])
        wts = np.concatenate([wts, np.ones(len(rp))])
        print(f"{len(missing)} sub-districts smaller than the 0.025 deg grid use their representative point")
    n = len(adm3) + 1

    def zonal(grid):
        v = grid.ravel()[k_cell]
        ok = ~np.isnan(v)
        num = np.bincount(k_uid[ok], weights=(v * wts)[ok], minlength=n)
        den = np.bincount(k_uid[ok], weights=wts[ok], minlength=n)
        with np.errstate(invalid="ignore"):
            return (num / den)[adm3["uid"].values]

    base = {v: read(fetch(v, "climatology", "historical", "median", REF), f"climatology-{v}-annual-mean")
            for v in ("rx5day", "rx1day", "pr")}
    assert base["rx5day"].shape == (nrow, ncol), base["rx5day"].shape
    zb = {v: zonal(a) for v, a in base.items()}
    # sanity: wettest and driest places
    for name, lon, lat in [("Mawsynram", 91.58, 25.30), ("Jaisalmer", 70.91, 26.92)]:
        print(f"  baseline Rx5day {name}: {base['rx5day'][int((N - lat) / RES), int((lon - W) / RES)]:.0f} mm")

    rows, view = [], {}
    h = g.set_index("adm3_id").reindex(adm3["adm3_id"])["pct_high_or_above"].values
    for ssp in SCENARIOS:
        for per in PERIODS:
            print(f"{ssp} {per}")
            f5 = {st: fetch("rx5day", "anomaly", ssp, st, per) for st in ("median", "p10", "p90")}
            d5 = {st: read(p, "anomaly-rx5day-annual-mean") for st, p in f5.items()}
            sig = read(f5["median"], "anomalysignificance-rx5day-annual-mean")
            robust = np.where(np.isnan(sig), 1.0, 0.0)  # robust is stored as missing
            robust[np.isnan(d5["median"])] = np.nan
            d1 = read(fetch("rx1day", "anomaly", ssp, "median", per), "anomaly-rx1day-annual-mean")
            dp = read(fetch("pr", "anomaly", ssp, "median", per), "anomaly-pr-annual-mean")

            pct = lambda anom, var: 100 * zonal(anom) / zb[var]
            t = pd.DataFrame({
                "adm3_id": adm3["adm3_id"].values, "scenario": SCENARIOS[ssp], "period": per,
                "base_rx5day_mm": zb["rx5day"].round(1),
                "d_rx5day_pct": pct(d5["median"], "rx5day").round(1),
                "d_rx5day_p10_pct": pct(d5["p10"], "rx5day").round(1),
                "d_rx5day_p90_pct": pct(d5["p90"], "rx5day").round(1),
                "robust_rx5day_pct": (100 * zonal(robust)).round(0),
                "d_rx1day_pct": pct(d1, "rx1day").round(1),
                "d_pr_pct": pct(dp, "pr").round(1),
            })
            t["flood_pressure"] = (h * (1 + t["d_rx5day_pct"] / 100)).round(1)
            rows.append(t)

            # 0.25 deg grid for the map: Rx5day change %, robust flag; outside India = -128
            gp = 100 * d5["median"] / base["rx5day"]
            q = np.clip(np.round(gp), -127, 127)
            q[~india_cell.reshape(nrow, ncol) | np.isnan(q)] = -128
            agree = (d5["p10"] > 0) & india_cell.reshape(nrow, ncol)  # >= 90% of models project an increase
            view[f"{ssp}_{per[:4]}"] = {
                "d5": t["d_rx5day_pct"].tolist(), "lo": t["d_rx5day_p10_pct"].tolist(), "hi": t["d_rx5day_p90_pct"].tolist(),
                "d1": t["d_rx1day_pct"].tolist(), "dp": t["d_pr_pct"].tolist(),
                "grid": base64.b64encode(q.astype(np.int8).tobytes()).decode(),
                "gagree": base64.b64encode(np.packbits(agree.ravel())).decode(),
            }

    out = pd.concat(rows, ignore_index=True)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)

    # viewer: arrays in TopoJSON feature order (scripts/03 keeps the non-sliver rows of india_adm3_gfsm)
    order = g.loc[g["qa_flag"] != "sliver", "adm3_id"]
    def clean(a):
        v = pd.Series(np.asarray(a, float), index=adm3["adm3_id"].values).reindex(order).values
        return [None if not np.isfinite(x) else round(float(x), 1) for x in v]
    payload = {"scenarios": SCENARIOS, "periods": PERIODS, "ref": REF,
               "grid": {"west": W, "north": N, "res": RES, "rows": nrow, "cols": ncol},
               "base5": clean(zb["rx5day"]),
               "data": {k: {**{f: clean(v[f]) for f in ("d5", "lo", "hi", "d1", "dp")},
                            "grid": v["grid"], "gagree": v["gagree"]} for k, v in view.items()}}
    OUT_VIEW.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")

    s = out[out["period"] == "2080-2099"].groupby("scenario")[["d_rx5day_pct", "d_rx1day_pct", "d_pr_pct"]].median()
    print("Median sub-district change by 2080-2099 (%):\n" + s.round(1).to_string())
    print(f"{OUT_VIEW} {OUT_VIEW.stat().st_size / 1e6:.2f} MB")



if __name__ == "__main__":
    main()
