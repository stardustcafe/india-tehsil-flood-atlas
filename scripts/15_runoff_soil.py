"""Flood-generating runoff and monsoon soil wetness under the SSPs, per sub-district.

Source: ISIMIP3b monthly output (data.isimip.org, CC BY 4.0), the same runs as scripts/13:
  total runoff (qtot)        CWatM, WaterGAP2-2e, MIROC-INTEG-LAND x 5 CMIP6 models = 15 runs
  soil moisture (soilmoist)  CWatM (all layers summed), WaterGAP2-2e x 5 CMIP6 models = 10 runs
  (MIROC-INTEG-LAND's soil files are 2 GB each and are left out.)
historical (histsoc) and SSP1-2.6 / 3-7.0 / 5-8.5 (2015soc-from-histsoc); SSP2-4.5 is the mean of SSP1-2.6 and
SSP3-7.0, as in scripts/13.

Indicators, change against 1985-2014 for 2020-39 ... 2080-99:
  wettest-month runoff   the year's largest monthly runoff, averaged over the period: water reaching streams in
                         the peak month, a driver of river floods that does not depend on routing
  monsoon soil wetness   June-September mean soil moisture: wetter soils turn more of the next storm into runoff
Each run's sub-district value is the area-weighted mean over its 0.5 deg cells, and the change is taken per run
before the median, 17th-83rd percentile range and share of runs with an increase are reported. Runoff change is
left out where the median run's wettest month yields under 10 mm (deserts), since percentages of almost nothing
mislead.

Outputs:
  data/raw/isimip/monthly/               yearly wettest-month and JJAS means per downloaded file (npz)
  data/processed/india_adm3_runoff.csv   one row per sub-district x SSP x period
  viewer/data/runoff.json                TopoJSON order
Usage: python scripts/15_runoff_soil.py [--workers 6]
"""
import argparse
import importlib.util
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("rp", ROOT / "scripts" / "13_return_period.py")
rp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rp)
aq = rp.aq

VARS = {"qtot": ["cwatm", "watergap2-2e", "miroc-integ-land"], "soilmoist": ["cwatm", "watergap2-2e"]}
OUT = rp.RAW / "monthly"
PROC, VIEW = rp.PROC, rp.VIEW
SSPS = ["ssp126", "ssp245", "ssp370", "ssp585"]
MIN_MM = 10  # mm of runoff in the median wettest month


def file_list(var):
    out = {}
    for ghm in VARS[var]:
        url = f"{rp.SEARCH}?simulation_round=ISIMIP3b&product=OutputData&sector=water_global&variable={var}" \
              f"&time_step=monthly&model={ghm}&page_size=500"
        while url:
            d = requests.get(url, timeout=120).json()
            for r in d["results"]:
                s = r["specifiers"]
                if s["climate_forcing"] in rp.GCMS and rp.SCENS.get(s["climate_scenario"]) == s["soc_scenario"]:
                    out[(ghm, s["climate_forcing"], s["climate_scenario"])] = [f["path"] for f in r["files"]]
            url = d.get("next")
    return out


def yearly(nc):
    """Per year: largest monthly value and June-September mean (soil layers summed)."""
    with h5py.File(nc, "r") as f:
        t = f["time"]
        units = t.attrs["units"]
        units = units.decode() if isinstance(units, bytes) else units
        y0, m0 = (int(x) for x in units.split("since")[1].strip().split()[0].split("-")[:2])
        v = t[:].astype("int64")
        if units.startswith("months"):
            k = (m0 - 1) + v
            years, months = y0 + k // 12, k % 12 + 1
        else:
            days = np.datetime64(f"{y0:04d}-{m0:02d}-01", "D") + v
            years = days.astype("datetime64[Y]").astype(int) + 1970
            months = days.astype("datetime64[M]").astype(int) % 12 + 1
        lat = f["lat"][:]
        name = next(k for k in f if k in ("qtot", "soilmoist"))
        a = f[name][:]
    a = np.where(a > 1e19, np.nan, a).astype(np.float32)
    if a.ndim == 4:  # soil layers
        a = np.nansum(a, axis=1) * np.where(np.isnan(a).all(axis=1), np.nan, 1)
    if lat[0] < lat[-1]:
        a = a[:, ::-1]
    ys = np.unique(years)
    ys = ys[[(years == y).sum() == 12 for y in ys]]
    wet = np.stack([np.nanmax(a[years == y], axis=0) for y in ys])
    jjas = np.stack([np.nanmean(a[(years == y) & (months >= 6) & (months <= 9)], axis=0) for y in ys])
    return {"years": ys, "wet": wet, "jjas": jjas}


def series(paths, key, lo, hi):
    ys, a = [], []
    for p in paths:
        z = np.load(OUT / (Path(p).stem + ".npz"))
        ys.append(z["years"]); a.append(z[key])
    ys, a = np.concatenate(ys), np.concatenate(a)
    return np.nanmean(a[(ys >= lo) & (ys <= hi)], axis=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()
    files = {v: file_list(v) for v in VARS}
    allp = sorted({p for v in files for k, ps in files[v].items() for p in ps
                   if all((k[0], k[1], s) in files[v] for s in rp.SCENS)})
    print(f"{len(allp)} files, {sum(not (OUT / (Path(p).stem + '.npz')).exists() for p in allp)} to fetch", flush=True)
    with ThreadPoolExecutor(args.workers) as ex:
        list(ex.map(lambda p: rp.cutout(p, yearly, OUT / (Path(p).stem + ".npz")), allp))

    # sub-district weights on the 0.5 deg grid
    adm3 = aq.adm3mod.load_adm3()
    nrow, ncol = round((aq.N - aq.S) / aq.RES), round((aq.E - aq.W) / aq.RES)
    uid, cell, frac = aq.weights(adm3, nrow, ncol)
    k = round(rp.RES / aq.RES)
    big = (cell // ncol // k) * rp.NC + (cell % ncol) // k
    lat = aq.N - (np.arange(nrow) + 0.5) * aq.RES
    area = np.repeat(np.cos(np.radians(lat)), ncol)[cell] * frac
    n = len(adm3) + 1
    order_uid = adm3["uid"].values

    def tehsil(grid):
        v = grid.ravel()[big]
        ok = np.isfinite(v)
        s = np.bincount(uid[ok], weights=v[ok] * area[ok], minlength=n)
        w = np.bincount(uid[ok], weights=area[ok], minlength=n)
        return np.where(w > 0, s / np.maximum(w, 1e-12), np.nan)[order_uid]

    res = {}
    for var, key, tag in (("qtot", "wet", "ro"), ("soilmoist", "jjas", "sm")):
        runs = [(h, g) for h in VARS[var] for g in rp.GCMS if all((h, g, s) in files[var] for s in rp.SCENS)]
        print(f"{tag}: {len(runs)} runs", flush=True)
        ch = {s: {y: [] for y in rp.PERIODS} for s in ("ssp126", "ssp370", "ssp585")}
        bases = []
        for h, g in runs:
            base = tehsil(series(files[var][(h, g, "historical")], key, *rp.BASE))
            bases.append(base)
            for s in ch:
                for y, (lo, hi) in rp.PERIODS.items():
                    fut = tehsil(series(files[var][(h, g, s)], key, lo, hi))
                    ch[s][y].append(np.where(base > 0, 100 * (fut / np.maximum(base, 1e-12) - 1), np.nan))
        # % changes of almost no runoff mean little (deserts): drop runoff where the median wettest month < MIN_MM
        dry = (np.nanmedian(bases, 0) * 86400 * 30.4 < MIN_MM) if var == "qtot" else np.zeros(len(adm3), bool)
        print(f"  {tag}: {dry.sum()} sub-districts below {MIN_MM} mm in the wettest month left out")
        for s in ch:
            for y in rp.PERIODS:
                a = np.stack(ch[s][y])
                a[:, dry] = np.nan
                res[(tag, s, y)] = {"med": np.nanmedian(a, 0), "lo": np.nanpercentile(a, 17, 0),
                                    "hi": np.nanpercentile(a, 83, 0), "agree": np.nanmean(a > 0, 0)}
        for y in rp.PERIODS:  # no SSP2-4.5 in ISIMIP3b: mean of SSP1-2.6 and SSP3-7.0
            a, b = res[(tag, "ssp126", y)], res[(tag, "ssp370", y)]
            res[(tag, "ssp245", y)] = {f: (a[f] + b[f]) / 2 for f in a}

    rows = []
    for s in SSPS:
        for y in rp.PERIODS:
            r = {"adm3_id": adm3["adm3_id"].values, "ssp": s, "period": f"{y}-{y + 19}"}
            for tag, name in (("ro", "wettest_month_runoff"), ("sm", "monsoon_soil_moisture")):
                d = res[(tag, s, y)]
                r.update({f"{name}_change_pct": d["med"].round(1), f"{name}_lo": d["lo"].round(1),
                          f"{name}_hi": d["hi"].round(1), f"{name}_share_up": d["agree"].round(2)})
            rows.append(pd.DataFrame(r))
    out = pd.concat(rows)
    out.to_csv(PROC / "india_adm3_runoff.csv", index=False)

    g = pd.read_csv(PROC / "india_adm3_gfsm.csv")
    pos = pd.Series(np.arange(len(adm3)), index=adm3["adm3_id"].values)[g.loc[g["qa_flag"] != "sliver", "adm3_id"]].values
    def arr(v, scale=1):
        x = v[pos]
        return [None if not np.isfinite(z) else int(round(z * scale)) for z in x]
    data = {}
    for s in SSPS:
        dd = {}
        for y in rp.PERIODS:
            ro, sm = res[("ro", s, y)], res[("sm", s, y)]
            dd.update({f"ro{y}": arr(ro["med"]), f"roLo{y}": arr(ro["lo"]), f"roHi{y}": arr(ro["hi"]),
                       f"roA{y}": arr(ro["agree"], 100), f"sm{y}": arr(sm["med"], 10), f"smA{y}": arr(sm["agree"], 100)})
        data[s] = dd
    payload = {"periods": list(rp.PERIODS), "base": list(rp.BASE),
               "note": "ro: % change in wettest-month runoff; sm: % change x10 in June-September soil moisture; A: % of runs up",
               "data": data}
    (VIEW / "runoff.json").write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")

    summ = out.groupby(["ssp", "period"])[["wettest_month_runoff_change_pct", "monsoon_soil_moisture_change_pct"]].median()
    print(summ.round(1).to_string(), "\n(median sub-district, % change against 1985-2014)")
    print(f"{VIEW / 'runoff.json'} {(VIEW / 'runoff.json').stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
