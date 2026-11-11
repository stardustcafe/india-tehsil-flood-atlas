"""How often today's 1-in-100-year river flood would come under the SSPs, per river cell and sub-district.

Source: ISIMIP3b daily river discharge (data.isimip.org, CC BY 4.0), 0.5 deg, from three global hydrological
models (CWatM, WaterGAP2-2e, MIROC-INTEG-LAND) each driven by five bias-adjusted CMIP6 models (GFDL-ESM4,
IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0, UKESM1-0-LL): 15 runs, historical (histsoc) and SSP1-2.6 / 3-7.0 / 5-8.5
(2015soc-from-histsoc). ISIMIP3b has no SSP2-4.5 runs. The ISIMIP files API masks each file to India and cuts
the India box on the server, so ~20 MB per decade is downloaded instead of ~1 GB; only annual maxima are kept.

Method (after Hirabayashi et al. 2013): for each run and river cell, a Gumbel distribution (L-moments) is fitted
to the annual maximum daily discharge of 1985-2014 to give today's 100-year flow, and to the annual maxima of each
future 20-year period. The future return period of today's 100-year flow is 1 / P(exceed) under the future fit,
capped at 1-1000 years. Each run is compared with its own historical run, so model bias largely cancels.
Reported: median of the 15 runs, the 17th-83rd percentile range, and the share of runs in which the flood comes
more often. SSP2-4.5 is shown as the geometric mean of SSP1-2.6 and SSP3-7.0 and labelled as such.

River cells: 0.5 deg cells draining at least 10,000 km2, from the DDM30 flow directions that the models route
on (Doell & Lehner 2002, ISIMIP). Each sub-district takes the river cell with the largest 100-year flow among the
cells covering at least 10% of its area; others get no value.

Outputs:
  data/raw/isimip/amax/                  annual maxima per downloaded file (npz)
  data/processed/india_adm3_returnperiod.csv  one row per sub-district with a river cell x SSP x period
  viewer/data/rp.json                    sub-district values, TopoJSON order
  viewer/data/rivers.geojson             river cells with their values, for the map layer
Usage: python scripts/13_return_period.py [--workers 8]
"""
import argparse
import importlib.util
import json
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("aq", ROOT / "scripts" / "10_aqueduct_floods.py")
aq = importlib.util.module_from_spec(spec)
spec.loader.exec_module(aq)

API = "https://files.isimip.org/api/v2"
SEARCH = "https://data.isimip.org/api/v1/datasets/"
GHMS = ["cwatm", "watergap2-2e", "miroc-integ-land"]
GCMS = ["gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll"]
SCENS = {"historical": "histsoc", "ssp126": "2015soc-from-histsoc", "ssp370": "2015soc-from-histsoc",
         "ssp585": "2015soc-from-histsoc"}
BASE = (1985, 2014)
PERIODS = {2020: (2020, 2039), 2040: (2040, 2059), 2060: (2060, 2079), 2080: (2080, 2099)}
BBOX = [68, 98, 6, 38]  # W E S N, the same window as scripts/10-12
RES = 0.5
NR, NC = 64, 60
RIVER_KM2 = 10000  # upstream area of a river cell
DDM30 = "ISIMIP3a/InputData/geo_conditions/river_routing/ddm30_flowdir_cru_neva.nc"
# DDM30 codes 1-8 run clockwise from east
MOVES = {1: (0, 1), 2: (-1, 1), 3: (-1, 0), 4: (-1, -1), 5: (0, -1), 6: (1, -1), 7: (1, 0), 8: (1, 1)}
RAW = ROOT / "data" / "raw" / "isimip"
AMAX = RAW / "amax"
PROC = ROOT / "data" / "processed"
VIEW = ROOT / "viewer" / "data"


def file_list():
    """(ghm, gcm, scen) -> ISIMIP paths of the daily discharge files overlapping the years needed."""
    out = {}
    for ghm in GHMS:
        url = f"{SEARCH}?simulation_round=ISIMIP3b&product=OutputData&sector=water_global&variable=dis" \
              f"&time_step=daily&model={ghm}&page_size=500"
        while url:
            d = requests.get(url, timeout=120).json()
            for r in d["results"]:
                s = r["specifiers"]
                if s["climate_forcing"] not in GCMS or SCENS.get(s["climate_scenario"]) != s["soc_scenario"]:
                    continue
                lo, hi = BASE if s["climate_scenario"] == "historical" else (2020, 2099)
                paths = [f["path"] for f in r["files"]
                         if int(f["path"][-12:-8]) <= hi and int(f["path"][-7:-3]) >= lo]
                out[(ghm, s["climate_forcing"], s["climate_scenario"])] = sorted(paths)
            url = d.get("next")
    return out


def _range(url, a, b):
    for attempt in range(8):
        try:
            r = requests.get(url, headers={"Range": f"bytes={a}-{b}"}, timeout=300)
            if r.status_code == 206 and len(r.content) == b - a + 1:
                return r.content
        except requests.RequestException:
            pass
        time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"{url}: bytes {a}-{b} failed")


def download(url, dest, conns=4):
    """The server throttles each connection, so fetch 2 MB ranges over several at once."""
    size = int(requests.head(url, timeout=120).headers["content-length"])
    step = 2 << 20
    with ThreadPoolExecutor(conns) as ex, open(dest, "wb") as f:
        for ch in ex.map(lambda a: _range(url, a, min(a + step, size) - 1), range(0, size, step)):
            f.write(ch)


def annual_max(nc):
    with h5py.File(nc, "r") as f:
        t = f["time"]
        units = t.attrs["units"].decode() if isinstance(t.attrs["units"], bytes) else t.attrs["units"]
        start = np.datetime64(units.split("since")[1].strip().split()[0].replace("-1-1", "-01-01"), "D")
        days = start + t[:].astype("int64")
        years = days.astype("datetime64[Y]").astype(int) + 1970
        lat = f["lat"][:]
        d = f["dis"][:]
    d = np.where(d > 1e19, np.nan, d)
    if lat[0] < lat[-1]:
        d = d[:, ::-1]
    ys = np.unique(years)
    return ys, np.stack([np.nanmax(d[years == y], axis=0) for y in ys]).astype(np.float32)


def fetch(path):
    out = AMAX / (Path(path).stem + ".npz")
    if out.exists():
        return out
    for attempt in range(4):
        try:
            j = requests.post(API, json={"paths": [path], "operations": [
                {"operation": "mask_country", "country": "ind"},
                {"operation": "select_bbox", "bbox": BBOX}]}, timeout=120).json()
            while j.get("status") not in ("finished", "failed"):
                time.sleep(10)
                j = requests.get(j["job_url"], timeout=120).json()
            if j["status"] == "failed" or j["meta"].get("errors"):
                raise RuntimeError(str(j.get("meta")))
            z = RAW / "tmp" / (Path(path).stem + ".zip")
            download(j["file_url"], z)
            with zipfile.ZipFile(z) as zf:
                name = next(n for n in zf.namelist() if n.endswith(".nc"))
                nc = zf.extract(name, RAW / "tmp")
            ys, a = annual_max(nc)
            np.savez_compressed(out, years=ys, amax=a)
            Path(nc).unlink(); z.unlink()
            print(f"  {Path(path).stem[:60]} {ys[0]}-{ys[-1]}", flush=True)
            return out
        except Exception as e:  # server queue hiccups: start the job again
            print(f"  retry {attempt + 1} {Path(path).name}: {e}", flush=True)
            time.sleep(30)
    raise RuntimeError(f"failed: {path}")


def upstream_area():
    """km2 draining through each 0.5 deg cell of the India window, accumulated over the global DDM30 network."""
    nc = RAW / "geo" / Path(DDM30).name
    if not nc.exists():
        nc.parent.mkdir(parents=True, exist_ok=True)
        nc.write_bytes(requests.get("https://files.isimip.org/" + DDM30, timeout=300).content)
    with h5py.File(nc, "r") as f:
        d, lat, lon = f["flowdirection"][:], f["lat"][:], f["lon"][:]  # lat ascending
    nr, ncol = d.shape
    acc = ((111.2 * RES) ** 2 * np.cos(np.radians(lat)))[:, None] * (d < 1e19)
    tgt = np.full(d.shape, -1)
    for k, (dr, dc) in MOVES.items():
        rr, cc = np.nonzero(d == k)
        tr, tc = rr + dr, (cc + dc) % ncol
        ok = (tr >= 0) & (tr < nr)
        tgt[rr[ok], cc[ok]] = tr[ok] * ncol + tc[ok]
    t, acc = tgt.ravel(), acc.ravel().copy()
    indeg = np.bincount(t[t >= 0], minlength=t.size)
    todo = list(np.nonzero(indeg == 0)[0])
    while todo:  # topological order: each cell passes its total on once everything upstream has arrived
        i = todo.pop()
        j = t[i]
        if j >= 0:
            acc[j] += acc[i]
            indeg[j] -= 1
            if indeg[j] == 0:
                todo.append(j)
    acc = acc.reshape(d.shape)[::-1]  # north first, like the discharge files
    r0, c0 = round((lat[-1] - (BBOX[3] - RES / 2)) / RES), round((BBOX[0] + RES / 2 - lon[0]) / RES)
    return acc[r0:r0 + NR, c0:c0 + NC]


def series(paths, lo, hi):
    ys, a = [], []
    for p in paths:
        z = np.load(AMAX / (Path(p).stem + ".npz"))
        ys.append(z["years"]); a.append(z["amax"])
    ys, a = np.concatenate(ys), np.concatenate(a)
    m = (ys >= lo) & (ys <= hi)
    _, first = np.unique(ys[m], return_index=True)
    return a[m][first]


def gumbel(x):
    """L-moment Gumbel fit along axis 0: (location, scale)."""
    x = np.sort(x, axis=0)
    n = x.shape[0]
    b0 = x.mean(0)
    b1 = (np.arange(n)[:, None, None] / (n - 1) * x).mean(0)
    l2 = 2 * b1 - b0
    scale = np.maximum(l2 / np.log(2), 1e-6)
    return b0 - 0.5772 * scale, scale


def return_period(q, loc, scale):
    p = 1 - np.exp(-np.exp(-(q - loc) / scale))
    return np.clip(1 / np.maximum(p, 1e-3), 1, 1000)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    AMAX.mkdir(parents=True, exist_ok=True)
    (RAW / "tmp").mkdir(exist_ok=True)
    files = file_list()
    runs = [(h, g) for h in GHMS for g in GCMS if all((h, g, s) in files for s in SCENS)]
    allp = sorted({p for k, v in files.items() if k[:2] in runs for p in v})
    print(f"{len(runs)} runs, {len(allp)} files, {sum(not (AMAX / (Path(p).stem + '.npz')).exists() for p in allp)} to fetch", flush=True)
    with ThreadPoolExecutor(args.workers) as ex:
        list(ex.map(fetch, allp))

    # per run: return period of its own historical 100-year flow in each future period
    q100, mean_peak, rp = [], [], {s: {y: [] for y in PERIODS} for s in SCENS if s != "historical"}
    for h, g in runs:
        base = series(files[(h, g, "historical")], *BASE)
        loc, sc = gumbel(base)
        q = loc - sc * np.log(-np.log(1 - 1 / 100))
        q100.append(q); mean_peak.append(base.mean(0))
        for s in rp:
            for y, (lo, hi) in PERIODS.items():
                rp[s][y].append(return_period(q, *gumbel(series(files[(h, g, s)], lo, hi))))
    q100, mean_peak = np.nanmedian(q100, 0), np.nanmedian(mean_peak, 0)
    river = np.isfinite(q100) & (upstream_area() >= RIVER_KM2)
    stats = {}
    for s in rp:
        for y in PERIODS:
            a = np.log(np.stack(rp[s][y]))
            stats[(s, y)] = {"med": np.exp(np.nanmedian(a, 0)), "lo": np.exp(np.nanpercentile(a, 17, 0)),
                             "hi": np.exp(np.nanpercentile(a, 83, 0)), "agree": np.nanmean(a < np.log(100), 0)}
    for y in PERIODS:  # SSP2-4.5 is not in ISIMIP3b: geometric mean of SSP1-2.6 and SSP3-7.0
        a, b = stats[("ssp126", y)], stats[("ssp370", y)]
        stats[("ssp245", y)] = {k: (np.sqrt(a[k] * b[k]) if k != "agree" else (a[k] + b[k]) / 2) for k in a}
    ssps = ["ssp126", "ssp245", "ssp370", "ssp585"]

    # sub-district -> its largest river cell
    adm3 = aq.adm3mod.load_adm3()
    nrow, ncol = round((aq.N - aq.S) / aq.RES), round((aq.E - aq.W) / aq.RES)
    uid, cell, frac = aq.weights(adm3, nrow, ncol)
    k = round(RES / aq.RES)
    big = (cell // ncol // k) * NC + (cell % ncol) // k
    pair = pd.DataFrame({"uid": uid, "big": big, "a": frac}).groupby(["uid", "big"], as_index=False)["a"].sum()
    pair["share"] = pair["a"] / pair.groupby("uid")["a"].transform("sum")
    pair = pair[(pair["share"] >= 0.1) & river.ravel()[pair["big"]]]
    pair["q100"] = q100.ravel()[pair["big"]]
    best = pair.sort_values("q100").groupby("uid").tail(1).set_index("uid")["big"]
    rows = []
    for s in ssps:
        for y in PERIODS:
            st = stats[(s, y)]
            c = best.values
            rows.append(pd.DataFrame({"uid": best.index, "ssp": s, "period": f"{y}-{y + 19}",
                                      "q100_today": q100.ravel()[c].round(), "rp_median": st["med"].ravel()[c].round(1),
                                      "rp_lo": st["lo"].ravel()[c].round(1), "rp_hi": st["hi"].ravel()[c].round(1),
                                      "share_more_often": st["agree"].ravel()[c].round(2)}))
    out = pd.concat(rows).merge(adm3[["uid", "adm3_id"]], on="uid").drop(columns="uid")
    out.insert(0, "adm3_id", out.pop("adm3_id"))
    out.to_csv(PROC / "india_adm3_returnperiod.csv", index=False)

    g = pd.read_csv(PROC / "india_adm3_gfsm.csv")
    order = adm3.set_index("adm3_id").loc[g.loc[g["qa_flag"] != "sliver", "adm3_id"], "uid"]
    has = order.isin(best.index).values
    cidx = best.reindex(order.values).fillna(0).astype(int).values
    def arr(v, scale=10):
        x = v.ravel()[cidx]
        return [int(round(x[i] * scale)) if has[i] and np.isfinite(x[i]) else 0 for i in range(len(x))]
    payload = {"periods": list(PERIODS), "base": list(BASE), "note": "return period x10 (years), agree x100, q m3/s",
               "q": [int(round(v)) if h else 0 for v, h in zip(q100.ravel()[cidx], has)],
               "data": {s: {f"{f}{y}": arr(stats[(s, y)][f], 100 if f == "agree" else 10)
                            for y in PERIODS for f in ("med", "lo", "hi", "agree")} for s in ssps}}
    (VIEW / "rp.json").write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")

    feats = []
    for r, c in zip(*np.nonzero(river)):
        props = {"q": int(round(q100[r, c]))}
        for s in ssps:
            for y in PERIODS:
                props[f"{s}_{y}"] = round(float(stats[(s, y)]["med"][r, c]), 1)
        x, y0 = BBOX[0] + (c + 0.5) * RES, BBOX[3] - (r + 0.5) * RES
        feats.append({"type": "Feature", "properties": props, "geometry": {"type": "Point", "coordinates": [x, y0]}})
    (VIEW / "rivers.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats},
                                                    separators=(",", ":")), encoding="utf-8")

    print(f"{river.sum()} river cells, {has.sum()} sub-districts with a river cell")
    summ = out.groupby(["ssp", "period"]).agg(median_rp=("rp_median", "median"),
                                              share_more_often=("rp_median", lambda v: (v < 100).mean()))
    print(summ.round(2).to_string(), "\n(sub-districts with a river cell: median return period of today's 100-year flood)")
    for f in ("rp.json", "rivers.geojson"):
        print(f"{VIEW / f} {(VIEW / f).stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
