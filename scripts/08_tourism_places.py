"""Notable places to visit in India from OpenStreetMap, categorised, ranked and tagged with GFSM flood class.

Sources:
  * OpenStreetMap via the Overpass API (ODbL): tourism, historic, religious, nature and leisure features
    that carry a `wikidata` tag. The Wikidata link is the notability filter; OSM has no visitor counts.
  * Wikidata (CC0): number of Wikipedia/Wikimedia sitelinks per place, used as a prominence score
    (how many language editions describe the place), which drives what the map shows first.
    OSM `wikidata` tags are sometimes wrong (a memorial linked to the person, a viewpoint to the peak it
    faces), so a link only counts when the item has coordinates within MATCH_KM of the OSM feature and is not
    a person; features linked to settlements or mountains are dropped. UNESCO status: Wikidata P1435 = Q9259
    or the OSM whc tags.
Each place gets its sub-district (spatial join), the GFSM class of the 30 m pixel at its point, and the share
of High + Very High land within RADIUS_M of it (GFSM NoData, e.g. water, excluded from the share).

Outputs:
  data/processed/india_places.csv, india_places.geojson
  viewer/data/places.geojson   compact properties for the viewer
Raw API responses are cached in data/raw/places/ (delete to refresh).

Usage: python scripts/08_tourism_places.py
"""
import importlib.util
import json
import math
import re
import time
from collections import defaultdict
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import requests
from pyproj import Transformer
from rasterio.windows import Window

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("adm3", ROOT / "scripts" / "02_adm3_flood_susceptibility.py")
adm3mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adm3mod)

RAW = ROOT / "data" / "raw" / "places"
OUT_DIR = ROOT / "data" / "processed"
VIEW = ROOT / "viewer" / "data" / "places.geojson"
UA = {"User-Agent": "india-tehsil-flood-atlas/0.1 (research; https://github.com/stardustcafe/india-tehsil-flood-atlas)"}
OVERPASS = "https://overpass-api.de/api/interpreter"
RADIUS_M = 250
MATCH_KM = 10
# Wikidata P31 classes: settlements and mountains are dropped (a "tourism" tag on a town or a peak);
# people keep the place but give no prominence
DROP_P31 = {"Q515", "Q1549591", "Q3957", "Q532", "Q5119", "Q8502", "Q207326", "Q131541", "Q1149652", "Q1187811",
            "Q486972", "Q1637706", "Q15284"}
PERSON_P31 = {"Q5"}

QUERY = """[out:json][timeout:900];
(
  nwr["tourism"~"^(attraction|museum|gallery|viewpoint|theme_park|zoo|aquarium)$"]["wikidata"]({bbox});
  nwr["historic"~"^(monument|memorial|fort|castle|palace|archaeological_site|ruins|tomb|city_gate|citywalls|temple|monastery|church|mosque|stepwell)$"]["wikidata"]({bbox});
  nwr["amenity"="place_of_worship"]["wikidata"]({bbox});
  nwr["natural"~"^(beach|waterfall|cave_entrance|hot_spring)$"]["wikidata"]({bbox});
  nwr["waterway"="waterfall"]["wikidata"]({bbox});
  nwr["leisure"~"^(park|garden|nature_reserve|water_park)$"]["wikidata"]({bbox});
  nwr["boundary"~"^(national_park|protected_area)$"]["wikidata"]({bbox});
);
out center tags;"""

# category key -> label; order is the assignment priority when a feature has several matching tags
CATEGORIES = {
    "heritage": "Forts, palaces and monuments",
    "religious": "Temples, mosques, churches and shrines",
    "museum": "Museums and galleries",
    "beach": "Beaches",
    "waterfall": "Waterfalls",
    "wildlife": "National parks, sanctuaries and zoos",
    "park": "Parks and gardens",
    "viewpoint": "Viewpoints and caves",
    "amusement": "Theme and water parks",
    "attraction": "Other attractions",
}
RELIGIOUS_HISTORIC = {"temple", "monastery", "church", "mosque"}
PARK_NAME = re.compile(r"national park|sanctuary|reserve|wls|np|biosphere|tiger|bird|wildlife|forest", re.I)


def categorise(t):
    h, tour, nat, lei, b = t.get("historic"), t.get("tourism"), t.get("natural"), t.get("leisure"), t.get("boundary")
    if t.get("amenity") == "place_of_worship" or h in RELIGIOUS_HISTORIC:
        return "religious", t.get("religion") or h or "place_of_worship"
    if h:
        return "heritage", h
    if tour in ("museum", "gallery"):
        return "museum", tour
    if nat in ("peak", "volcano", "glacier"):
        return None, None
    if nat == "beach":
        return "beach", nat
    if nat == "waterfall" or t.get("waterway") == "waterfall":
        return "waterfall", "waterfall"
    if b == "protected_area" and lei != "nature_reserve" and not PARK_NAME.search(t.get("name:en") or t.get("name") or ""):
        return None, None  # generic protected areas (heritage zones, restricted islands, sites under legal protection)
    if b in ("national_park", "protected_area") or lei == "nature_reserve" or tour in ("zoo", "aquarium"):
        return "wildlife", b or lei or tour
    if lei in ("park", "garden"):
        return "park", lei
    if tour == "viewpoint" or nat in ("cave_entrance", "hot_spring"):
        return "viewpoint", tour or nat
    if tour == "theme_park" or lei == "water_park":
        return "amusement", tour or lei
    if tour == "attraction":
        return "attraction", "attraction"
    return None, None


def overpass(bbox, cache):
    if cache.exists():
        return json.loads(cache.read_text(encoding="utf-8"))
    q = QUERY.replace("{bbox}", ",".join(map(str, bbox)))
    for attempt in range(6):
        try:
            r = requests.post(OVERPASS, data={"data": q}, headers=UA, timeout=1000)
            if r.status_code == 200:
                cache.write_text(r.text, encoding="utf-8")
                return r.json()
            print(f"  overpass {r.status_code}, retrying")
        except requests.RequestException as e:
            print(f"  overpass {type(e).__name__}, retrying")
        time.sleep(30 * (attempt + 1))
    raise SystemExit(f"Overpass failed for {bbox}")


def wikidata_items(qids):
    cache = RAW / "wikidata_items.json"
    known = json.loads(cache.read_text(encoding="utf-8")) if cache.exists() else {}
    todo = [q for q in qids if q not in known]
    for i in range(0, len(todo), 50):
        batch = todo[i:i + 50]
        for attempt in range(5):
            try:
                r = requests.get("https://www.wikidata.org/w/api.php", headers=UA, timeout=60, params={
                    "action": "wbgetentities", "ids": "|".join(batch), "props": "sitelinks|labels|claims",
                    "languages": "en", "format": "json"})
                r.raise_for_status()
                break
            except requests.RequestException:
                time.sleep(5 * (attempt + 1))
        else:
            raise SystemExit("Wikidata API failed")
        for q, e in r.json().get("entities", {}).items():
            cl = e.get("claims", {})
            ids = lambda prop: [c["mainsnak"].get("datavalue", {}).get("value", {}).get("id") for c in cl.get(prop, [])]
            coord = next((c["mainsnak"]["datavalue"]["value"] for c in cl.get("P625", []) if "datavalue" in c["mainsnak"]), None)
            known[q] = {"n": len(e.get("sitelinks", {})), "en": e.get("labels", {}).get("en", {}).get("value"),
                        "p31": [i for i in ids("P31") if i], "whc": "Q9259" in ids("P1435"),
                        "xy": [coord["longitude"], coord["latitude"]] if coord else None}
        for q in batch:  # invalid or deleted ids
            known.setdefault(q, {"n": 0, "en": None, "p31": [], "whc": False, "xy": None})
        if (i // 50) % 40 == 0:
            print(f"  wikidata {i + len(batch)}/{len(todo)}")
            cache.write_text(json.dumps(known), encoding="utf-8")
    cache.write_text(json.dumps(known), encoding="utf-8")
    return known


def sample_gfsm(df):
    """GFSM class at each point and share of High+ within RADIUS_M, reading only small windows."""
    to_merc = Transformer.from_crs("EPSG:4326", adm3mod.MERC, always_xy=True)
    x, y = to_merc.transform(df["lon"].values, df["lat"].values)
    cls = np.zeros(len(df), np.uint8)
    share = np.full(len(df), np.nan)
    by_tile = defaultdict(list)
    for k, (lon, lat) in enumerate(zip(df["lon"], df["lat"])):
        la, lo = math.floor(lat), math.floor(lon)
        by_tile[f"{'N' if la >= 0 else 'S'}{abs(la):02d}E{lo:03d}"].append(k)
    for tid, ks in by_tile.items():
        tif = adm3mod.TILE_DIR / f"{tid}_fsm_ei5.tif"
        if not tif.exists():
            continue
        with rasterio.open(tif) as src:
            for k in ks:
                # Mercator pixels shrink on the ground by cos(lat); widen the window to keep RADIUS_M true
                rpx = int(math.ceil(RADIUS_M / (abs(src.res[0]) * math.cos(math.radians(df["lat"].iat[k])))))
                row, col = src.index(x[k], y[k])
                win = Window(col - rpx, row - rpx, 2 * rpx + 1, 2 * rpx + 1)
                a = src.read(1, window=win, boundless=True, fill_value=0)
                yy, xx = np.ogrid[-rpx:rpx + 1, -rpx:rpx + 1]
                disk = a[(yy ** 2 + xx ** 2) <= rpx ** 2]
                cls[k] = a[rpx, rpx]
                valid = disk[disk > 0]
                if valid.size:
                    share[k] = 100 * (valid >= 4).mean()
    return cls, share


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    rows = []
    for i, (s, w) in enumerate([(lat, lon) for lat in (6, 22) for lon in (68, 83)]):
        bbox = (s, w, s + 16, w + 15)  # south, west, north, east
        print(f"Overpass chunk {i + 1}/4 {bbox}")
        for e in overpass(bbox, RAW / f"overpass_{s}_{w}.json")["elements"]:
            t = e.get("tags", {})
            cat, sub = categorise(t)
            name = t.get("name:en") or t.get("name")
            lat, lon = (e["lat"], e["lon"]) if e["type"] == "node" else (e.get("center", {}).get("lat"), e.get("center", {}).get("lon"))
            if not cat or lat is None:
                continue
            rows.append({"osm": f"{e['type'][0]}{e['id']}", "wikidata": t["wikidata"].split(";")[0].strip(),
                         "name": name, "name_local": t.get("name"), "category": cat, "subtype": sub,
                         "unesco": int(t.get("heritage:operator", "").lower() == "whc" or "whc:inscription_date" in t
                                       or "ref:whc" in t),
                         "lat": lat, "lon": lon})
    df = pd.DataFrame(rows)
    prio = {c: i for i, c in enumerate(CATEGORIES)}
    # one row per Wikidata item: OSM often maps a site as both a node and its outline
    df = (df.assign(_p=df["category"].map(prio), _n=df["osm"].str[0].map({"n": 0, "w": 1, "r": 2}))
            .sort_values(["wikidata", "_p", "_n"]).drop_duplicates("wikidata").drop(columns=["_p", "_n"]))
    print(f"{len(df)} places with Wikidata links in the India window")

    wd = wikidata_items(sorted(df["wikidata"].unique()))
    item = df["wikidata"].map(lambda q: wd.get(q, {}))
    p31 = item.map(lambda it: set(it.get("p31", [])))
    drop = p31.map(lambda c: bool(c & DROP_P31))
    print(f"Dropping {drop.sum()} features linked to settlements or mountains")

    def km(it, lon, lat):
        if not it.get("xy"):
            return np.inf
        x, y = it["xy"]
        return 6371 * math.hypot(math.radians(x - lon) * math.cos(math.radians(lat)), math.radians(y - lat))
    dist = np.array([km(it, lo, la) for it, lo, la in zip(item, df["lon"], df["lat"])])
    df["wd_match"] = ((dist <= MATCH_KM) & ~p31.map(lambda c: bool(c & PERSON_P31))).astype(int).values
    df["sitelinks"] = np.where(df["wd_match"] == 1, item.map(lambda it: it.get("n", 0)), 0)
    df["unesco"] = (df["unesco"].astype(bool) | ((df["wd_match"] == 1) & item.map(lambda it: it.get("whc", False)))).astype(int)
    print(f"Wikidata link not confirmed (no coordinates within {MATCH_KM} km, or a person) for "
          f"{(df['wd_match'] == 0).sum()} places: kept, prominence 0")
    df = df[~drop.values]
    df["name"] = df["name"].fillna(df["wikidata"].map(lambda q: wd.get(q, {}).get("en")))
    df = df[df["name"].notna()]

    pts = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df["lon"], df["lat"]), crs="EPSG:4326")
    adm3 = adm3mod.load_adm3()
    pts = gpd.sjoin(pts, adm3[["adm3_id", "adm3_name", "district", "state", "geometry"]], predicate="within", how="inner")
    pts = pts.drop(columns="index_right").reset_index(drop=True)
    pts["gfsm_class"], pts["pct_high_250m"] = sample_gfsm(pts)
    pts["pct_high_250m"] = pts["pct_high_250m"].round(1)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pts = pts.sort_values("sitelinks", ascending=False).reset_index(drop=True)
    pts.drop(columns="geometry").to_csv(OUT_DIR / "india_places.csv", index=False)
    pts.to_file(OUT_DIR / "india_places.geojson", driver="GeoJSON")

    feats = [{"type": "Feature", "geometry": {"type": "Point", "coordinates": [round(r.lon, 5), round(r.lat, 5)]},
              "properties": {"n": r.name, "c": r.category, "t": r.subtype, "w": r.wikidata, "k": int(r.sitelinks),
                             "u": int(r.unesco), "g": int(r.gfsm_class),
                             "h": None if pd.isna(r.pct_high_250m) else float(r.pct_high_250m),
                             "a": r.adm3_name, "d": r.district, "s": r.state}}
             for r in pts.itertuples()]
    VIEW.write_text(json.dumps({"type": "FeatureCollection", "features": feats}, ensure_ascii=False, separators=(",", ":")),
                    encoding="utf-8")

    summ = pts.groupby("category").agg(n=("osm", "size"), high=("gfsm_class", lambda c: (c >= 4).mean() * 100))
    print(summ.reindex(list(CATEGORIES)).round(1).to_string())
    print(f"{len(pts)} places in India; {VIEW} {VIEW.stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
