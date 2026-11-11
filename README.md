# India Tehsil Flood Atlas

**Live map:** https://stardustcafe.github.io/india-tehsil-flood-atlas/

Exploratory research coupling the Global Flood Susceptibility Map (GFSM v1, 30 m) with
socio-economic and demographic indicators to build an open flood vulnerability map for India.

## Data

| Dataset | Source | Local path |
| --- | --- | --- |
| GFSM v1 tiles (1°×1° GeoTIFF, EPSG:3395, 30 m) | Zenodo [10.5281/zenodo.20568218](https://doi.org/10.5281/zenodo.20568218) | `data/raw/gfsm/tiles_india/` |
| GFSM tile index | same Zenodo record | `data/raw/gfsm/tile_index/` |
| India boundary (India point-of-view) | Natural Earth 10m `admin_0_countries_ind` | `data/raw/boundaries/` |
| Sub-district (ADM3), district, state boundaries | [geoBoundaries](https://www.geoboundaries.org) gbOpen IND (LGD-derived, ODbL 1.0) | `data/raw/boundaries/geoboundaries/` |
| Population 2020, GPW v4.11 UN WPP-adjusted count, 30″ | [NASA SEDAC / CIESIN](https://doi.org/10.7927/H4PN93PB) (Earthdata login, CC BY 4.0) | `data/raw/population/gpw/` |
| Climate projections, CMIP6 0.25° ensemble (Rx5day, Rx1day, precipitation) | [World Bank CCKP](https://climateknowledgeportal.worldbank.org/) (`wbg-cckp` S3 bucket, CC BY 4.0) | `data/raw/climate/cckp/` |
| Places to visit (tourism, heritage, religious, nature) | OpenStreetMap via Overpass (ODbL) + Wikidata sitelinks (CC0) | `data/raw/places/` (API cache) |
| Population 2020, GHS-POP R2023A, EPSG:4326 30″ | [EC JRC GHSL](https://ghsl.jrc.ec.europa.eu/download.php?ds=pop) (no login, CC BY 4.0) | `data/raw/population/ghs/` |

## Outputs

`scripts/02_adm3_flood_susceptibility.py` → `data/processed/`

- `india_adm3_gfsm.geojson`: 6,824 sub-districts with GFSM class areas (km², latitude-corrected for
  Mercator distortion), class shares (%), `pct_high_or_above`, `mean_class`, `coverage`, `qa_flag`.
- `india_adm3_gfsm.csv`: the same attributes without geometry.
- `india_adm3_gfsm_web.geojson`: geometry simplified to 0.002° for web maps.

`qa_flag`: `ok`; `no_gfsm` (GFSM has no data, mostly Lakshadweep islands); `sliver` (degenerate source polygon).
District/state names come from a spatial join and may differ from official parent units where boundaries changed.

GFSM classes: 0 NoData, 1 Very Low, 2 Low, 3 Moderate, 4 High, 5 Very High.

### Population exposure

`scripts/06_population_exposure.py --source {gpw,ghs}` → `data/processed/`

- `india_adm3_pop_{source}.csv`: people per sub-district by GFSM class, `pop_high_or_above`, `pct_pop_high_or_above`.
- `india_adm3_flood_exposure.csv`: the GFSM table joined with every source, columns prefixed `gpw_` / `ghs_`.
  (`india_adm3_flood_exposure.geojson` adds geometry; it is not committed because of its size. Join the CSV
  to `india_adm3_gfsm.geojson` on `adm3_id`.)
- `{source}2020_india_pop.tif`, `{source}2020_india_pop_high.tif`: 30″ grids of people and of people on
  High + Very High land, India only.

Method: every 1°×1° GFSM tile holds exactly 120×120 population cells. Each 30 m pixel is assigned to its
(sub-district, population cell, class), and a cell's people are shared out in proportion to the true area of
each part over the cell's GFSM land area (classes 1–5, including land across the border, so border cells are
not over-counted). GFSM NoData pixels, mostly water, get no people.

The two sources answer different questions. GPW spreads each Census 2011 sub-district's count evenly over its
area, so its exposure is close to population × High+ area share. GHS-POP distributes people by built-up
volume, so it follows where settlements actually are. Report both.

### Climate projections

`scripts/09_climate_projections.py` → `data/processed/india_adm3_climate.csv`, `viewer/data/climate.json`

Four CMIP6 scenarios: SSP1-2.6, SSP2-4.5, SSP3-7.0 and SSP5-8.5, the successors of RCP2.6/4.5/7.0/8.5.
Four periods: 2020–39, 2040–59, 2060–79 and 2080–99. Changes are against 1995–2014, from the World Bank
CCKP multi-model ensemble (median and 10th/90th percentile across about 30 models). Per sub-district,
area-weighted:

- `d_rx5day_pct`: change in the wettest 5-day rainfall of the year, the main driver of river and widespread flooding.
  The model range is in `_p10_pct` / `_p90_pct`, and `robust_rx5day_pct` is the % of area passing the
  IPCC AR6 robustness test.
- `d_rx1day_pct`: change in the wettest single day (flash and urban flooding).
- `d_pr_pct`: change in mean precipitation.
- `flood_pressure` = `pct_high_or_above` × (1 + `d_rx5day_pct`/100).

**What this is and is not.** GFSM is a static susceptibility map. The flood-pressure index scales today's
High+ share by the projected intensification of extreme rainfall, so it ranks where heavier rain meets
flood-prone land. It is not a projection of flood probability or depth: runoff, river discharge, upstream
rainfall, dams, embankments, land-use change and sea-level rise are not modelled. Models disagree widely
over the Indian monsoon. "Models agree" in the viewer means at least 90% of models project an increase.

By 2080–99 the median sub-district's wettest 5-day rainfall rises by about 10% (SSP1-2.6), 13% (SSP2-4.5),
20% (SSP3-7.0) and 26% (SSP5-8.5).

### Places to visit

`scripts/08_tourism_places.py` → `data/processed/india_places.{csv,geojson}`, `viewer/data/places.geojson`

OSM features tagged as tourism (attraction, museum, viewpoint, zoo, theme park), historic (fort, palace,
monument, archaeological site, ...), places of worship, beaches, waterfalls, parks and protected areas, kept
only when they carry a `wikidata` tag. OSM has no visitor numbers, so the Wikidata link is the notability
filter and the number of Wikipedia/Wikimedia editions about a place (`sitelinks`) is the prominence score.
Each place is put in one of ten categories and gets its sub-district, the GFSM class of the 30 m pixel at its
point (`gfsm_class`) and the share of High + Very High land within 250 m (`pct_high_250m`). Large sites
(parks, forts) are represented by their centre point.

## Reproduce

```bash
python -m venv .venv
.venv/Scripts/pip install geopandas rasterio pyogrio requests
.venv/Scripts/python scripts/01_download_gfsm_india.py
.venv/Scripts/python scripts/02_adm3_flood_susceptibility.py   # ~5 min on 11 workers
.venv/Scripts/python scripts/03_build_viewer_data.py           # TopoJSON for the viewer (needs Node/npx)
bash scripts/04_build_gfsm_tiles.sh 10                          # GFSM raster tiles z4-10 (GDAL >= 3.x, ~1.5 min)
.venv/Scripts/python scripts/05_bundle_gfsm_tiles.py           # pack tiles into bundles for the viewer
# population: unzip GPW into data/raw/population/gpw/ and GHS-POP into data/raw/population/ghs/, then
.venv/Scripts/python scripts/06_population_exposure.py --source gpw --workers 3   # ~4 min each; workers limited by RAM
.venv/Scripts/python scripts/06_population_exposure.py --source ghs --workers 3
.venv/Scripts/python scripts/07_build_pop_tiles.py             # 1 km population tiles z4-8 for the viewer
.venv/Scripts/python scripts/03_build_viewer_data.py           # rerun to add population to the TopoJSON
.venv/Scripts/python scripts/08_tourism_places.py              # places to visit from OSM + Wikidata
.venv/Scripts/python scripts/09_climate_projections.py         # CMIP6 rainfall extremes (downloads ~660 MB once)
```

The raster overlay is built by warping the 30 m tiles once to EPSG:3857 at the z10 pixel size
(152.87 m) with `-r mode`, clipped to the India ADM3 outline, then cut with `gdal2tiles --xyz`.
Display tiles only: class shares shift by about ±1 percentage point from the 30 m data. Use
`data/processed` for any analysis. The 4,045 tiles are stored as 6-colour palette PNGs, concatenated into
82 bundles (`viewer/data/gfsm/*.png` + `index.json`), which the viewer reads through a custom MapLibre protocol.

## Viewer

`viewer/index.html` is a static MapLibre (WebGL) map of the sub-district results. It loads
`viewer/data/india_adm3.topo.json` (2.9 MB; boundaries simplified to 6% with mapshaper, shared borders kept).
Serve the folder rather than opening the file directly: `python -m http.server -d viewer 8000`.
It can be hosted as-is on GitHub Pages or any static host.

Tiles are selected by intersecting the GFSM tile index with India's boundary buffered by 0.1°,
then extracted from the four regional archives (N00E060, N00E080, N20E060, N20E080).

## Citation

Waleed, M., Sajjad, M., Al-Ghamdi, S. G., & Gao, M. (2026). A high-resolution global flood
susceptibility dataset derived from multi-source earth observation and geospatial data.
*Scientific Data*. https://doi.org/10.1038/s41597-026-08338-1 — Data: CC BY 4.0.
