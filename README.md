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

### Flood zones (1-in-100-year)

`scripts/10_aqueduct_floods.py` → `data/processed/india_adm3_floodzone.csv`, `viewer/data/flood.json`, `viewer/data/flood_*/`

WRI Aqueduct Floods v2 inundation maps (Ward et al. 2020, CC BY 4.0), about 1 km, without flood defences, for the
1-in-100-year flood: river (GLOFRIS; baseline 1960–99 climate, and 5 CMIP5 models under RCP4.5 and RCP8.5
for 2030, 2050 and 2080) and coast (median sea-level rise, with land subsidence). A cell is in the flood
zone when river or coastal depth is above zero. Its 2020 population is shared among the sub-districts it
overlaps (10" rasterisation). River results are worked out per model; the CSV gives the median and the
lowest/highest of the 5 models. Population is held at 2020.

| People in the flood zone (crore) | GPW | GHS-POP |
|---|---|---|
| Today | 18.2 | 14.6 |
| RCP4.5, 2050 / 2080 | 19.4 / 20.9 | 15.3 / 16.6 |
| RCP8.5, 2050 / 2080 | 20.7 / 22.5 | 16.7 / 18.4 |

The viewer matches these to the CMIP6 scenarios by nearest RCP (SSP1-2.6 and SSP2-4.5 → RCP4.5; SSP3-7.0
and SSP5-8.5 → RCP8.5) and period (2020–39 → 2030, 2040–59 → 2050, 2060–79 and 2080–99 → 2080), and labels
which map is shown. Without defences the figures are an upper bound where embankments hold; the 1 km cells
also count everyone in a cell that is only partly flooded.

### Projected population (SSPs)

`scripts/11_ssp_population.py` → `data/processed/india_adm3_ssp_pop.csv`, `viewer/data/popssp.json`

1 km population grids for SSP1, SSP2, SSP3 and SSP5 (Wang, Meng & Long 2022, Scientific Data, CC BY 4.0),
read year by year out of the 2.5 GB figshare zips with HTTP range requests. Our 2020 counts (GPW / GHS-POP)
stay the anchor and the SSP grids supply growth: each sub-district's population growth since 2020, and
separately the growth of people inside the future flood zone (cells flooded by at least 3 of 5 river models,
or by the coastal map), so that growth concentrated on floodplains counts. SSPs are paired with the flood
maps by nearest RCP (SSP1/2 → RCP4.5, SSP3/5 → RCP8.5).

| People in the flood zone, GPW (crore) | 2030 | 2050 | 2080 |
|---|---|---|---|
| SSP1-2.6 | 20.4 | 22.6 | 22.0 |
| SSP2-4.5 | 20.9 | 24.9 | 27.4 |
| SSP3-7.0 | 22.0 | 29.7 | 40.1 |
| SSP5-8.5 | 20.7 | 24.1 | 23.8 |

Today: 18.2 crore. Under SSP2-4.5 by 2080 the increase of 9.2 crore splits into 2.6 crore from the larger
flood zone (climate) and 6.6 crore from population growth. The viewer shows this split for any scenario.

### Places to visit

`scripts/08_tourism_places.py` → `data/processed/india_places.{csv,geojson}`, `viewer/data/places.geojson`

OSM features tagged as tourism (attraction, museum, viewpoint, zoo, theme park), historic (fort, palace,
monument, archaeological site, ...), places of worship, beaches, waterfalls, parks and protected areas, kept
only when they carry a `wikidata` tag. OSM has no visitor numbers, so the Wikidata link is the notability
filter and the number of Wikipedia/Wikimedia editions about a place (`sitelinks`) is the prominence score.
Each place is put in one of ten categories and gets its sub-district, the GFSM class of the 30 m pixel at its
point (`gfsm_class`) and the share of High + Very High land within 250 m (`pct_high_250m`). Large sites
(parks, forts) are represented by their centre point.

### Sea-level rise

`scripts/12_sea_level.py` → `data/processed/india_adm3_sealevel.csv`, `viewer/data/slr.json`, `viewer/data/lowland/`

IPCC AR6 relative sea-level projections (Garner et al. 2021, CC BY 4.0; median and likely range, medium
confidence, against 1995–2014) at the nearest 1° grid point, combined with the 30 m DeltaDTM coastal terrain
model (Pronk et al. 2024, CC BY 4.0). Every 30 m pixel below 3 m is binned in 5 cm steps per 30″ cell, so each
coastal sub-district (482 with at least 1 km² below 3 m) gets the land and the 2020 population below any water
level. Tide-gauge points are left out because their land-motion terms disagree within a few km. Only the
needed members are read out of the 9 GB AR6 and 17 GB DeltaDTM archives with HTTP range requests
(`scripts/remote_zip.py`).

| Median coastal sub-district | 2050 people below (GPW, lakh) | 2100 rise | 2100 people below (lakh) |
|---|---|---|---|
| SSP1-2.6 | 7.5 | 0.39 m | 10.6 |
| SSP2-4.5 | 7.7 | 0.49 m | 12.7 |
| SSP3-7.0 | 7.8 | 0.61 m | 15.3 |
| SSP5-8.5 | 8.1 | 0.72 m | 18.4 |

Today 5.5 lakh people (GPW) live on 708 km² below mean sea level, led by Kuttanad in Kerala. This is a
bathtub estimate with 2020 population: no check that land connects to the sea, and most of it is behind
embankments, so it shows land that depends on protection rather than land that will flood. Tides and surges
are in the Aqueduct coastal flood zone. The viewer scales these counts by SSP growth.

### River flood frequency (return-period shift)

`scripts/13_return_period.py` → `data/processed/india_adm3_returnperiod.csv`, `viewer/data/rp.json`, `viewer/data/rivers.geojson`

ISIMIP3b daily river discharge (CC BY 4.0) from three global hydrological models (CWatM, WaterGAP2-2e,
MIROC-INTEG-LAND), each driven by five bias-adjusted CMIP6 models: 15 runs at 0.5°. The ISIMIP files API masks
each file to India on the server, so about 20 MB per decade is downloaded instead of about 1 GB, and only annual
maxima are kept. For every run and river cell, Gumbel distributions (L-moments) are fitted to the annual maximum
daily flow of 1985–2014 and of each future 20-year period. The future return period of today's 100-year flow is
then compared with the run's own past, so model bias largely cancels (after Hirabayashi et al. 2013). River cells
drain at least 10,000 km² on the DDM30 network the models route on. Each sub-district takes the largest river
cell covering at least 10% of its area, which gives 3,613 sub-districts. ISIMIP3b has no SSP2-4.5 runs, so
SSP2-4.5 is the geometric mean of SSP1-2.6 and SSP3-7.0, labelled as an estimate.

| Median river sub-district: today's 100-year flood comes every … years | 2020–39 | 2040–59 | 2060–79 | 2080–99 |
|---|---|---|---|---|
| SSP1-2.6 | 89 | 57 | 42 | 48 |
| SSP2-4.5 (estimated) | 80 | 61 | 43 | 44 |
| SSP3-7.0 | 73 | 75 | 38 | 41 |
| SSP5-8.5 | 58 | 78 | 30 | 13 |

By 2080–99, today's 100-year flood comes more often in 83–100% of river sub-districts, depending on the
scenario. Values for individual 20-year periods are noisy: a 100-year event estimated from 20 years of maxima
varies, so read trends across periods and scenarios rather than single cells. The models leave out dams and
embankments, and the indicator covers river flow only, not exposure.

### Cyclone wind

`scripts/14_cyclone_wind.py` → `data/processed/india_adm3_cyclone.csv`, `viewer/data/cyclone.json`, `viewer/data/wind_today/`, `viewer/data/wind_future/`

STORM synthetic tropical cyclones (Bloemendaal et al. 2020, 2022, CC0): 10 km wind speeds at fixed return
periods for the North Indian basin, from 10,000 years of synthetic cyclones resampled from 1980–2017 tracks. The
future version adds the 2015–2050 climate signal of four high-resolution HighResMIP models under high
emissions. The atlas takes the median of the four models and their range. Winds are 10-minute sustained, slightly
below IMD's 3-minute winds, so the IMD categories applied here are conservative. Only the North Indian basin
GeoTIFFs are read out of the archives.

| People where the 1-in-100-year wind reaches very severe cyclone strength (≥ 33 m/s, 118 km/h) | Today | 2015–2050 |
|---|---|---|
| GPW (crore) | 16.8 | 14.8 |
| GHS-POP (crore) | 17.3 | 15.0 |

The models project fewer intense Bay of Bengal cyclones, so the strongest Odisha winds fall from about 44 to
39 m/s, while the 1-in-10-year wind rises in many places. Future values are tied to one scenario and period, so
they do not follow the viewer's scenario selector.

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
.venv/Scripts/python scripts/10_aqueduct_floods.py             # Aqueduct flood zones + depth tiles (streams ~3 GB, keeps India crops)
.venv/Scripts/python scripts/11_ssp_population.py            # SSP population growth (streams ~3 GB of the figshare zips)
.venv/Scripts/pip install h5py
.venv/Scripts/python scripts/12_sea_level.py --workers 8     # AR6 sea level + DeltaDTM (~1.4 GB download, ~1 h)
.venv/Scripts/python scripts/13_return_period.py --workers 8 # ISIMIP river flow, 465 India cut-outs (~10 GB, ~2 h)
.venv/Scripts/python scripts/14_cyclone_wind.py             # STORM cyclone wind return periods (~5 MB of GeoTIFFs)
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

The sidebar has four tabs. **Layers** lists every map layer in drawing order (places, labels and rivers,
boundaries, rainfall change, 1-in-100-year flood zone, river flood frequency, cyclone wind, low-lying coastal land, population, GFSM 30 m, sub-district statistics, and a
Sentinel-2 cloudless satellite basemap by EOX for 2024 or 2016), each with its own switch, opacity and settings.
**Findings** gives headline numbers as plain sentences for India or the chosen state, today and under the chosen
scenario, plus rankings. **Glossary** defines the terms; underlined terms anywhere in the viewer link to it.
**About** holds the purpose, caveats and sources. Clicking a sub-district opens a profile with Today and Future tabs.

For speed on GitHub Pages the first view loads only the boundaries (0.9 MB gzipped). Climate (0.9 MB),
flood-zone (0.6 MB), projected population (0.5 MB), river flood frequency (0.4 MB), cyclone wind (0.05 MB), sea-level (0.2 MB) and places (0.1 MB) data load the first time a layer, metric, tab or profile needs them, and
raster layers fetch only the tile bundles in view.

Tiles are selected by intersecting the GFSM tile index with India's boundary buffered by 0.1°,
then extracted from the four regional archives (N00E060, N00E080, N20E060, N20E080).

## Citation

Waleed, M., Sajjad, M., Al-Ghamdi, S. G., & Gao, M. (2026). A high-resolution global flood
susceptibility dataset derived from multi-source earth observation and geospatial data.
*Scientific Data*. https://doi.org/10.1038/s41597-026-08338-1 — Data: CC BY 4.0.
