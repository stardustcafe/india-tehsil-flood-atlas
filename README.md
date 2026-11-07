# gfsmxpop — Flood vulnerability mapping for India

Exploratory research coupling the Global Flood Susceptibility Map (GFSM v1, 30 m) with
socio-economic and demographic indicators to build an open flood vulnerability map for India.

## Data

| Dataset | Source | Local path |
| --- | --- | --- |
| GFSM v1 tiles (1°×1° GeoTIFF, EPSG:3395, 30 m) | Zenodo [10.5281/zenodo.20568218](https://doi.org/10.5281/zenodo.20568218) | `data/raw/gfsm/tiles_india/` |
| GFSM tile index | same Zenodo record | `data/raw/gfsm/tile_index/` |
| India boundary (India point-of-view) | Natural Earth 10m `admin_0_countries_ind` | `data/raw/boundaries/` |
| Sub-district (ADM3), district, state boundaries | [geoBoundaries](https://www.geoboundaries.org) gbOpen IND (LGD-derived, ODbL 1.0) | `data/raw/boundaries/geoboundaries/` |

## Outputs

`scripts/02_adm3_flood_susceptibility.py` → `data/processed/`

- `india_adm3_gfsm.geojson`: 6,824 sub-districts with GFSM class areas (km², latitude-corrected for
  Mercator distortion), class shares (%), `pct_high_or_above`, `mean_class`, `coverage`, `qa_flag`.
- `india_adm3_gfsm.csv`: the same attributes without geometry.
- `india_adm3_gfsm_web.geojson`: geometry simplified to 0.002° for web maps.

`qa_flag`: `ok`; `no_gfsm` (GFSM has no data, mostly Lakshadweep islands); `sliver` (degenerate source polygon).
District/state names come from a spatial join and may differ from official parent units where boundaries changed.

GFSM classes: 0 NoData, 1 Very Low, 2 Low, 3 Moderate, 4 High, 5 Very High.

## Reproduce

```bash
python -m venv .venv
.venv/Scripts/pip install geopandas rasterio pyogrio requests
.venv/Scripts/python scripts/01_download_gfsm_india.py
.venv/Scripts/python scripts/02_adm3_flood_susceptibility.py   # ~5 min on 11 workers
.venv/Scripts/python scripts/03_build_viewer_data.py           # TopoJSON for the viewer (needs Node/npx)
bash scripts/04_build_gfsm_tiles.sh 10                          # GFSM raster tiles z4-10 (GDAL >= 3.x, ~1.5 min)
.venv/Scripts/python scripts/05_bundle_gfsm_tiles.py           # pack tiles into bundles for the viewer
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
