# Data licences

The code in this repository is MIT-licensed (see `LICENSE`). The data products are derived from
third-party datasets and carry their licences:

| Files | Derived from | Licence |
| --- | --- | --- |
| `data/processed/*`, `viewer/data/india_adm3.topo.json` | geoBoundaries gbOpen IND ADM3 boundaries + GFSM v1 | **ODbL 1.0** (share-alike, from geoBoundaries); GFSM attribution required |
| `viewer/data/gfsm/*` (raster tiles) | GFSM v1 | **CC BY 4.0** |
| `data/processed/india_adm3_climate.csv`, `viewer/data/climate.json` | World Bank CCKP CMIP6 0.25° ensemble, with GFSM v1 | **CC BY 4.0** |
| `data/processed/india_places.*`, `viewer/data/places.geojson` | OpenStreetMap + Wikidata, with GFSM v1 class | **ODbL 1.0** (from OpenStreetMap) |
| `data/processed/*2020_india_pop*.tif`, `viewer/data/pop_*/*` | GPW v4.11 or GHS-POP R2023A, split by GFSM v1 | **CC BY 4.0** |

## Required attribution

- Waleed, M., Sajjad, M., Al-Ghamdi, S. G., & Gao, M. (2026). A high-resolution global flood
  susceptibility dataset derived from multi-source earth observation and geospatial data.
  *Scientific Data*. https://doi.org/10.1038/s41597-026-08338-1 (data: https://doi.org/10.5281/zenodo.20568218, CC BY 4.0)
- geoBoundaries: Runfola, D. et al. (2020) geoBoundaries: A global database of political administrative
  boundaries. *PLoS ONE* 15(4): e0231866. Source: Pathways Data Pvt. Ltd., lgdirectory.gov.in. ODbL 1.0.
- GPW: Center for International Earth Science Information Network (CIESIN), Columbia University. 2018.
  Gridded Population of the World, Version 4 (GPWv4): Population Count Adjusted to Match 2015 Revision of
  UN WPP Country Totals, Revision 11. NASA SEDAC. https://doi.org/10.7927/H4PN93PB. CC BY 4.0.
- GHS-POP: Schiavina, M., Freire, S., Carioli, A., MacManus, K. (2023). GHS-POP R2023A, GHS population grid
  multitemporal (1975-2030). European Commission, Joint Research Centre (JRC).
  https://doi.org/10.2905/2FF68A52-5B5B-4A22-8F40-C41DA8332CFE. CC BY 4.0.
- Climate: World Bank Group, Climate Change Knowledge Portal (CCKP), CMIP6 0.25° collection,
  https://climateknowledgeportal.worldbank.org/. CC BY 4.0. Based on CMIP6 model output (WCRP).
- Places to visit: © OpenStreetMap contributors (ODbL); prominence from Wikidata (CC0).
- Flood zones: WRI Aqueduct Floods v2 hazard maps, CC BY 4.0. Ward, P.J., Winsemius, H.C., Kuzma, S., et al. (2020),
  Aqueduct Floods Methodology, Technical Note, World Resources Institute. www.wri.org/publication/aqueduct-floods-methodology
- Projected population: Wang, X., Meng, X., Long, Y. (2022), Projecting 1 km-grid population distributions from 2020
  to 2100 globally under shared socioeconomic pathways, Scientific Data 9, 563. https://doi.org/10.6084/m9.figshare.19608594, CC BY 4.0.
- Viewer place names and rivers: © OpenStreetMap contributors (ODbL), served by OpenFreeMap using the OpenMapTiles schema.
- Viewer satellite imagery (streamed, not redistributed): Sentinel-2 cloudless (https://s2maps.eu) by EOX IT Services GmbH,
  contains modified Copernicus Sentinel data. The 2016 mosaic is CC BY 4.0; the 2018 onwards mosaics (the viewer uses 2024) are CC BY-NC-SA 4.0, non-commercial use only.

Under ODbL, if you publicly use an adapted version of the sub-district database, you must offer that
adapted database under ODbL as well.

Boundaries are for analysis and do not represent an authoritative depiction of India's external boundaries.
