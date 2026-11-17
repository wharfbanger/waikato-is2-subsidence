# Raw data (not tracked)

Place under `data/raw/` (paths in `config.yaml`):

| Folder / file | Source |
|---|---|
| `icesat2/<Area>_??_atl06sr.parquet` | SlideRule ATL06-SR (step 2 and step 8 download them if missing) |
| `boundaries/Nth_Waikato.shp`, `Southland.shp`, `Canterbury.shp` | Study areas of Harvey and Rissmann (2026) and Harvey et al. (2026) |
| `soils/smap_organic_waikato.shp`, `smap_soilorder_waikato.gpkg` | S-map soil order (Manaaki Whenua – Landcare Research, August 2024) |
| `soils/Waikato_peat_bodies_stats.gpkg` | Organic Soil bodies (Kopuatai, Torehape, …) |
| `copernicus/Copernicus_DSM_COG_10_*_DEM.tif` | Copernicus GLO-30 (s3://copernicus-dem-30m) |
| `worldcover/ESA_WorldCover_10m_2021_v200_S39E174_Map.tif` | ESA WorldCover 2021 v200 |
| `proj/nz_linz_nzgeoid2016.tif`, `us_nga_egm08_25.tif` | PROJ-data (LINZ NZGeoid2016; NGA EGM2008) |
| `tile_index.json` | LINZ nz-elevation tile index (collection, tile names, bounds, base URL) |
| `weather/*.json` | Open-Meteo archive (step 7 downloads if missing) |
| `south_island/*_UD_corrected_using_regression.tif` | Published InSAR of Harvey et al. (2026) |
| `linz_dem/<collection>/` (path `linz_dem`) | LINZ 1 m DEM tiles (step 1 downloads) |
| `laz/NZ15_Huntly`, `laz/NZ24_Waikato` (path `laz`) | LINZ point clouds via OpenTopography |
