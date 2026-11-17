# waikato-is2-subsidence

Processing for Harvey (in prep.), *ICESat-2 confirms peat subsidence in the Waikato, New Zealand*. Every table, plot and map of the paper is produced by one master script from public data.

```
python run_all.py                    # all steps
python run_all.py --start 3 --end 4  # a range of steps
```

All parameters are in `config.yaml`. Steps (Table S1 of the paper):

| Step | Script | Does |
|---|---|---|
| 1 | `steps/s01_prepare.py` | LINZ 1 m LiDAR tiles, 30 m grids, Copernicus to NZVD2016, WorldCover and S-map masks, LiDAR mosaics |
| 2 | `steps/s02_icesat2.py` | SlideRule ATL06-SR segments, screening, NZVD2016 heights, labels, LiDAR heights |
| 3 | `steps/s03_icesat2_rates.py` | Per-pass target minus reference, trends, periods, segmented fit, seasonal amplitude, LiDAR survey bias |
| 4 | `steps/s04_dem_rates.py` | LiDAR − Copernicus; repeat LiDAR (capture dates from point clouds, same-campaign rates, tie to ICESat-2) |
| 5 | `steps/s05_insar.py` | LiCSBAS (WSL) velocity, vertical, local reference, 240 m map, zones (no noise mask, for reproducibility) |
| 6 | `steps/s06_compare.py` | Hexagons and maximum rates with is2rate; map hexagons; comparisons |
| 7 | `steps/s07_weather.py` | Open-Meteo rainfall and soil moisture |
| 8 | `steps/s08_south_island.py` | Canterbury and Southland |
| 9 | `steps/s09_outputs.py` | Tables (Excel), plot data (Excel), maps (GeoTIFF/GeoPackage), QGIS project |

Requirements: Python 3.12 with numpy, pandas, scipy, geopandas, rasterio, pyproj, laspy (lazrs), sliderule, openpyxl, pyyaml;
[is2rate](https://github.com/wharfbanger/is2rate) v0.3.0; LiCSBAS in WSL (step 5); QGIS 3 (optional, for the project file).
`VALIDATION.md` records the comparison of every step with draft v84.
