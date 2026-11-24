# waikato-is2-subsidence

Processing for Harvey (in prep.), *ICESat-2 confirms peat subsidence in the Waikato, New Zealand*. Every table, plot and map of the paper is produced by one master script from public data.

```
python run_all.py                    # all steps
python run_all.py --start 3 --end 4  # a range of steps
```

All parameters are in `config.yaml`; put machine-specific paths (LiDAR tiles, point clouds; LiCSBAS job, only to rerun LiCSBAS) in an untracked `config.local.yaml`, which `run_all.py` merges over it. Raw data are not included. The ICESat-2 segments and the LiCSBAS outputs used by the paper are archived on Zenodo (https://doi.org/10.5281/zenodo.24385729): copy the `*_atl06sr.parquet` and `*_failures.json` files to `data/raw/icesat2/` and unzip `insar_LiCSBAS_073D_12785_121717_2018-2022.zip` in `data/raw/` (it holds `insar/`); steps 2, 5 and 8 then use them instead of SlideRule and LiCSBAS. Otherwise step 1 downloads the LINZ DEM tiles, step 2 the ICESat-2 segments (SlideRule) and step 7 the Open-Meteo series; S-map, WorldCover, Copernicus GLO-30, NZGeoid2016, the LINZ point clouds, the LiCSAR interferograms and the published South Island InSAR are listed in `docs/data_raw.md`. Steps (Table S1 of the paper):

| Step | Script | Does |
|---|---|---|
| 1 | `steps/s01_prepare.py` | LINZ 1 m LiDAR tiles, 30 m grids, Copernicus to NZVD2016, WorldCover and S-map masks, LiDAR mosaics |
| 2 | `steps/s02_icesat2.py` | SlideRule ATL06-SR segments, screening, NZVD2016 heights, labels, LiDAR heights |
| 3 | `steps/s03_icesat2_rates.py` | Per-pass target minus reference, trends, periods, change of rate at July 2022, LiDAR survey bias |
| 4 | `steps/s04_dem_rates.py` | LiDAR − Copernicus; repeat LiDAR (capture dates from point clouds; same-flight rates of the Torehape site, no ICESat-2) |
| 5 | `steps/s05_insar.py` LiCSBAS velocity (archived outputs in `data/raw/insar/`, else LiCSBAS run with bash, in WSL on Windows), vertical, local reference, 240 m map, zones (no noise mask, for reproducibility) |
| 6 | `steps/s06_compare.py` | ICESat-2 hexagon rates and maximum rates (pipeline ICESat-2 code, lib/is2.py); map hexagons; comparisons |
| 7 | `steps/s07_weather.py` | Open-Meteo rainfall and soil moisture |
| 8 | `steps/s08_south_island.py` | Canterbury and Southland |
| 9 | `steps/s09_outputs.py` | Tables (Excel), plot data (Excel), maps (GeoTIFF/GeoPackage), QGIS project |

Requirements: Python 3.12 with numpy, pandas, scipy, geopandas, rasterio, pyproj, laspy (lazrs), sliderule, openpyxl, pyyaml;
QGIS 3 (optional, for the project file). LiCSBAS (Linux/macOS, or WSL on Windows) is needed only without the archived InSAR outputs, or with `insar.force_licsbas: true`.
`VALIDATION.md` records the comparison of every step with draft v84.
