# Validation of the pipeline against draft v84

Each step is checked against the values printed in draft v84 before the next step is built. Differences are recorded here, never tuned away.

## Step 2 — ICESat-2 segments (17 Nov 2026)
- Screening of the cached SlideRule output gives 3,074,432 segments, identical to the v84 working file (`Waikato_sampled.parquet`).

## Step 3 — ICESat-2 rates (17 Nov 2026; first run on the v84 segment file, rerun on step-2 output below)

Reproduced exactly (value, SE, passes):

| Quantity | v84 | Pipeline |
|---|---|---|
| Organic-Soil farmland vs mineral-soil farmland, 2018–2026 | −19.3 ± 1.56, 134 passes | −19.3 ± 1.56, 134 |
| Kopuatai bog vs mineral-soil farmland within 5 km | +14.5 ± 4.2, 20 | +14.5 ± 4.18, 20 |
| Mineral-soil farmland within 5 km of Kopuatai, absolute | +3.2 ± 2.3 | +3.2 ± 2.34 |
| Organic, Oct 2018 – Jun 2022 / Jul 2022 – Oct 2026 | −27.5 ± 4.1 / −16.4 ± 5.0 | −27.5 ± 4.11 / −16.4 ± 4.95 |
| Segmented: before / after / difference CI | −28.1 / −15.0 / −25.1 to −1.2 | −28.1 / −15.0 / 1.2 to 25.1 (after − before) |
| Bog segmented: before / step / after | −10.9 / 52 / +15.5 | −10.9 / 52.1 / +15.5 |
| Seasonal amplitude, Organic / mineral | 29 (15–44) / 10 | 28.6 (14.6–44.3) / 9.7 |
| 2021 survey bias: regional / farm − built-up; 2024 farm − built-up | +3.6 / 4.2 / 0.8 cm | +3.6 / 4.2 / 0.8 cm |

Not reproduced, because v84 took these from an older script (`waikato_check.py`) with different class definitions from the rest of the paper:

| Quantity | v84 | Pipeline | Cause |
|---|---|---|---|
| Mineral-soil farmland, absolute | +1.4 ± 1.38, 264 passes | +1.2 ± 1.40, 228 | v84: Nth_Waikato polygon (extends beyond the study rectangle), Organic Soils from smap_organic shapefile, outlier median over all land; pipeline: study rectangle and S-map soil order, as the Methods state |
| Mineral-soil farmland vs built-up land | +2.2 ± 2.08 | +2.5 ± 2.06 | same |
| Built-up land, absolute | −1.8 | −1.5 | same |
| Hauraki Plains mineral-soil farmland, absolute | +0.9 ± 2.15 | +2.9 ± 1.69 | v84: lon/lat box 175.35–175.65° E, 37.12–37.45° S; pipeline: NZTM box of Section 2 (1805–1846 km E, 5840–5886 km N) |

Consequence: the repeat-LiDAR tie uses the mineral-soil absolute rate; +1.2 instead of +1.4 shifts the tied rates by −0.2.

## Text/code discrepancies found while reading the v84 code (to resolve in step 4)
- LiDAR − Copernicus: v84 text says the survey-bias correction is applied to the 2021 LiDAR only; `dem_hex.py` also corrects the 2024 LiDAR by its own 10 km bias grid (`C24 = L24 + bias24 − COP`).
- LiDAR − Copernicus: v84 text says both LiDAR DEMs are LINZ 1 m tiles averaged to 30 m; the 2024 surface in `dem_hex.py` is a 10 m mosaic (`LDS_10m_DEM_merged_downsampled.tif`) averaged to 30 m.

## Steps 2–3 end to end (17 Nov 2026)
- Step 2 run from the raw SlideRule download, LINZ tiles and WorldCover (256 s): every column of `segments.parquet` is identical to the v84 working file (3,074,432 rows; heights, LiDAR heights and slopes differ by 0.0).
- Step 3 on the step-2 output gives results identical to the test run above.
