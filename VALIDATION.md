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

## Decision (user, 17 Nov 2026): option 1
The paper takes the pipeline's mineral-soil values (study rectangle, S-map soil order, NZTM Hauraki box): mineral-soil farmland +1.2 ± 1.4 absolute, +2.5 ± 2.1 relative to built-up land, Hauraki Plains +2.9 ± 1.7. The tie uses +1.2. To be applied in v85 with the step-4 results.
- Section 3.3 compares unlike quantities: "ICESat-2 heights exceed the 2021 LiDAR by 3.6 cm on average [absolute, grassland grid] ... against 0.8 cm for the 2024 LiDAR [farmland minus built-up, same passes]". Like for like: absolute 2021 +3.6 cm vs 2024 −2.5 cm (grassland, 10 km grid regional; farmland median −1.8 cm); farmland minus built-up 2021 4.2 cm vs 2024 0.8 cm. v84 (dem_hex.py) corrects both surveys with their own 10 km grids; the pipeline does the same and the text must say so.

## Step 1 — grids (17 Nov 2026, 362 s)
- Copernicus NZVD2016 30 m grid and 2021 LiDAR 30 m grid: identical to the v84 inputs (difference 0.000 m everywhere, same coverage).
- 2024 LiDAR 30 m grid: v84 used a ~10 m mosaic in geographic coordinates (`LDS_10m_DEM_merged_downsampled.tif`, EPSG:4167), resampled twice; the pipeline averages the LINZ 1 m tiles, as the Methods state. On flat farmland: median difference −2 mm, median absolute difference 2.7 cm (about 2.4 mm/yr per cell over 11.2 years), 1.8% of cells > 0.5 m.

## Step 4a — LiDAR − Copernicus (17 Nov 2026)
| Quantity | v84 | Pipeline |
|---|---|---|
| Organic-Soil farmland (median of referenced grassland cells) | −14.3 (95% CI −17.6 to −10.8) | −14.1 (−16.8 to −10.7) |
| Unreferenced mineral reference | +24.3 | +24.3 |
- Correcting the 2024 survey for its bias as well (as the v84 code did) gives the same −14.1: the correction is nearly constant within each reference window. The pipeline follows the Methods (2021 survey only).
- Remaining 0.2 difference: 2024 surface (above) and bootstrap seed.

## Step 4b — repeat LiDAR (17 Nov 2026, 461 s incl. capture dates from 2,891 LAZ tiles)
96 m grid from the study-area corner (v84: from tile bounds, 16 m offset); capture dates = median GPS time per cell (every 20th return).

18 Nov 2026 (paper v93): moved to 30 m cells (median of 3 x 3 cells of 10 m; the 30 m grid of the other maps) to use one analysis grid. 96 m vs 30 m: northern Torehape -15.3 vs -15.4, southern Torehape -28.7 vs -28.4, Kopuatai margins -12.5 vs -12.5 mm/yr (same campaign); tied Organic-Soil farmland -18.5 vs -17.4 (16-placement range -20.3 to -15.6). The table below is the 96 m run.
| Quantity | v84 | Pipeline |
|---|---|---|
| Area differenced | 335 km² | 337 km² |
| Same campaign, northern Torehape (Feb 2015 → 2024) | −13.8 cm, 9.0 yr, −15.4 | −13.7 cm, 9.0 yr, −15.3 |
| Same campaign, southern Torehape (2017–18 → 2024) | −17.3 cm, 6.1 yr, −28.8 | −17.2 cm, 6.1 yr, −28.7 |
| Same campaign, drained margins of Kopuatai | −7.2 cm, 5.9 yr, −12.1 | −7.4 cm, 5.9 yr, −12.5 |
| Mineral-soil farmland untied, 2015 / 2017–18 flights | +2.7 / +10.8 | +2.8 / +10.7 |
| Tied (5 km), Organic-Soil farmland | −19.0 | −19.0 (tie +1.2, option 1) |
| Tied, mineral / Hauraki Plains mineral | +1.1 / +1.1 | +1.1 / +1.1 |
| Tied, Kopuatai bog | −10.3 | −10.7 |
| Grid-placement spread, Organic: 5 km / 10 km | −20.0 to −16.5 / −22.1 to −14.3 | −20.5 to −16.0 / −21.6 to −15.3 |
All within 0.5 mm/yr; differences from the 16 m grid offset and the tie value (+1.2 instead of +1.4).

## Step 5 — InSAR (17 Nov 2026)
- LiCSBAS step 15 rerun on a copy of the 14 Jan 2025 time series with the Methods thresholds (fewer than 141 unwrapped ifgs, vstd > 2 mm/yr, n_gap 30, n_loop_err 100, resid_rms 10, others default) keeps 1,157,797 pixels, exactly the published new_params mask. (The job log's step-15 run used LiCSBAS defaults for n_unw and vstd; that was an earlier run, not the map used. The Methods text is correct.)
- The pipeline's map-referencing gives cell-for-cell the same result as the v84 code on the same 30 m grid.
| Quantity | v84 | Pipeline |
|---|---|---|
| InSAR zones ≤ −40 / −50 / −60 (km²) | 114.6 / 64.0 / 36.6 | 114.6 / 64.0 / 36.6 |
| Organic-Soil farmland (240 m map at segment locations) | −20.1 | see insar.json (soil-order classes; v84 used the peat shapefile and Nth_Waikato polygon) |
- Mask sensitivity on the pipeline chain: no mask gives Organic −20.6 (masked −20.5) and zones 1–4% larger (−10 to −60 mm/yr); LiCSBAS default n_unw/vstd thresholds give zones 2–20% smaller.
- Kopuatai bog InSAR: pipeline −10.5 (240 m farmland map at bog segments) vs v84 −26.0 (peat_bodies.py: wetland cells referenced to mineral farmland); definition to be ported in step 6.

## Decision (user, 17 Nov 2026): no InSAR mask
LiCSBAS step 15 is not applied (config `insar.mask: false`), for ease of reproduction; to be stated in Table S1 (step 5) and the Methods.
Unmasked pipeline values: Organic-Soil farmland −20.6 (masked −20.5; v84 −20.1 on v84 classes), mineral +0.5, Hauraki Plains +0.6;
zones ≤ −10/−20/−30/−40/−50/−60: 636.5 / 331.9 / 202.3 / 119.1 / 67.1 / 37.8 km² (v84 masked: 641 / 328 / 195 / 114 / 63 / 36, Table 4).

## Step 6a — hexagons (17 Nov 2026)
- is2rate hexagon rates from the step-2 download are identical to the paper's (572 hexagons, median difference 0.0, 95th percentile 0.1 mm/yr), although the paper's run used the reproduction's separate SlideRule fetch.
- v84 kept only hexagons that fell on an older hexagon layer (`is2_hex_smap_3km.gpkg`, v54_chain.py): an undocumented extent rule. On that layer the pipeline gives 453 hexagons, r = 0.61 (LiDAR − Copernicus) and 0.57 (InSAR, unmasked) vs v84 454, 0.62, 0.59.
- Pipeline rule (Methods): hexagons at least 25% farmland with their centre in the study rectangle: 471 hexagons, r = 0.58 and 0.55.

## Step 6b — maximum rates, Table 4 (17 Nov 2026, 54 is2rate runs, 55 min)
- Distance bands: ICESat-2 rates identical to v84 in all three periods (largest difference 0.3 mm/yr, 2022–2026, 250–500 m); areas identical.
- LiDAR − Copernicus band medians: −8.5, −12.1, −18.3, −22.2, −23.5, −25.7 vs v84 −8.5, −12.4, −18.5, −22.4, −23.9, −26.1 (2024 surface from 1 m tiles).
- InSAR zones (unmasked pipeline map; areas in km²): ≤ −40 mm/yr 119.1 km² (v84 114.5); 2018–2026 rates within 2 mm/yr at ≤ −25 and stricter;
  at loose thresholds the InSAR-period rate moves more: ≤ −10 −22.4 ± 4.2 (v84 −15.1 ± 3.8), ≤ −20 −28.6 ± 4.9 (v84 −31.9 ± 4.8).
  Sensitivity run with the masked pipeline map: tests/zones_masked.py (pending).
- Masked-map sensitivity (tests/zones_masked.py): masked vs unmasked pipeline zones agree within about 2 mm/yr (largest: ≤ −60, InSAR period, −53.0 vs −47.1, 0.6 SE).
  The larger differences from v84 (e.g. ≤ −10, InSAR period: v84 −15.1, masked pipeline −20.7; ≤ −40: −58.8 vs −51.9) come from the v84 zone map,
  `Waikato_InSAR_rate_rel_local_pasture_240m.tif`, an older map version that differs cell by cell from the stated chain (5th/95th percentile ±6 mm/yr), not from the mask.

## Steps 4–5 additions (17 Nov 2026)
- Repeat LiDAR 240 m map (Figure 3b) built in step 4; Section 4.2 InSAR quarters (Hauraki Plains mineral-soil farmland, 240 m cells): InSAR −8.4 to +7.3, repeat LiDAR −1.5 to +1.3, 1,543 cells (v84: −9.1 to +7.3, −2.0 to +0.8, 1,708).
- Kopuatai bog InSAR (wetland cells on the Kopuatai body, referenced locally): −28.4 unmasked (v84 −26.0, masked).
- Mineral-soil farmland InSAR relative to Cambridge built-up land (3 km of the centre): −18.8 unmasked (v84 −18.1, published map without rescaling).

## Step 7 — weather: reproduced (rain 1,190 / 1,632 mm/yr; soil moisture 0.360 / 0.394; water years identical; ERA5 1,006 / 1,323 / 1,176; LiDAR − Copernicus interval 1,159 vs 1,160 from the 2024.24 vs 2024.2 date).

## Step 8 — Canterbury and Southland
| Quantity | v84 | Pipeline |
|---|---|---|
| Southland farmland absolute | −1.6 (−3.4 to +0.1), 244 | −1.5 (−3.4 to +0.4), 246 |
| Southland rel. built-up | +1.1 (−5.2 to +6.9), 62 | +3.6 (−1.1 to +8.3), 53 |
| Southland fastest-InSAR farmland rel. built-up | +0.1 | +2.3 |
| Canterbury farmland absolute | −1.7 (−4.0 to +0.8), 222 | −1.6 (−3.3 to +0.2), 222 |
| Canterbury rel. built-up | +2.3 (−2.1 to +6.9), 83 | +1.4 (−2.6 to +5.4), 72 |
| InSAR on farmland / built-up | −16.3/−0.2, −6.7/−0.1 | identical |
v84 used 15 segments per pass for farmland and 8 for built-up land, and bootstrap CIs; the pipeline uses the Methods' 10 and t-distribution CIs.
Step 3 aligned to 10 for built-up land: Waikato mineral-soil farmland relative to built-up land +2.9 ± 1.9 (100 passes); built-up absolute −1.9.

## Step 9 — outputs: Table3/Table4.xlsx, paper_values.json, plot workbooks, GeoTIFF/GeoPackage maps, waikato_is2.qgz (checked by rendering).

## Decision (user, 17 Nov 2026): tied repeat LiDAR = median over 16 placements of the 5 km grid
Default placement alone lay near the fast end of the spread (independent check of v83). Medians: Organic-Soil farmland −18.5 (placements −20.5 to −16.0), mineral +0.9, Hauraki Plains mineral +0.9, Kopuatai bog −12.1 (default placement −19.0, +1.1, +1.1, −10.7). Draft v86.

## Decision (user, 17 Nov 2026): LiDAR − Copernicus on farmland; no grassland category (draft v87)
Fifth independent check: Section 3.4 (farmland) and Table S1 (grassland) disagreed. The grassland rule was inherited from dem_hex.py, not chosen.
- Crop-canopy test: mineral-soil cropland +4.3 vs grassland +3.3 (no false subsidence on cropland), so no bias from crops in the Copernicus DSM.
- Cropped Organic Soils subside faster: ICESat-2 −21.6 ± 2.2 (cropland) vs −17.9 ± 1.6 (grassland); LiDAR − Copernicus −20.2 vs −15.2. Grassland-only excluded them.
- Now farmland cells and farmland reference everywhere (Table value, bands, hexagons with widened reference, cap test); survey bias measured on mineral-soil farmland (+3.5 cm; grassland +3.6).
- LiDAR − Copernicus −17.2 (95% CI −20.2 to −14.0), about 11% less than ICESat-2 and within its uncertainty (was −14.1, "about 25% less"). Hexagons 472, r = 0.63 / 0.58.

## Sixth independent check (v87) and decisions (17 Nov 2026; drafts v88, v89)
- Block bootstrap: the pipeline counted 10 km squares from round 10 km lines, Table S1 from the study-area corner (5,785,000 N). Fixed to the corner: 42 squares (as the check), CI widens from about 7 to 10 mm/yr; also the bias grid. The CI depends on square placement (5.6–9.6 mm/yr wide for 5 km shifts).
- Survey bias: the absolute 2021 bias depends on the statistic (segment median +5.0 cm; median over passes +2.1; mean over passes +0.6). Correcting or not changes LiDAR − Copernicus by 0.1 (−17.2 uncorrected, −17.1 with the 4.2 cm same-pass offset, −17.3 with the old grid). Decision (user): no survey-bias correction; Section 3.3 removed, one sentence in the LiDAR − Copernicus paragraph; offsets reported as farmland minus built-up on the same passes (4.2 / 0.8 cm).
- LiDAR − Copernicus (v89): −17.2 (95% CI −21.9 to −12.1), unreferenced mineral reference +20.0; hexagons r = 0.64 / 0.58 (472).

## Every number from the pipeline (draft v90, 17 Nov 2026)
Ported into the pipeline (previously read from pre-pipeline files or written in the text): period rates, segmented fits and the separate-fit difference CI;
seasonal amplitudes; bog and local mineral values; weather periods and water years; ERA5 periods; segment count (2.2 million); S-map Organic area (770 km²);
land area of the study rectangle (8,868 km²); Table 1 overlaps (2,725 / 6,003 km²; v84 2,681 / 6,046; check 2,693 / 5,971); band areas in the study rectangle;
LiCSBAS network (505 interferograms, 118 acquisitions); annual-cycle omission (1.0 mm/yr); SE ratios; mask sensitivity (7% of Organic-Soil farmland cells; 1.0 mm/yr; zones 1–5% smaller; zone rates up to 2.7 mm/yr).
Corrected in the paper: 3 km bog limit changes the rate by 0.6 (was "< 0.1"; check 0.7); SE ratio median 0.96, range 0.79–1.02 (was 1.00, 0.81–1.00; check 0.98, 0.82–1.08).
Bugs fixed: WorldCover window for the 30 m masks began at 174.85° E, missing a 1,342 km² western strip of the study rectangle (1.5 km² of Organic Soils); hexagons 472 -> 486.
is2rate runs are now reused only when the target polygon is unchanged (md5 in the parameter file).
Differences from the sixth check explained: survey offsets (median vs mean over passes, and scatter between passes); InSAR class values (now medians of 30 m cells, as LiDAR − Copernicus: −20.3 vs check −20.5);
bog passes (reference within 5 km of the Kopuatai body vs of the bog); Hauraki quarters (240 m cells entirely > 500 m from Organic Soils vs mostly mineral). Hexagons 486 vs 477: edge rules, not pursued.
A provenance scan of v90 finds no number without a pipeline or config source other than the study-area coordinates and instrument facts.

## Seventh independent check (v90) and draft v91 (17 Nov 2026)
Reproduced: LiDAR − Copernicus CI (now ~10% wider in the check), unreferenced reference (+20.0 both), InSAR class values, InSAR-period agreement (no longer borderline), Hauraki quarters, bog, all of Table 5.
- Mask sensitivity: the check applied the mask pixel by pixel; the pipeline masked before averaging (7.1% removed). Pixel by pixel: 10.0% removed, Organic change 1.3 (check 10.0%, 1.39). Decision (user): Section 5.3 keeps one sentence with these two numbers; mask thresholds dropped from Table S1; zone test (tests/zones_masked.py) removed.
- Survey offsets: not reproducible to better than a few cm (pass selection); stated qualitatively.
- SE ratios: "within about 20% (median 0.96)"; bands: "about 2–3 mm/yr further inside".
