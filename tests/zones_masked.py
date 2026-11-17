"""Sensitivity: Table 4 InSAR zones from the pipeline map WITH the LiCSBAS step-15 mask (Harvey et al. 2026 thresholds), to separate the
effect of dropping the mask from differences between the pipeline map and the v84 zone map. Writes outputs/tables/max_rates_zones_masked.json."""
import json, os, sys, numpy as np, yaml, rasterio, geopandas as gpd
sys.path.insert(0, os.getcwd())
from rasterio.transform import from_origin
from rasterio.warp import reproject, Resampling
from rasterio.features import shapes
from shapely.geometry import shape
from lib.maps import read, rel, blocks240
from steps.s06_compare import is2rate
cfg = yaml.safe_load(open("config.yaml")); W = "data/work/grid30"; D = "data/work/is2rate/max_masked"; os.makedirs(f"{D}/polys", exist_ok=True)
rd = lambda p: np.fromfile(p, dtype=np.float32).reshape(1144, 1225)
with rasterio.open(r"C:\Users\wharf\Documents\Waikato_IS2_raw_data\02_insar\licsbas_073D_2018_22\job_record\73_20180709_20220630_new_params.vel.mskd.geo.tif") as r: m14 = r.read(1)
v = np.where(np.isfinite(m14) & (m14 != 0), rd("data/work/insar/vel") / rd("data/work/insar/U.geo"), np.nan).astype(np.float32)
cop, tf = read(f"{W}/cop30_nzvd2016.tif"); g = np.full(cop.shape, np.nan, np.float32)
reproject(v, g, src_transform=from_origin(174.7164445 - 0.0005, -37.094 + 0.0005, 0.001, 0.001), src_crs="EPSG:4326", dst_transform=tf, dst_crs="EPSG:2193", resampling=Resampling.average, src_nodata=np.nan, dst_nodata=np.nan)
FARM = (read(f"{W}/farmfrac.tif")[0] >= 0.9) & (read(f"{W}/slope_deg.tif")[0] < 3); near = read(f"{W}/near_organic.tif")[0].astype(bool)
M = blocks240(rel(g, FARM, FARM & ~near, cfg)[0], cfg); tf240 = from_origin(tf.c, tf.f, 240, 240); R = {}
P = {"Oct 2018 - Jun 2022": ("2018-10-13", "2022-06-30"), "Jul 2022 - Oct 2026": ("2022-07-01", "2026-10-20"), "2018-2026": ("2018-10-13", "2026-10-20")}
for t in cfg["insar"]["zones_mm_yr"]:
    nm = f"insar_le_{-t}"; m = np.isfinite(M) & (M <= t)
    gpd.GeoDataFrame(geometry=[gpd.GeoSeries([shape(s) for s, k in shapes(m.astype("uint8"), mask=m, transform=tf240) if k == 1]).unary_union], crs=2193).to_file(f"{D}/polys/{nm}.gpkg")
    R[nm] = {"km2": round(m.sum() * 0.0576, 1), "InSAR map median": round(float(np.nanmedian(M[m])), 1)}
    for pk, (t0, t1) in P.items():
        tag = f"{nm}_{t0[:4]}_{t1[:4]}"
        is2rate(cfg, tag, dict(start=t0, end=t1, seasonal="yes", target_classes="30, 40", ref_classes="30, 40", reference="local",
                               ref_area=os.path.abspath("data/work/is2rate/max/polys/study.gpkg").replace("\\", "/"), target_area=os.path.abspath(f"{D}/polys/{nm}.gpkg").replace("\\", "/"),
                               cache_dir=os.path.abspath("data/work/is2rate/max/cache").replace("\\", "/")), D)
        j = json.load(open(f"{D}/{tag}_target.json")); R[nm][pk] = dict(rate=j["rate"], se=j["se_formula"], passes=j["passes"])
    print(nm, R[nm], flush=True)
json.dump(R, open("outputs/tables/max_rates_zones_masked.json", "w"), indent=1)
