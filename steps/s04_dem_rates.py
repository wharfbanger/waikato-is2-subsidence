"""Step 4 (Table S1): LiDAR - Copernicus and repeat LiDAR.
(a) LiDAR - Copernicus: the 2021 LiDAR corrected by its ICESat-2 survey bias (10 km grid, step 3; correcting the 2024 survey too changes nothing after referencing), minus the NZVD2016 Copernicus DEM,
    divided by the interval; the 2024 pair where it exists, else the 2021 pair; each farmland cell referenced locally; 240 m map;
    Organic-Soil farmland value (median of grassland cells) with a 10 km block-bootstrap CI.
Outputs: data/work/grid30/lidcop_rate.tif, lidcop_rel.tif; outputs/maps/lidcop_rel_240m.tif; outputs/tables/lidcop.json."""
import json, os, numpy as np, rasterio
from lib.maps import read, rel, ref_adapt, blocks240, block_id, block_boot, bias_grid


def write(p, a, tf, res=30.0):
    from rasterio.transform import from_origin
    t = from_origin(tf.c, tf.f, res, res)
    with rasterio.open(p, "w", driver="GTiff", height=a.shape[0], width=a.shape[1], count=1, dtype="float32", crs="EPSG:2193", transform=t, nodata=np.nan, compress="deflate") as d:
        d.write(a.astype("float32"), 1)


def lidar_copernicus(cfg, correct_2024=True):
    W, work = f"{cfg['paths']['work']}/grid30", cfg["paths"]["work"]; c0 = cfg["copernicus"]["date"]; G = cfg["grids"]["bias_square_m"]
    cop, tf = read(f"{W}/cop30_nzvd2016.tif"); L21, _ = read(f"{W}/dem30_waikato_2021.tif"); L24, _ = read(f"{W}/dem30_waikato_2024.tif")
    b21 = bias_grid(f"{work}/bias2021_10km.csv", float(open(f"{work}/bias2021_regional.txt").read()), cop.shape, tf, G)
    b24 = bias_grid(f"{work}/bias2024_10km.csv", float(open(f"{work}/bias2024_regional.txt").read()), cop.shape, tf, G) if correct_2024 else 0.0
    d = cfg["lidar_dates"]
    C21 = (L21 + b21 - cop) * 1000 / (d["waikato_2021"] - c0); C24 = (L24 + b24 - cop) * 1000 / (d["waikato_2024"] - c0)
    return np.where(np.isfinite(C24), C24, C21).astype(np.float32), np.where(np.isfinite(C24), 2, np.where(np.isfinite(C21), 1, 0)).astype(np.int8), tf


def run(cfg):
    W, out = f"{cfg['paths']['work']}/grid30", cfg["paths"]["out"]; os.makedirs(f"{out}/maps", exist_ok=True); os.makedirs(f"{out}/tables", exist_ok=True)
    m = cfg["masks"]; slope, _ = read(f"{W}/slope_deg.tif"); grass, _ = read(f"{W}/grassfrac.tif"); farmf, _ = read(f"{W}/farmfrac.tif")
    org = read(f"{W}/organic.tif")[0].astype(bool); near = read(f"{W}/near_organic.tif")[0].astype(bool)
    flat = slope < m["slope_max_deg"]; PAST = (grass >= m["farmland_fraction_min"]) & flat; FARM = (farmf >= m["farmland_fraction_min"]) & flat
    PEAT, NONP, NONPF = PAST & org, PAST & ~near, FARM & ~near
    R = {}
    for tag, c24 in (("2021 survey bias-corrected (Methods; used)", False), ("both surveys bias-corrected (as the v84 code)", True)):
        C, src, tf = lidar_copernicus(cfg, c24); cap = cfg["maps"]["ref_cap_mm_yr"]
        ok = PEAT & np.isfinite(C) & (np.abs(C) < cap); ref = ref_adapt(C, NONP, cfg); okr = ok & np.isfinite(ref); r = (C - ref)[okr]
        R[tag] = dict(cells=int(okr.sum()), median_rel=round(float(np.median(r)), 1),
                      median_rel_CI=block_boot(r, block_id(C.shape, tf, cfg["maps"]["block_bootstrap_km"] * 1000)[okr], cfg["maps"]["n_block_boot"]),
                      unreferenced_reference_median=round(float(np.nanmedian(ref[okr])), 1), cells_2024_pair=int((src == 2).sum()), cells_2021_pair=int((src == 1).sum()))
        print(tag, R[tag], flush=True)
        if not c24:
            write(f"{W}/lidcop_rate.tif", C, tf); relf, _ = rel(C, FARM, NONPF, cfg); write(f"{W}/lidcop_rel.tif", relf, tf)
            M = blocks240(relf, cfg); write(f"{out}/maps/lidcop_rel_240m.tif", M, tf, 240.0)
    json.dump(R, open(f"{out}/tables/lidcop.json", "w"), indent=1)
