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


# ---------------------------------------------------------------- (b) repeat LiDAR
def repeat_lidar(cfg):
    import pandas as pd, geopandas as gpd
    from rasterio.features import geometry_mask
    from rasterio.warp import reproject, Resampling
    from lib.repeat import diff96, capture_dates, cells_to_grid
    raw, work, out = cfg["paths"]["raw"], cfg["paths"]["work"], cfg["paths"]["out"]; rl = cfg["repeat_lidar"]; os.makedirs(f"{work}/repeat", exist_ok=True)
    dh, tf = diff96(cfg, *cfg["lidar"]["repeat_pair"]); sh = dh.shape
    t1 = cells_to_grid(capture_dates(cfg, rl["laz_earlier"], f"{work}/repeat/capture_{rl['laz_earlier']}.parquet"), sh)
    t2 = cells_to_grid(capture_dates(cfg, rl["laz_later"], f"{work}/repeat/capture_{rl['laz_later']}.parquet"), sh)
    dt = t2 - t1; dt_coll = rl["interval_collection_yr"]; dtu = np.where(np.isfinite(dt) & (dt > 4), dt, dt_coll); rate = dh / dtu * 1000
    with rasterio.open(f"{raw}/worldcover/ESA_WorldCover_10m_2021_v200_S39E174_Map.tif") as w:
        win = rasterio.windows.from_bounds(175.2, -37.6, 175.8, -37.0, w.transform).round_offsets().round_lengths(); s_ = w.read(1, window=win); wt = w.window_transform(win)
    fa = np.zeros(sh, np.float32)
    reproject(np.isin(s_, cfg["masks"]["worldcover_farmland"]).astype(np.float32), fa, src_transform=wt, src_crs="EPSG:4326", dst_transform=tf, dst_crs="EPSG:2193", resampling=Resampling.average)
    org = gpd.read_file(f"{raw}/soils/smap_organic_waikato.shp").to_crs(2193); bod = gpd.read_file(f"{raw}/soils/Waikato_peat_bodies_stats.gpkg").to_crs(2193)
    peat = ~geometry_mask(org.geometry, sh, tf); near = ~geometry_mask(org.geometry.buffer(cfg["masks"]["buffer_m"]), sh, tf)
    body = {nm: ~geometry_mask(bod[bod["name"] == nm].geometry, sh, tf) for nm in ("Torehape", "Kopuatai")}
    for nm, a in (("repeat_dh_96m", dh), ("repeat_rate_96m", rate), ("capture_earlier_96m", t1), ("capture_later_96m", t2)): write(f"{work}/repeat/{nm}.tif", a, tf, 96.0)
    R = {"cells_differenced": int(np.isfinite(dh).sum()), "area_km2": round(float(np.isfinite(dh).sum() * 0.096 ** 2), 0)}
    # (1) same campaign: Organic-Soil farmland minus mineral-soil farmland flown at the same time; untied, no ICESat-2
    ok = np.isfinite(dh) & (fa > rl["farm_fraction_min"]) & (np.abs(dh) < 1) & np.isfinite(dt) & (dt > 4); P, MIN = ok & peat, ok & ~near
    CAP = {"Feb 2015": t1 < 2016, "Dec 2017 - Mar 2018": t1 > 2017.5}
    SITES = {"northern Torehape": (P & body["Torehape"], "Feb 2015"), "southern Torehape": (P & body["Torehape"], "Dec 2017 - Mar 2018"),
             "drained margins of Kopuatai": (P & body["Kopuatai"], "Dec 2017 - Mar 2018")}
    R["same campaign"] = {}
    for nm, (m, cp) in SITES.items():
        tm, rm = m & CAP[cp], MIN & CAP[cp]
        R["same campaign"][nm] = dict(campaign=cp, dh_cm=round(float(np.median(dh[tm]) - np.median(dh[rm])) * 100, 1), interval_yr=round(float(np.median(dt[tm])), 1),
                                      rate=round(float(np.median(rate[tm]) - np.median(rate[rm])), 1), target_cells=int(tm.sum()), reference_cells=int(rm.sum()))
    R["mineral-soil farmland, untied"] = {cp: round(float(np.median(rate[MIN & c])), 1) for cp, c in CAP.items()}
    # (3) tie to ICESat-2 on a square grid; class rates are medians over 16 grid placements
    TIE = json.load(open(f"{out}/tables/icesat2_rates.json"))["absolute"]["mineral-soil farmland"]["rate"]
    npk = (fa > rl["farm_fraction_min"]) & ~near & np.isfinite(rate)
    YY, XX = np.mgrid[0:sh[0], 0:sh[1]]; X = tf.c + (XX + .5) * 96; Y = tf.f - (YY + .5) * 96
    def tie(G, x0, y0):
        bid = ((X - x0) // G).astype(np.int64) * 100000 + ((Y - y0) // G).astype(np.int64); adj = np.full(sh, np.nan); n = 0
        for b in np.unique(bid[npk]):
            m = (bid == b) & npk
            if m.sum() >= rl["tie_min_cells"]: adj[bid == b] = np.median(rate[m]) - TIE; n += 1
        adj[~np.isfinite(adj)] = np.median(rate[npk]) - TIE; return rate - adj, n
    S = pd.read_parquet(f"{work}/segments_classified.parquet").drop_duplicates(subset=["X", "Y"]); a = cfg["study_area"]
    F = S[S.land == "farmland"]; MINs = F[(F.order != "Organic") & (F.order != "not mapped") & ~F.near_peat]; hx0, hy0, hx1, hy1 = a["hauraki_box"]
    CL = {"Organic-Soil farmland": F[F.order == "Organic"], "mineral-soil farmland": MINs, "Hauraki Plains mineral-soil farmland": MINs[MINs.X.between(hx0, hx1) & MINs.Y.between(hy0, hy1)],
          "Kopuatai bog": S[(S.land == "natural wetland") & (S.order == "Organic") & (S.body == "Kopuatai")]}
    def samp(t, d):
        r_ = ((tf.f - d.Y.values) // 96).astype(int); c_ = ((d.X.values - tf.c) // 96).astype(int); ok_ = (r_ >= 0) & (r_ < sh[0]) & (c_ >= 0) & (c_ < sh[1])
        v = np.full(len(d), np.nan); v[ok_] = t[r_[ok_], c_[ok_]]; v = v[np.isfinite(v)]; return round(float(np.median(v)), 1), int(v.size)
    ox, oy = cfg["grids"]["square_origin"]; t5, nsq = tie(cfg["grids"]["tie_square_m"], ox, oy)
    R["tie"] = dict(ICESat2_mineral_rate=TIE, squares_with_own_offset=nsq, default_placement={k: dict(zip(("median", "points"), samp(t5, d))) for k, d in CL.items()}, spread={})
    for G in (5000, 10000):
        v = {k: [] for k in CL}
        for dx in np.arange(0, G, G / 4):
            for dy in np.arange(0, G, G / 4):
                t, _ = tie(G, ox + dx, oy + dy)
                for k, d in CL.items(): v[k].append(samp(t, d)[0])
        R["tie"]["spread"][f"{G // 1000} km"] = {k: dict(min=min(x), median=float(np.median(x)), max=max(x)) for k, x in v.items()}
    # tied class rates = median over the 16 placements of the 5 km grid (the default placement alone lies near the fast end of the spread)
    R["tie"]["classes"] = {k: dict(median=round(v["median"], 1), points=R["tie"]["default_placement"][k]["points"]) for k, v in R["tie"]["spread"][f"{cfg['grids']['tie_square_m'] // 1000} km"].items()}
    write(f"{work}/repeat/repeat_tied_96m.tif", t5.astype(np.float32), tf, 96.0)
    # (2) map: 96 m rates on the 30 m grid (nearest), each farmland cell referenced locally, 240 m grid (Figure 3b)
    from lib.maps import rel, blocks240
    W = f"{work}/grid30"; cop, tf30 = read(f"{W}/cop30_nzvd2016.tif"); r30 = np.full(cop.shape, np.nan, np.float32)
    reproject(rate.astype(np.float32), r30, src_transform=tf, src_crs="EPSG:2193", dst_transform=tf30, dst_crs="EPSG:2193", resampling=Resampling.nearest, src_nodata=np.nan, dst_nodata=np.nan)
    FARM = (read(f"{W}/farmfrac.tif")[0] >= cfg["masks"]["farmland_fraction_min"]) & (read(f"{W}/slope_deg.tif")[0] < cfg["masks"]["slope_max_deg"]); nr = read(f"{W}/near_organic.tif")[0].astype(bool)
    rr_, _ = rel(r30, FARM, FARM & ~nr, cfg); write(f"{W}/repeat_rel.tif", rr_, tf30); write(f"{out}/maps/repeat_rel_240m.tif", blocks240(rr_, cfg), tf30, 240.0)
    json.dump(R, open(f"{out}/tables/repeat_lidar.json", "w"), indent=1); print(json.dumps(R, indent=1), flush=True)


_run_a = run
def run(cfg):
    _run_a(cfg); repeat_lidar(cfg)
