"""Step 4 (Table S1): LiDAR - Copernicus and repeat LiDAR.
(a) LiDAR - Copernicus: LiDAR minus the NZVD2016 Copernicus DEM (no survey-offset correction: offsets cancel in the local reference),
    divided by the interval; the 2024 pair where it exists, else the 2021 pair; each farmland cell referenced locally; 240 m map;
    Organic-Soil farmland value (median of farmland cells, farmland reference) with a 10 km block-bootstrap CI; grassland-only as a sensitivity.
Outputs: data/work/grid30/lidcop_rate.tif, lidcop_rel.tif; outputs/maps/lidcop_rel_240m.tif; outputs/tables/lidcop.json."""
import json, os, numpy as np, rasterio
from lib.maps import read, rel, ref_adapt, blocks240, block_id, block_boot, bias_grid


def write(p, a, tf, res=30.0):
    from rasterio.transform import from_origin
    t = from_origin(tf.c, tf.f, res, res)
    with rasterio.open(p, "w", driver="GTiff", height=a.shape[0], width=a.shape[1], count=1, dtype="float32", crs="EPSG:2193", transform=t, nodata=np.nan, compress="deflate") as d:
        d.write(a.astype("float32"), 1)


def lidar_copernicus(cfg, offset_2021_m=0.0):
    """(LiDAR - Copernicus) / interval; the 2024 pair where it exists, else the 2021 pair. offset_2021_m is added to the 2021 LiDAR (sensitivity only)."""
    W = f"{cfg['paths']['work']}/grid30"; c0 = cfg["copernicus"]["date"]; d = cfg["lidar_dates"]
    cop, tf = read(f"{W}/cop30_nzvd2016.tif"); L21, _ = read(f"{W}/dem30_waikato_2021.tif"); L24, _ = read(f"{W}/dem30_waikato_2024.tif")
    C21 = (L21 + offset_2021_m - cop) * 1000 / (d["waikato_2021"] - c0); C24 = (L24 - cop) * 1000 / (d["waikato_2024"] - c0)
    return np.where(np.isfinite(C24), C24, C21).astype(np.float32), np.where(np.isfinite(C24), 2, np.where(np.isfinite(C21), 1, 0)).astype(np.int8), tf


def areas(cfg, cop_shape, tf):
    """S-map Organic Soils in the study rectangle; Table 1 overlaps: LINZ tile footprints of each survey on land in the study rectangle (land: WorldCover not water)."""
    import geopandas as gpd, pyproj
    from shapely.geometry import box
    from rasterio.features import geometry_mask
    from rasterio.warp import reproject, Resampling
    raw = cfg["paths"]["raw"]; a = cfg["study_area"]; rect = box(a["xmin"], a["ymin"], a["xmax"], a["ymax"])
    org = gpd.read_file(f"{raw}/soils/smap_organic_waikato.shp").to_crs(2193).unary_union
    TI = json.load(open(f"{raw}/tile_index.json")); tr = pyproj.Transformer.from_crs(4326, 2193, always_xy=True)
    def foot(col):
        from shapely.ops import transform, unary_union
        return unary_union([transform(tr.transform, box(*b)) for n, b, *_ in TI[col]["tiles"]]).intersection(rect)
    with rasterio.open(f"{raw}/worldcover/ESA_WorldCover_10m_2021_v200_S39E174_Map.tif") as w:
        win = rasterio.windows.from_bounds(174.6, -38.2, 175.9, -36.95, w.transform).round_offsets().round_lengths(); s_ = w.read(1, window=win); wt = w.window_transform(win)
    land = np.zeros(cop_shape, np.float32); reproject(((s_ != 80) & (s_ != 0)).astype(np.float32), land, src_transform=wt, src_crs="EPSG:4326", dst_transform=tf, dst_crs="EPSG:2193", resampling=Resampling.average)
    f24, f21 = foot("waikato_2024"), foot("waikato_2021")
    km2 = lambda g: round(float((land * ~geometry_mask([g], cop_shape, tf)).sum() * 0.0009), 0) if not g.is_empty else 0.0
    return {"land in the study rectangle (km2)": round(float(land.sum() * 0.0009), 0), "S-map Organic Soils in the study rectangle (km2)": round(org.intersection(rect).area / 1e6, 0),
            "Copernicus -> 2024 overlap (km2)": km2(f24), "Copernicus -> 2021 overlap, outside 2024 (km2)": km2(f21.difference(f24))}


def run(cfg):
    W, out = f"{cfg['paths']['work']}/grid30", cfg["paths"]["out"]; os.makedirs(f"{out}/maps", exist_ok=True); os.makedirs(f"{out}/tables", exist_ok=True)
    m = cfg["masks"]; slope, _ = read(f"{W}/slope_deg.tif"); grass, _ = read(f"{W}/grassfrac.tif"); farmf, _ = read(f"{W}/farmfrac.tif")
    org = read(f"{W}/organic.tif")[0].astype(bool); near = read(f"{W}/near_organic.tif")[0].astype(bool)
    flat = slope < m["slope_max_deg"]; PAST = (grass >= m["farmland_fraction_min"]) & flat; FARM = (farmf >= m["farmland_fraction_min"]) & flat
    PEAT, NONP, NONPF = PAST & org, PAST & ~near, FARM & ~near
    R = {}; cap = cfg["maps"]["ref_cap_mm_yr"]; O = cfg["grids"]["square_origin"]; nb = cfg["maps"]["n_block_boot"]; G = cfg["maps"]["block_bootstrap_km"] * 1000
    off = json.load(open(f"{out}/tables/icesat2_rates.json"))["survey bias"]["survey_2021"]["farm_minus_built_same_pass_cm"] / 100
    for tag, o in (("Organic-Soil farmland", 0.0), ("sensitivity: 2021 LiDAR corrected by its same-pass offset", off)):
        C, src, tf = lidar_copernicus(cfg, o)
        ok = FARM & org & np.isfinite(C) & (np.abs(C) < cap); ref = ref_adapt(C, NONPF, cfg); okr = ok & np.isfinite(ref); r = (C - ref)[okr]
        R[tag] = dict(cells=int(okr.sum()), median_rel=round(float(np.median(r)), 1), median_rel_CI=block_boot(r, block_id(C.shape, tf, G, O)[okr], nb),
                      unreferenced_reference_median=round(float(np.nanmedian(ref[okr])), 1), cells_2024_pair=int((src == 2).sum()), cells_2021_pair=int((src == 1).sum()), offset_2021_cm=round(o * 100, 1))
        print(tag, R[tag], flush=True)
        if o == 0.0:
            write(f"{W}/lidcop_rate.tif", C, tf); relf, _ = rel(C, FARM, NONPF, cfg); write(f"{W}/lidcop_rel.tif", relf, tf)
            M = blocks240(relf, cfg); write(f"{out}/maps/lidcop_rel_240m.tif", M, tf, 240.0)
    # Organic-Soil farmland by land cover (farmland reference): grassland >= 90% and cropland >= 50% of the 30 m cell
    from rasterio.warp import reproject, Resampling
    C, _, tf = lidar_copernicus(cfg); ref = ref_adapt(C, NONPF, cfg); ok = np.isfinite(C) & (np.abs(C) < cap) & np.isfinite(ref) & org
    with rasterio.open(f"{cfg['paths']['raw']}/worldcover/ESA_WorldCover_10m_2021_v200_S39E174_Map.tif") as w:
        win = rasterio.windows.from_bounds(174.6, -38.2, 175.9, -36.95, w.transform).round_offsets().round_lengths(); s_ = w.read(1, window=win); wt = w.window_transform(win)
    crop = np.zeros(C.shape, np.float32); reproject((s_ == 40).astype(np.float32), crop, src_transform=wt, src_crs="EPSG:4326", dst_transform=tf, dst_crs="EPSG:2193", resampling=Resampling.average)
    R["Organic-Soil farmland by land cover"] = {k: round(float(np.median((C - ref)[ok & m])), 1) for k, m in (("grassland", PAST), ("cropland", (crop >= 0.5) & FARM))}
    R["areas"] = areas(cfg, cop_shape=C.shape, tf=tf)
    json.dump(R, open(f"{out}/tables/lidcop.json", "w"), indent=1)


# ---------------------------------------------------------------- (b) repeat LiDAR
def repeat_lidar(cfg):
    import pandas as pd, geopandas as gpd
    from rasterio.features import geometry_mask
    from rasterio.warp import reproject, Resampling
    from lib.repeat import diff_grid, capture_dates, cells_to_grid
    raw, work, out = cfg["paths"]["raw"], cfg["paths"]["work"], cfg["paths"]["out"]; rl = cfg["repeat_lidar"]; os.makedirs(f"{work}/repeat", exist_ok=True)
    RS = rl["cell_m"]; dh, tf = diff_grid(cfg, *cfg["lidar"]["repeat_pair"], cfg["lidar"]["resolutions_m"]["repeat"], RS); sh = dh.shape
    t1 = cells_to_grid(capture_dates(cfg, rl["laz_earlier"], f"{work}/repeat/capture_{rl['laz_earlier']}_{RS}m.parquet", RS), sh)
    t2 = cells_to_grid(capture_dates(cfg, rl["laz_later"], f"{work}/repeat/capture_{rl['laz_later']}_{RS}m.parquet", RS), sh)
    dt = t2 - t1; dt_coll = rl["interval_collection_yr"]; dtu = np.where(np.isfinite(dt) & (dt > 4), dt, dt_coll); rate = dh / dtu * 1000
    with rasterio.open(f"{raw}/worldcover/ESA_WorldCover_10m_2021_v200_S39E174_Map.tif") as w:
        win = rasterio.windows.from_bounds(175.2, -37.6, 175.8, -37.0, w.transform).round_offsets().round_lengths(); s_ = w.read(1, window=win); wt = w.window_transform(win)
    fa = np.zeros(sh, np.float32)
    reproject(np.isin(s_, cfg["masks"]["worldcover_farmland"]).astype(np.float32), fa, src_transform=wt, src_crs="EPSG:4326", dst_transform=tf, dst_crs="EPSG:2193", resampling=Resampling.average)
    org = gpd.read_file(f"{raw}/soils/smap_organic_waikato.shp").to_crs(2193); bod = gpd.read_file(f"{raw}/soils/Waikato_peat_bodies_stats.gpkg").to_crs(2193)
    peat = ~geometry_mask(org.geometry, sh, tf); near = ~geometry_mask(org.geometry.buffer(cfg["masks"]["buffer_m"]), sh, tf)
    body = {nm: ~geometry_mask(bod[bod["name"] == nm].geometry, sh, tf) for nm in ("Torehape", "Kopuatai")}
    for nm, a in (("repeat_dh", dh), ("repeat_rate", rate), ("capture_earlier", t1), ("capture_later", t2)): write(f"{work}/repeat/{nm}_{RS}m.tif", a, tf, float(RS))
    R = {"cells_differenced": int(np.isfinite(dh).sum()), "area_km2": round(float(np.isfinite(dh).sum() * (RS / 1000) ** 2), 0)}
    # (1) same campaign: Organic-Soil farmland minus mineral-soil farmland flown at the same time; untied, no ICESat-2
    ok = np.isfinite(dh) & (fa > rl["farm_fraction_min"]) & (np.abs(dh) < 1) & np.isfinite(dt) & (dt > 4); P, MIN = ok & peat, ok & ~near
    CAP = {"Feb 2015": t1 < 2016, "Dec 2017 - Mar 2018": t1 > 2017.5}
    SITES = {"northern Torehape": (P & body["Torehape"], "Feb 2015"), "southern Torehape": (P & body["Torehape"], "Dec 2017 - Mar 2018"),
             "drained margins of Kopuatai": (P & body["Kopuatai"], "Dec 2017 - Mar 2018")}
    # flights of the earlier survey (capture dates): each has its own height offset, so each site cell is referenced to mineral-soil farmland of its own flight
    # (v96; the 2017-18 campaign pooled the December 2017 and March 2018 flights, whose mineral-soil farmland differs by 3.3 mm/yr)
    FL = {"Feb 2015": t1 < 2016, "Dec 2017": (t1 > 2017.5) & (t1 < rl["flight_split_yr"]), "Mar 2018": t1 >= rl["flight_split_yr"]}
    REFR = {f: float(np.median(rate[MIN & c])) for f, c in FL.items()}; REFD = {f: float(np.median(dh[MIN & c])) for f, c in FL.items()}
    rref, dref = np.full(sh, np.nan), np.full(sh, np.nan)
    for f, c in FL.items(): rref[c] = REFR[f]; dref[c] = REFD[f]
    R["same campaign"] = {}
    for nm, (m, cp) in SITES.items():
        tm = m & CAP[cp]
        R["same campaign"][nm] = dict(campaign=cp, dh_cm=round(float(np.median(dh[tm] - dref[tm])) * 100, 1), interval_yr=round(float(np.median(dt[tm])), 1),
                                      rate=round(float(np.median(rate[tm] - rref[tm])), 1), target_cells=int(tm.sum()),
                                      cells_by_flight={f: int((tm & c).sum()) for f, c in FL.items() if (tm & c).sum()},
                                      rate_by_flight={f: round(float(np.median(rate[tm & c]) - REFR[f]), 1) for f, c in FL.items() if (tm & c).sum() >= 20},
                                      reference_cells={f: int((MIN & c).sum()) for f, c in FL.items() if (tm & c).sum()})
    # site masks (Table 4 rows): written on the 30 m grid for steps 5 (InSAR) and here (LiDAR - Copernicus)
    if RS == 30:
        code = np.zeros(sh, np.uint8)
        for i, (nm, (m, cp)) in enumerate(SITES.items(), 1): code[m & CAP[cp]] = i
        with rasterio.open(f"{work}/repeat/sites_30m.tif", "w", driver="GTiff", height=sh[0], width=sh[1], count=1, dtype="uint8", crs="EPSG:2193", transform=tf, nodata=0, compress="deflate") as d_: d_.write(code, 1)
        lr = read(f"{cfg['paths']['work']}/grid30/lidcop_rel.tif")[0]
        for i, nm in enumerate(SITES, 1): R["same campaign"][nm]["LiDAR - Copernicus median"] = round(float(np.nanmedian(lr[code == i])), 1)
    R["mineral-soil farmland, untied"] = {cp: round(float(np.median(rate[MIN & c])), 1) for cp, c in CAP.items()}
    R["mineral-soil farmland, untied, by flight"] = {f: round(v, 1) for f, v in REFR.items()}
    # v97: capture dates of each flight (farmland cells used; 1st-99th percentile), cells between flights, and farmland cells without a date
    def _d(y): return str((np.datetime64("1970-01-01T00:00:00", "s") + np.timedelta64(int((float(y) - 1970) * 365.2425 * 86400 + 13 * 3600), "s")).astype("datetime64[D]"))   # inverse of lib/repeat.py decimal year; NZ daylight time (UTC+13)
    R["flights, capture dates of farmland cells (1st-99th percentile, NZDT)"] = {f: [_d(np.percentile(t1[ok & c], 1)), _d(np.percentile(t1[ok & c], 99))] for f, c in FL.items()}
    FA = np.isfinite(dh) & (fa > rl["farm_fraction_min"]) & (np.abs(dh) < 1)
    R["farmland cells between flights (2016 - mid 2017; not used)"] = int((FA & (t1 >= 2016) & (t1 <= 2017.5)).sum())
    R["farmland cells without a capture date (interval of the collections used; not in same-flight rates)"] = {"cells": int((FA & ~np.isfinite(dt)).sum()), "share": round(float((FA & ~np.isfinite(dt)).sum() / FA.sum()), 3)}
    # 10 km squares (bootstrap grid) holding the cells of the bog and of each site (Table 4 note)
    YYc, XXc = np.mgrid[0:sh[0], 0:sh[1]]; Xc = tf.c + (XXc + .5) * RS; Yc = tf.f - (YYc + .5) * RS; bx0, by0 = cfg["grids"]["square_origin"]; Gb = cfg["maps"]["block_bootstrap_km"] * 1000
    sq = ((Xc - bx0) // Gb).astype(np.int64) * 100000 + ((Yc - by0) // Gb).astype(np.int64)
    def nsq(m):
        u, c = np.unique(sq[m], return_counts=True); c = np.sort(c)[::-1]; return dict(squares=int(u.size), share_of_two_largest=round(float(c[:2].sum() / c.sum()), 2))
    R["10 km squares"] = {nm: nsq(m & CAP[cp]) for nm, (m, cp) in SITES.items()}
    # (3) tie to ICESat-2 on a square grid; class rates are medians over 16 grid placements
    TIE = json.load(open(f"{out}/tables/icesat2_rates.json"))["absolute"]["mineral-soil farmland"]["rate"]
    npk = (fa > rl["farm_fraction_min"]) & ~near & np.isfinite(rate)
    YY, XX = np.mgrid[0:sh[0], 0:sh[1]]; X = tf.c + (XX + .5) * RS; Y = tf.f - (YY + .5) * RS
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
        r_ = ((tf.f - d.Y.values) // RS).astype(int); c_ = ((d.X.values - tf.c) // RS).astype(int); ok_ = (r_ >= 0) & (r_ < sh[0]) & (c_ >= 0) & (c_ < sh[1])
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
    write(f"{work}/repeat/repeat_tied_{RS}m.tif", t5.astype(np.float32), tf, float(RS))
    # (2) map: rates on the 30 m grid (nearest; identical when cell_m is 30), each farmland cell referenced locally, 240 m grid (Figure 3b)
    from lib.maps import rel, blocks240
    W = f"{work}/grid30"; cop, tf30 = read(f"{W}/cop30_nzvd2016.tif"); r30 = np.full(cop.shape, np.nan, np.float32)
    reproject(rate.astype(np.float32), r30, src_transform=tf, src_crs="EPSG:2193", dst_transform=tf30, dst_crs="EPSG:2193", resampling=Resampling.nearest, src_nodata=np.nan, dst_nodata=np.nan)
    FARM = (read(f"{W}/farmfrac.tif")[0] >= cfg["masks"]["farmland_fraction_min"]) & (read(f"{W}/slope_deg.tif")[0] < cfg["masks"]["slope_max_deg"]); nr = read(f"{W}/near_organic.tif")[0].astype(bool)
    rr_, _ = rel(r30, FARM, FARM & ~nr, cfg); write(f"{W}/repeat_rel.tif", rr_, tf30); write(f"{out}/maps/repeat_rel_240m.tif", blocks240(rr_, cfg), tf30, 240.0)
    json.dump(R, open(f"{out}/tables/repeat_lidar.json", "w"), indent=1); print(json.dumps(R, indent=1), flush=True)


_run_a = run
def run(cfg):
    _run_a(cfg); repeat_lidar(cfg)
