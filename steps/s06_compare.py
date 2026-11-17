"""Step 6 (Table S1): hexagons, maximum rates and sensitivity.
(a) ICESat-2 hexagon rates with is2rate (3 km, regional seasonal cycle), from the step-2 SlideRule download and the step-1 LiDAR mosaics;
(b) LiDAR - Copernicus and InSAR on the same hexagons: median of at least 30 referenced 30 m grassland cells; hexagons at least 25% farmland;
    Pearson r with the ICESat-2 hexagon rates;
(c) maximum rates with is2rate: distance bands inside Organic Soils and InSAR zones, three periods (Table 4).
Outputs: outputs/maps/hex_3km.gpkg, outputs/tables/hexagons.json, outputs/tables/max_rates.json."""
import glob, json, os, subprocess, sys, numpy as np, pandas as pd, geopandas as gpd, rasterio
from rasterio.features import rasterize
from lib.maps import read, local_ref


def is2rate(cfg, name, extra, out_dir):
    """Runs is2rate with the pipeline's inputs; returns its output folder."""
    i, w, raw = cfg["is2rate"], cfg["paths"]["work"], cfg["paths"]["raw"]; os.makedirs(out_dir, exist_ok=True)
    P = dict(bbox=", ".join(str(x) for x in i["bbox"]), start=i["start"], end=i["end"], crs=2193,
             dem=", ".join(os.path.abspath(f"{w}/vrt/{c}.vrt").replace("\\", "/") for c in cfg["lidar"]["waikato"]),
             dem_geoid=os.path.abspath(f"{raw}/proj/nz_linz_nzgeoid2016.tif").replace("\\", "/"), landcover="worldcover",
             exclude=os.path.abspath(f"{raw}/soils/smap_organic_waikato.shp").replace("\\", "/"), exclude_buffer_m=cfg["masks"]["buffer_m"],
             segments_file=os.path.abspath(f"{raw}/icesat2").replace("\\", "/") + "/Waikato_??_atl06sr.parquet", out_dir=os.path.abspath(out_dir).replace("\\", "/"), name=name)
    P.update(extra); pf = f"{out_dir}/{name}_params.txt"; open(pf, "w").write("\n".join(f"{k} = {v}" for k, v in P.items()) + "\n")
    r = subprocess.run([sys.executable, i["script"], pf], capture_output=True, text=True); open(f"{out_dir}/{name}.log", "w").write(r.stdout + r.stderr)
    if r.returncode: raise RuntimeError(f"is2rate {name} failed; see {out_dir}/{name}.log")
    return out_dir


def hexagons(cfg):
    W, out = f"{cfg['paths']['work']}/grid30", cfg["paths"]["out"]; h = cfg["hexagons"]; m = cfg["masks"]
    d = is2rate(cfg, "hex_3km", dict(hex_width_m=cfg["grids"]["hex_width_m"], grid_origin=", ".join(str(x) for x in cfg["grids"]["hex_origin"]), seasonal="regional"),
                f"{cfg['paths']['work']}/is2rate/hex_3km")
    hx = gpd.read_file(f"{d}/hex_3km_hexagons.gpkg"); nh = int(hx.hid.max()) + 1
    slope, tf = read(f"{W}/slope_deg.tif"); grass = read(f"{W}/grassfrac.tif")[0]; farmf = read(f"{W}/farmfrac.tif")[0]; near = read(f"{W}/near_organic.tif")[0].astype(bool)
    flat = slope < m["slope_max_deg"]; PAST = (grass >= m["farmland_fraction_min"]) & flat; FARM = (farmf >= m["farmland_fraction_min"]) & flat; NONP = PAST & ~near
    hid = rasterize([(g, int(k)) for g, k in zip(hx.geometry, hx.hid)], out_shape=slope.shape, transform=tf, fill=-1, dtype="int32")
    inn = hid >= 0; nall = np.bincount(hid[inn], minlength=nh); nf = np.bincount(hid[inn & FARM], minlength=nh); hx["farm_share"] = (nf / np.maximum(nall, 1))[hx.hid.values]
    a = cfg["study_area"]; keep = (hx.farm_share >= h["min_farm_share"]) & hx.x.between(a["xmin"], a["xmax"]) & hx.y.between(a["ymin"], a["ymax"]);   # centre in the study rectangle
    hx.loc[~keep, ["rate", "se"]] = np.nan; cap = cfg["maps"]["ref_cap_mm_yr"]; R = {}
    for nm, f in (("LiDAR - Copernicus", "lidcop_rate.tif"), ("InSAR", "insar_up.tif")):
        a = read(f"{W}/{f}")[0]; v = np.where(np.isfinite(a) & (np.abs(a) < cap), a, np.nan)
        ref = local_ref(v, NONP, cfg["maps"]["ref_window_km"][0] * 1000, cap, cfg["maps"]["ref_min_cells"]); P3 = PAST & np.isfinite(v) & np.isfinite(ref)
        g = pd.DataFrame(dict(hid=hid[P3], v=(v - ref)[P3])).groupby("hid").v.agg(["median", "size"]); g = g[(g["size"] >= h["min_cells"]) & g.index.isin(hx.hid[keep])]
        col = "dem_rate" if nm.startswith("LiDAR") else "insar_rate"; hx[col] = hx.hid.map(g["median"])
        both = hx[col].notna() & hx.rate.notna()
        R[nm] = dict(hexagons=int(hx[col].notna().sum()), with_ICESat2=int(both.sum()), r=round(float(np.corrcoef(hx[col][both], hx.rate[both])[0, 1]), 2))
    R["ICESat-2"] = dict(hexagons=int(hx.rate.notna().sum()), median_passes=float(hx.loc[hx.rate.notna(), "n_pass"].median()) if "n_pass" in hx else None)
    hx.to_file(f"{out}/maps/hex_3km.gpkg"); json.dump(R, open(f"{out}/tables/hexagons.json", "w"), indent=1); print(json.dumps(R, indent=1), flush=True)


def run(cfg):
    os.makedirs(f"{cfg['paths']['out']}/maps", exist_ok=True); hexagons(cfg)


def max_rates(cfg):
    """Table 4: ICESat-2 rate of each subset (one target, regional mineral-soil reference on the same passes) per period, with is2rate;
    InSAR map median of each zone; LiDAR - Copernicus median of the referenced 30 m grassland cells of each band."""
    from rasterio.features import shapes, geometry_mask
    from shapely.geometry import shape
    from lib.maps import ref_adapt
    work, out, raw = cfg["paths"]["work"], cfg["paths"]["out"], cfg["paths"]["raw"]; D = f"{work}/is2rate/max"; os.makedirs(f"{D}/polys", exist_ok=True); a = cfg["study_area"]
    org = gpd.read_file(f"{raw}/soils/smap_organic_waikato.shp").to_crs(2193).unary_union.buffer(0); T = {}
    E = cfg["max_rates"]["distance_bands_m"]
    for lo, hi in zip(E[:-1], E[1:]):
        g = org.buffer(-lo).difference(org.buffer(-hi)) if lo > 0 else org.difference(org.buffer(-hi)); T[f"band_{lo}_{hi}"] = g
    with rasterio.open(f"{out}/maps/insar_rel_240m.tif") as d: I = d.read(1); tf = d.transform
    for t in cfg["insar"]["zones_mm_yr"]:
        m = np.isfinite(I) & (I <= t); T[f"insar_le_{-t}"] = gpd.GeoSeries([shape(s) for s, v in shapes(m.astype("uint8"), mask=m, transform=tf) if v == 1], crs=2193).unary_union
    study = gpd.GeoDataFrame(geometry=[__import__("shapely").geometry.box(a["xmin"], a["ymin"], a["xmax"], a["ymax"])], crs=2193); study.to_file(f"{D}/polys/study.gpkg")
    P = {"Oct 2018 - Jun 2022": ("2018-10-13", "2022-06-30"), "Jul 2022 - Oct 2026": ("2022-07-01", "2026-10-20"), "2018-2026": ("2018-10-13", "2026-10-20")}
    R = {}
    for nm, g in T.items():
        gpd.GeoDataFrame(geometry=[g], crs=2193).to_file(f"{D}/polys/{nm}.gpkg"); R[nm] = {"km2": round(g.area / 1e6, 1)}
        for pk, (t0, t1) in P.items():
            tag = f"{nm}_{t0[:4]}_{t1[:4]}"
            is2rate(cfg, tag, dict(start=t0, end=t1, seasonal="yes", target_classes="30, 40", ref_classes="30, 40", reference="local",
                                   ref_area=os.path.abspath(f"{D}/polys/study.gpkg").replace("\\", "/"), target_area=os.path.abspath(f"{D}/polys/{nm}.gpkg").replace("\\", "/"),
                                   cache_dir=os.path.abspath(f"{D}/cache").replace("\\", "/")), D)
            j = json.load(open(f"{D}/{tag}_target.json")); R[nm][pk] = dict(rate=j.get("rate"), se=j.get("se_formula", j.get("se")), passes=j.get("n_pass", j.get("passes")))
        print(nm, R[nm], flush=True)
    W = f"{work}/grid30"; m = cfg["masks"]; slope, tf30 = read(f"{W}/slope_deg.tif"); grass = read(f"{W}/grassfrac.tif")[0]
    near = read(f"{W}/near_organic.tif")[0].astype(bool); orgm = read(f"{W}/organic.tif")[0].astype(bool)
    PAST = (grass >= m["farmland_fraction_min"]) & (slope < m["slope_max_deg"]); C = read(f"{W}/lidcop_rate.tif")[0]; cap = cfg["maps"]["ref_cap_mm_yr"]
    ref = ref_adapt(C, PAST & ~near, cfg); ok = PAST & orgm & np.isfinite(C) & (np.abs(C) < cap) & np.isfinite(ref); relc = C - ref
    for nm, g in T.items():
        if nm.startswith("band"):
            mk = ~geometry_mask([g], C.shape, tf30); R[nm]["LiDAR - Copernicus median"] = round(float(np.median(relc[mk & ok])), 1)
        else: R[nm]["InSAR map median"] = round(float(np.nanmedian(I[np.isfinite(I) & (I <= -int(nm.split('_')[-1]))])), 1)
    json.dump(R, open(f"{out}/tables/max_rates.json", "w"), indent=1)


_run_hex = run
def run(cfg):
    _run_hex(cfg); max_rates(cfg)


def repeat_vs_insar(cfg):
    """Section 4.2: Hauraki Plains mineral-soil farmland 240 m cells grouped into quarters by InSAR rate; group medians of InSAR and repeat LiDAR."""
    out = cfg["paths"]["out"]; a = cfg["study_area"]
    I, tf = read(f"{out}/maps/insar_rel_240m.tif"); L = read(f"{out}/maps/repeat_rel_240m.tif")[0]
    nr = read(f"{cfg['paths']['work']}/grid30/near_organic.tif")[0].astype(float); ny, nx = I.shape
    near240 = nr[:ny * 8, :nx * 8].reshape(ny, 8, nx, 8).mean(axis=(1, 3))
    YY, XX = np.mgrid[0:ny, 0:nx]; X = tf.c + (XX + .5) * 240; Y = tf.f - (YY + .5) * 240; hx0, hy0, hx1, hy1 = a["hauraki_box"]
    m = np.isfinite(I) & np.isfinite(L) & (near240 == 0) & (X >= hx0) & (X < hx1) & (Y >= hy0) & (Y < hy1)
    q = pd.qcut(I[m], 4, labels=False); R = dict(cells=int(m.sum()), quarters=[dict(InSAR=round(float(np.median(I[m][q == k])), 1), repeat_LiDAR=round(float(np.median(L[m][q == k])), 1)) for k in range(4)])
    print("repeat vs InSAR", R, flush=True); return R


_run_max = run
def run(cfg):
    _run_max(cfg); json.dump(repeat_vs_insar(cfg), open(f"{cfg['paths']['out']}/tables/repeat_vs_insar.json", "w"), indent=1)


def sensitivity(cfg):
    """Section 5.3: (1) buffer-zone width for the ICESat-2 Organic-Soil farmland rate; (2) reference cap for LiDAR - Copernicus."""
    from scipy.ndimage import distance_transform_edt
    from lib.is2 import diff_series, fit
    from lib.maps import ref_adapt
    work, out = cfg["paths"]["work"], cfg["paths"]["out"]; W = f"{work}/grid30"; rc = cfg["icesat2"]["rates"]; R = {"buffer_m": {}, "reference_cap_mm_yr": {}}
    org, tf = read(f"{W}/organic.tif"); dist = distance_transform_edt(org == 0) * 30.0   # distance to Organic Soils on the 30 m grid
    S = pd.read_parquet(f"{work}/segments_classified.parquet"); r_ = ((tf.f - S.Y.values) // 30).astype(int); c_ = ((S.X.values - tf.c) // 30).astype(int)
    ok = (r_ >= 0) & (r_ < dist.shape[0]) & (c_ >= 0) & (c_ < dist.shape[1]); S["dist"] = np.nan; S.loc[ok, "dist"] = dist[r_[ok], c_[ok]]
    F = S[S.land == "farmland"]; ORG = F[F.order == "Organic"]
    for b in cfg["sensitivity"]["buffers_m"]:
        MIN = F[(F.order != "Organic") & (F.order != "not mapped") & (F.dist > b)] if b != cfg["masks"]["buffer_m"] else F[(F.order != "Organic") & (F.order != "not mapped") & ~F.near_peat]
        f = fit(diff_series(ORG, MIN, rc["min_segments_per_pass"]), n_boot=0); R["buffer_m"][b] = dict(rate=f["rate"], se=f["se"], passes=f["passes"])
    m = cfg["masks"]; slope = read(f"{W}/slope_deg.tif")[0]; grass = read(f"{W}/grassfrac.tif")[0]; near = read(f"{W}/near_organic.tif")[0].astype(bool)
    PAST = (grass >= m["farmland_fraction_min"]) & (slope < m["slope_max_deg"]); C = read(f"{W}/lidcop_rate.tif")[0]; PEAT = PAST & (org > 0)
    for cap in cfg["sensitivity"]["ref_caps_mm_yr"]:
        c2 = dict(cfg, maps=dict(cfg["maps"], ref_cap_mm_yr=cap if cap else 1e9)); ref = ref_adapt(C, PAST & ~near, c2)
        okc = PEAT & np.isfinite(C) & (np.abs(C) < (cap if cap else 1e9)) & np.isfinite(ref); R["reference_cap_mm_yr"][str(cap)] = round(float(np.median((C - ref)[okc])), 1)
    json.dump(R, open(f"{out}/tables/sensitivity.json", "w"), indent=1); print("sensitivity", R, flush=True); return R


_run_cmp = run
def run(cfg):
    _run_cmp(cfg); sensitivity(cfg)
