"""Step 6 (Table S1): hexagons, maximum rates and sensitivity.
(a) ICESat-2 hexagon rates with the pipeline's own ICESat-2 code (lib/is2.py, as step 3): 3 km hexagons, farmland target, mineral-soil farmland
    reference within +-10 km of the centre (hexagon excluded) on the same passes, one regional annual cycle per hexagon class removed before the fit;
    v120: was is2rate (user decision: one ICESat-2 code path);
(b) LiDAR - Copernicus and InSAR on the same hexagons: median of at least 30 referenced 30 m farmland cells; hexagons at least 25% farmland;
    Pearson r with the ICESat-2 hexagon rates;
(c) maximum rates: InSAR zones, three periods, with their distance inside Organic Soils and LiDAR - Copernicus median (Table 5).
Outputs: outputs/maps/hex_3km.gpkg, outputs/tables/hexagons.json, outputs/tables/max_rates.json."""
import glob, json, os, subprocess, sys, numpy as np, pandas as pd, geopandas as gpd, rasterio
from rasterio.features import rasterize
from lib.maps import read, local_ref


def hex_grid(cfg):
    """Flat-topped hexagons of width hex_width_m (flat to flat) from hex_origin, covering the study rectangle plus one hexagon."""
    from shapely.geometry import Polygon
    a, g = cfg["study_area"], cfg["grids"]; w = g["hex_width_m"]; R_ = w / np.sqrt(3); x0, y0 = g["hex_origin"]; dx = 1.5 * R_
    rows = []
    for i in range(int(np.floor((a["xmin"] - w - x0) / dx)), int(np.ceil((a["xmax"] + w - x0) / dx)) + 1):
        xc = x0 + i * dx; off = (w / 2) * (i % 2)
        for k in range(int(np.floor((a["ymin"] - w - y0 - off) / w)), int(np.ceil((a["ymax"] + w - y0 - off) / w)) + 1):
            yc = y0 + off + k * w
            rows.append(dict(x=xc, y=yc, geometry=Polygon([(xc + R_ * np.cos(np.pi / 3 * m), yc + R_ * np.sin(np.pi / 3 * m)) for m in range(6)])))
    hx = gpd.GeoDataFrame(rows, crs=2193); hx.insert(0, "hid", np.arange(len(hx))); return hx


def regional_cycle(series):
    """One annual cycle (sin, cos; m) shared by several per-pass series, each with its own trend and one constant per DEM (pooled least squares)."""
    from lib.is2 import design
    blocks = [design(pp, seasonal=False)[0] for pp in series]; ncol = sum(b.shape[1] for b in blocks); nr = sum(len(pp) for pp in series)
    A = np.zeros((nr, ncol + 2)); y = np.zeros(nr); r0 = c0 = 0
    for pp, B in zip(series, blocks):
        n, m = B.shape; A[r0:r0 + n, c0:c0 + m] = B; t = pp["t"].values; A[r0:r0 + n, -2] = np.sin(2 * np.pi * t); A[r0:r0 + n, -1] = np.cos(2 * np.pi * t)
        y[r0:r0 + n] = pp["r"].values; r0 += n; c0 += m
    return np.linalg.lstsq(A, y, rcond=None)[0][-2:]


def hex_rates(cfg, F, MIN, hx, t1=None):
    """ICESat-2 rate of each hexagon: farmland segments in it minus mineral-soil farmland within +-reference_km of its centre (outside the
    hexagon) on the same pass and DEM; per class (hexagons with at least half their target segments on Organic Soils, and the rest) one annual
    cycle is fitted to all series together and removed; then a trend with one constant per DEM. SE: bootstrap of whole passes (and formula)."""
    from lib.is2 import diff_series, fit
    from scipy.spatial import cKDTree
    h = cfg["hexagons"]; ms = cfg["icesat2"]["rates"]["min_segments_per_pass"]; Rr = h["reference_km"] * 1000.0
    if t1 is not None: F, MIN = F[F.t < t1], MIN[MIN.t < t1]
    tree = cKDTree(np.column_stack([hx.x.values, hx.y.values]))   # nearest centre = containing cell of a regular hexagon grid
    F = F.assign(hid=hx.hid.values[tree.query(np.column_stack([F.X.values, F.Y.values]))[1]])
    MIN = MIN.assign(hid=hx.hid.values[tree.query(np.column_stack([MIN.X.values, MIN.Y.values]))[1]])
    mx, my, mh = MIN.X.values, MIN.Y.values, MIN.hid.values; xc_, yc_ = dict(zip(hx.hid, hx.x)), dict(zip(hx.hid, hx.y)); SER = {}; frac = {}
    for hid_, tg in F.groupby("hid"):
        xc, yc = xc_[hid_], yc_[hid_]
        ref = MIN[(np.abs(mx - xc) <= Rr) & (np.abs(my - yc) <= Rr) & (mh != hid_)]
        ser = diff_series(tg, ref, ms)
        if ser["pass"].nunique() >= 8: SER[hid_] = ser; frac[hid_] = float((tg.order == "Organic").mean())
    cyc = {}
    for c in (True, False):
        hs = [k for k in SER if (frac[k] >= 0.5) == c]
        if len(hs) >= 5: cyc[c] = regional_cycle([SER[k] for k in hs])
    rows = []
    for k, ser in SER.items():
        a_, b_ = cyc.get(frac[k] >= 0.5, (0.0, 0.0)); s2 = ser.assign(r=ser.r - a_ * np.sin(2 * np.pi * ser.t) - b_ * np.cos(2 * np.pi * ser.t))
        f_ = fit(s2, n_boot=h["n_boot"], seasonal=False, seed=int(k))
        if f_: rows.append(dict(hid=k, rate=f_["rate"], se=f_["boot_se"], se_formula=f_["se"], n_pass=f_["passes"], organic_frac=round(frac[k], 3)))
    return pd.DataFrame(rows), {("Organic Soils" if c else "other"): [round(float(v) * 1000, 2) for v in ab] for c, ab in cyc.items()}


def hexagons(cfg):
    W, out = f"{cfg['paths']['work']}/grid30", cfg["paths"]["out"]; h = cfg["hexagons"]; m = cfg["masks"]
    df = pd.read_parquet(f"{cfg['paths']['work']}/segments_classified.parquet")   # step 3 (screened, labelled)
    F = df[df.land == "farmland"]; MIN = F[(F.order != "Organic") & (F.order != "not mapped") & ~F.near_peat]
    hx = hex_grid(cfg); nh = len(hx)
    H1, cyc1 = hex_rates(cfg, F, MIN, hx); hx = hx.merge(H1, on="hid", how="left")
    slope, tf = read(f"{W}/slope_deg.tif"); farmf = read(f"{W}/farmfrac.tif")[0]; near = read(f"{W}/near_organic.tif")[0].astype(bool)
    flat = slope < m["slope_max_deg"]; FARM = (farmf >= m["farmland_fraction_min"]) & flat
    hid = rasterize([(g, int(k)) for g, k in zip(hx.geometry, hx.hid)], out_shape=slope.shape, transform=tf, fill=-1, dtype="int32")
    inn = hid >= 0; nall = np.bincount(hid[inn], minlength=nh); nf = np.bincount(hid[inn & FARM], minlength=nh); hx["farm_share"] = (nf / np.maximum(nall, 1))[hx.hid.values]
    a = cfg["study_area"]; keep = (hx.farm_share >= h["min_farm_share"]) & hx.x.between(a["xmin"], a["xmax"]) & hx.y.between(a["ymin"], a["ymax"])   # centre in the study rectangle
    hx.loc[~keep, ["rate", "se"]] = np.nan; cap = cfg["maps"]["ref_cap_mm_yr"]; R = {}
    from lib.maps import ref_adapt
    for nm, f in (("LiDAR - Copernicus", "lidcop_rate.tif"), ("InSAR", "insar_up.tif")):
        a_ = read(f"{W}/{f}")[0]; v = np.where(np.isfinite(a_) & (np.abs(a_) < cap), a_, np.nan)
        ref = ref_adapt(v, FARM & ~near, cfg); P3 = FARM & np.isfinite(v) & np.isfinite(ref)   # farmland cells, farmland reference widened where sparse (Section 3.2)
        g = pd.DataFrame(dict(hid=hid[P3], v=(v - ref)[P3])).groupby("hid").v.agg(["median", "size"]); g = g[(g["size"] >= h["min_cells"]) & g.index.isin(hx.hid[keep])]
        col = "dem_rate" if nm.startswith("LiDAR") else "insar_rate"; hx[col] = hx.hid.map(g["median"])
        both = hx[col].notna() & hx.rate.notna()
        R[nm] = dict(hexagons=int(hx[col].notna().sum()), with_ICESat2=int(both.sum()), r=round(float(np.corrcoef(hx[col][both], hx.rate[both])[0, 1]), 2))
    # the same comparison against ICESat-2 hexagons over the InSAR period (to 30 June 2022), with OLS slopes (map on ICESat-2)
    te = pd.Timestamp(h["insar_period_end"]); t1 = te.year + te.dayofyear / 365.25
    H2, cyc2 = hex_rates(cfg, F, MIN, hx[["hid", "x", "y", "geometry"]], t1=t1)
    hx = hx.merge(H2[["hid", "rate", "se", "n_pass"]].rename(columns={"rate": "rate_2018_2022", "se": "se_2018_2022", "n_pass": "n_pass_2018_2022"}), on="hid", how="left")
    hx.loc[~keep.values, ["rate_2018_2022", "se_2018_2022"]] = np.nan
    for nm, col in (("LiDAR - Copernicus", "dem_rate"), ("InSAR", "insar_rate")):
        for per, ic in (("2018-2026", "rate"), ("2018-2022", "rate_2018_2022")):
            b_ = hx[col].notna() & hx[ic].notna(); x, y = hx[ic][b_].values, hx[col][b_].values
            R[nm][f"vs ICESat-2 {per}"] = dict(hexagons=int(b_.sum()), r=round(float(np.corrcoef(x, y)[0, 1]), 2), slope=round(float(np.polyfit(x, y, 1)[0]), 2))
    ok = hx.rate.notna()
    R["ICESat-2"] = dict(hexagons=int(ok.sum()), median_passes=float(hx.loc[ok, "n_pass"].median()), median_se=round(float(hx.loc[ok, "se"].median()), 1),
                         regional_annual_cycle_mm={"2018-2026": cyc1, "2018-2022": cyc2})
    hx = hx[hx.rate.notna() | hx.dem_rate.notna() | hx.insar_rate.notna()]
    hx.to_file(f"{out}/maps/hex_3km.gpkg"); json.dump(R, open(f"{out}/tables/hexagons.json", "w"), indent=1); print(json.dumps(R, indent=1), flush=True)


def run(cfg):
    os.makedirs(f"{cfg['paths']['out']}/maps", exist_ok=True); hexagons(cfg)


def max_rates(cfg):
    """Table 5: ICESat-2 rate of each InSAR zone (one target, regional mineral-soil reference on the same passes) per period;
    InSAR map median, LiDAR - Copernicus median, share on Organic Soils and distance inside the S-map boundary of the zone's 30 m farmland cells.
    v107: distance bands removed (user decision, 22 Nov 2026); the zones carry the interior result, with their distance inside Organic Soils."""
    from rasterio.features import shapes, geometry_mask
    from shapely.geometry import shape
    from lib.maps import ref_adapt
    work, out, raw = cfg["paths"]["work"], cfg["paths"]["out"], cfg["paths"]["raw"]; D = f"{work}/zones"; os.makedirs(f"{D}/polys", exist_ok=True); a = cfg["study_area"]
    org = gpd.read_file(f"{raw}/soils/smap_organic_waikato.shp").to_crs(2193).unary_union.buffer(0); T = {}
    with rasterio.open(f"{out}/maps/insar_rel_240m.tif") as d: I = d.read(1); tf = d.transform
    for t in cfg["insar"]["zones_mm_yr"]:
        m = np.isfinite(I) & (I <= t); T[f"insar_le_{-t}"] = gpd.GeoSeries([shape(s) for s, v in shapes(m.astype("uint8"), mask=m, transform=tf) if v == 1], crs=2193).unary_union
    study = gpd.GeoDataFrame(geometry=[__import__("shapely").geometry.box(a["xmin"], a["ymin"], a["xmax"], a["ymax"])], crs=2193); study.to_file(f"{D}/polys/study.gpkg")
    P = {"Oct 2018 - Jun 2022": ("2018-10-13", "2022-06-30"), "Jul 2022 - Oct 2026": ("2022-07-01", "2026-10-20"), "2018-2026": ("2018-10-13", "2026-10-20")}
    R = {}
    # v105: zone and band rates with the pipeline's own ICESat-2 code (step 3: pure 3 x 3 farmland, mineral-soil reference, pass-and-DEM residual screen),
    # as Table 4 and Table 6 (before v105 is2rate, which classes land cover at the segment point only). Target: farmland segments inside the polygon;
    # reference: regional mineral-soil farmland on the same passes, without the target's own segments.
    import shapely
    from steps.s03_icesat2_rates import classify
    from lib.is2 import diff_series, fit, step_fit
    rc = cfg["icesat2"]["rates"]; ms = rc["min_segments_per_pass"]
    df, _ = classify(cfg); F = df[df.land == "farmland"]; MIN = F[(F.order != "Organic") & (F.order != "not mapped") & ~F.near_peat]
    (a0_, a1_), (u0_, u1_), (v0_, v1_) = rc["periods"]   # decimal-year limits of step 3
    P2 = {"Oct 2018 - Jun 2022": (u0_, u1_), "Jul 2022 - Oct 2026": (v0_, v1_), "2018-2026": (a0_, a1_), "Jul 2022 - Jul 2026 (sensitivity)": (v0_, 2026.0 + 212 / 365.25)}
    g30, tf30g = read(f"{cfg['paths']['work']}/grid30/slope_deg.tif")
    def cell_in(mask, d):
        r_ = ((tf30g.f - d.Y.values) // 30).astype(int); c_ = ((d.X.values - tf30g.c) // 30).astype(int); ok = (r_ >= 0) & (r_ < mask.shape[0]) & (c_ >= 0) & (c_ < mask.shape[1])
        o = np.zeros(len(d), bool); o[ok] = mask[r_[ok], c_[ok]]; return o
    FIG6 = []   # v120: per-pass series for Figure 6 (Organic-Soil farmland and the -40 zone)
    for nm, g in T.items():
        gpd.GeoDataFrame(geometry=[g], crs=2193).to_file(f"{D}/polys/{nm}.gpkg"); R[nm] = {"km2": round(g.area / 1e6, 1), "km2_in_study_rectangle": round(g.intersection(study.geometry.iloc[0]).area / 1e6, 1)}
        gm = rasterize([(g, 1)], out_shape=g30.shape, transform=tf30g, fill=0, dtype="uint8") > 0   # polygon on the 30 m grid (fast membership)
        inside = cell_in(gm, F); tgt = F[inside]; ref = MIN[~cell_in(gm, MIN)]
        R[nm]["target segments"] = int(len(tgt))
        for pk, (t0, t1) in P2.items():
            ser = diff_series(tgt, ref, ms); ser = ser[(ser.t >= t0) & (ser.t < t1)]
            f_ = fit(ser, n_boot=0); R[nm][pk] = dict(rate=f_["rate"], se=f_["se"], passes=f_["passes"]) if f_ else None
        ser_all = diff_series(tgt, ref, ms); st = step_fit(ser_all, rc["break_date"], n_boot=rc["n_boot"]); R[nm]["two trends and a step"] = st   # v120: with bootstrap CI
        if nm == "insar_le_40": FIG6.append(ser_all.assign(series="InSAR zone -40"))   # v114 (review A1)
        print(nm, R[nm], flush=True)
    # v114 (review A3): Organic-Soil farmland inside and outside the -10 zone, and mineral-soil farmland inside it (2018-2022 and 2018-2026)
    g10 = rasterize([(T["insar_le_10"], 1)], out_shape=g30.shape, transform=tf30g, fill=0, dtype="uint8") > 0
    ORGF = F[F.order == "Organic"]; MINr = MIN[~cell_in(g10, MIN)]
    FIG6.append(diff_series(ORGF, MIN, ms).assign(series="Organic-Soil farmland"))   # as step 3
    os.makedirs(f"{out}/plots", exist_ok=True); pd.concat(FIG6)[["series", "pass", "col", "t", "r"]].to_csv(f"{out}/plots/fig6_series.csv", index=False)
    sp = {"Organic-Soil farmland inside the -10 zone": ORGF[cell_in(g10, ORGF)], "Organic-Soil farmland outside the -10 zone": ORGF[~cell_in(g10, ORGF)],
          "mineral-soil farmland inside the -10 zone": MIN[cell_in(g10, MIN)]}
    R["Organic-Soil farmland and the -10 zone"] = {k: {"segments": int(len(v)), **{pk: fit(diff_series(v, MINr, ms)[lambda s: (s.t >= P2[pk][0]) & (s.t < P2[pk][1])], n_boot=0)
                                                    for pk in ("Oct 2018 - Jun 2022", "2018-2026")}} for k, v in sp.items()}
    print(R["Organic-Soil farmland and the -10 zone"], flush=True)
    # v120: the is2rate check of the -40 zone removed (hexagons now use the pipeline code; one ICESat-2 code path)
    W = f"{work}/grid30"; m = cfg["masks"]; slope, tf30 = read(f"{W}/slope_deg.tif"); grass = read(f"{W}/grassfrac.tif")[0]
    near = read(f"{W}/near_organic.tif")[0].astype(bool); orgm = read(f"{W}/organic.tif")[0].astype(bool)
    PAST = (read(f"{W}/farmfrac.tif")[0] >= m["farmland_fraction_min"]) & (slope < m["slope_max_deg"]); C = read(f"{W}/lidcop_rate.tif")[0]; cap = cfg["maps"]["ref_cap_mm_yr"]   # farmland
    ref = ref_adapt(C, PAST & ~near, cfg); ok = PAST & orgm & np.isfinite(C) & (np.abs(C) < cap) & np.isfinite(ref); relc = C - ref
    from scipy.ndimage import distance_transform_edt
    dist_in = distance_transform_edt(orgm) * 30.0   # distance of an Organic Soil cell to the S-map boundary (m)
    okc = PAST & np.isfinite(C) & (np.abs(C) < cap) & np.isfinite(ref)
    R["Organic-Soil farmland, median distance inside Organic Soils (m)"] = round(float(np.median(dist_in[PAST & orgm])), -1)   # baseline for the zones
    for nm, g in T.items():
        mk = rasterize([(g, 1)], out_shape=C.shape, transform=tf30, fill=0, dtype="uint8") > 0; fz = mk & PAST   # zone farmland cells (cell centre inside)
        R[nm]["InSAR map median"] = round(float(np.nanmedian(I[np.isfinite(I) & (I <= -int(nm.split('_')[-1]))])), 1)
        R[nm]["LiDAR - Copernicus median"] = round(float(np.median(relc[mk & okc])), 1)
        R[nm]["farmland share on Organic Soils"] = round(float(orgm[fz].mean()), 3)
        R[nm]["median distance inside Organic Soils (m), Organic farmland cells"] = round(float(np.median(dist_in[fz & orgm])), -1)
    json.dump(R, open(f"{out}/tables/max_rates.json", "w"), indent=1)


_run_hex = run
def run(cfg):
    _run_hex(cfg); max_rates(cfg)


_run_max = run
def run(cfg):
    _run_max(cfg)   # v115: repeat_vs_insar (Hauraki Plains quarters) removed with the repeat LiDAR map


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
    PAST = (read(f"{W}/farmfrac.tif")[0] >= m["farmland_fraction_min"]) & (slope < m["slope_max_deg"]); C = read(f"{W}/lidcop_rate.tif")[0]; PEAT = PAST & (org > 0)   # farmland
    for cap in cfg["sensitivity"]["ref_caps_mm_yr"]:
        c2 = dict(cfg, maps=dict(cfg["maps"], ref_cap_mm_yr=cap if cap else 1e9)); ref = ref_adapt(C, PAST & ~near, c2)
        okc = PEAT & np.isfinite(C) & (np.abs(C) < (cap if cap else 1e9)) & np.isfinite(ref); R["reference_cap_mm_yr"][str(cap)] = round(float(np.median((C - ref)[okc])), 1)
    json.dump(R, open(f"{out}/tables/sensitivity.json", "w"), indent=1); print("sensitivity", R, flush=True); return R


_run_cmp = run
def run(cfg):
    _run_cmp(cfg); sensitivity(cfg)
