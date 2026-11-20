"""Step 3 (Table S1): ICESat-2 rates of land classes, and the survey bias of the 2021 LiDAR.
Input: data/work/segments.parquet (step 2). Outputs: data/work/segments_classified.parquet, outputs/tables/icesat2_rates.json and
icesat2_rates.xlsx, outputs/plots/fig6_series.xlsx, data/work/bias2021_10km.csv."""
import json, os, numpy as np, pandas as pd, geopandas as gpd
from lib.is2 import per_pass, fit, diff_series, design, segmented, amplitude, separate_difference

FARM, BUILT, WET = (30, 40), 50, 90


def classify(cfg):
    """Residual against the 2024 LiDAR where it has a height, else the 2021 LiDAR; then the slope test on that DEM (as is2rate and Table S1); screening; land and soil labels.
    v99: before, a segment steep (>= 2 deg) in the 2024 LiDAR fell back to the 2021 LiDAR (drains/earthworks dug after 2021 compared with the earlier surface)."""
    raw, work = cfg["paths"]["raw"], cfg["paths"]["work"]; a = cfg["study_area"]; q = cfg["icesat2"]["screening"]
    df = pd.read_parquet(f"{work}/segments.parquet")
    u24 = df.L2024.notna()
    df["L"] = np.where(u24, df.L2024, df.L2021); df["S"] = np.where(u24, df.S2024, df.S2021); df["col"] = np.where(u24, "2024", "2021")
    df = df[df.L.notna() & (df.S < q["terrain_slope_max_deg"]) & df.X.between(a["xmin"], a["xmax"], inclusive="left") & df.Y.between(a["ymin"], a["ymax"], inclusive="left")].copy()
    df["r"] = df.H - df.L
    df = df[(df.r - df.groupby(["pass", "col"]).r.transform("median")).abs() <= q["residual_max_m"]].copy()   # v99+: within 0.5 m of the median of its pass and DEM (all land), as is2rate; was the record-wide farmland median
    farm = df.wc_solid & df.wc.isin(FARM); built = df.wc_solid & (df.wc == BUILT); wet = df.wc == WET
    df = df[farm | built | wet].copy(); df["land"] = np.select([df.wc == WET, df.wc == BUILT], ["natural wetland", "built-up"], "farmland")
    sm = gpd.read_file(f"{raw}/soils/smap_soilorder_waikato.gpkg")[["NZSCOrder", "geometry"]]
    j = gpd.sjoin(gpd.GeoDataFrame(df[["X"]], geometry=gpd.points_from_xy(df.X, df.Y), crs=2193), sm, how="left", predicate="within"); j = j[~j.index.duplicated()]
    df["order"] = j.NZSCOrder.fillna("not mapped").values
    b = gpd.read_file(f"{raw}/soils/Waikato_peat_bodies_stats.gpkg").to_crs(2193).reset_index(drop=True)
    j = gpd.sjoin(gpd.GeoDataFrame(df[["X"]], geometry=gpd.points_from_xy(df.X, df.Y), crs=2193), b[["name", "geometry"]], how="left", predicate="within"); j = j[~j.index.duplicated()]
    df["body"] = j["name"].values
    df.to_parquet(f"{work}/segments_classified.parquet"); return df, b


def survey_bias(cfg, raw_df):
    """Survey offset of each LiDAR survey on farmland: ICESat-2 minus LiDAR on mineral-soil farmland minus the same on built-up land, per pass
    (passes within one year of the survey), median over passes. Not applied: offsets cancel in the local reference (Section 3.3)."""
    out = {}
    farm = raw_df.wc_solid & raw_df.wc.isin(FARM); built = raw_df.wc_solid & (raw_df.wc == BUILT)
    for yr, tm in (("2021", 2021.14), ("2024", 2024.24)):
        e = raw_df[raw_df["L" + yr].notna() & (raw_df["S" + yr] < 2) & ((raw_df.t - tm).abs() < 1.0)].copy(); e["r"] = e.H - e["L" + yr]; e = e[(e.r - e.r.median()).abs() < 0.5]
        f_, b_ = e[farm.loc[e.index] & ~e.near_peat], e[built.loc[e.index]]
        jj = pd.concat([f_.groupby("pass").r.median(), b_.groupby("pass").r.median()], axis=1, keys=["f", "b"]).dropna()
        out[f"survey_{yr}"] = dict(farm_median_cm=round(float(f_.r.median()) * 100, 1), farm_minus_built_same_pass_cm=round(float(np.median(jj.f - jj.b)) * 100, 1), paired_passes=len(jj))
    return out


def run(cfg):
    work, out = cfg["paths"]["work"], cfg["paths"]["out"]; os.makedirs(f"{out}/tables", exist_ok=True); os.makedirs(f"{out}/plots", exist_ok=True)
    rc = cfg["icesat2"]["rates"]; ms, nb = rc["min_segments_per_pass"], rc["n_boot"]; a = cfg["study_area"]
    df, bodies = classify(cfg)
    F = df[df.land == "farmland"]; ORG = F[F.order == "Organic"]; MIN = F[(F.order != "Organic") & (F.order != "not mapped") & ~F.near_peat]
    BU = df[df.land == "built-up"]
    hx0, hy0, hx1, hy1 = a["hauraki_box"]; HAU = MIN[MIN.X.between(hx0, hx1) & MIN.Y.between(hy0, hy1)]
    kp = bodies.geometry.iloc[int(bodies.index[bodies.name == "Kopuatai"][0])]
    BOG = df[(df.land == "natural wetland") & (df.order == "Organic") & (df.body == "Kopuatai")]
    LOC = MIN[gpd.GeoSeries(gpd.points_from_xy(MIN.X, MIN.Y), crs=2193).within(kp.buffer(rc["bog_reference_km"] * 1000)).values]
    S = {"Organic-Soil farmland": diff_series(ORG, MIN, ms), "Kopuatai bog": diff_series(BOG, LOC, ms)}
    A = {"mineral-soil farmland": per_pass(MIN, ms), "Organic-Soil farmland": per_pass(ORG, ms), "built-up land": per_pass(BU, ms),
         "Hauraki Plains mineral-soil farmland": per_pass(HAU, ms), "mineral-soil farmland within 5 km of Kopuatai": per_pass(LOC, ms)}
    R = {"last pass in study area": str(sorted(df["pass"].unique())[-1]), "segments": {k: int(len(v)) for k, v in (("Organic-Soil farmland", ORG), ("mineral-soil farmland", MIN), ("built-up land", BU), ("Kopuatai bog", BOG))}}
    R["relative to mineral-soil farmland"] = {k: fit(v, n_boot=nb) for k, v in S.items()}
    R["Organic-Soil farmland by land cover"] = {k: fit(diff_series(ORG[ORG.wc == c], MIN, ms), n_boot=nb) for k, c in (("grassland", 30), ("cropland", 40))}
    R["absolute"] = {k: fit(v, n_boot=nb) for k, v in A.items()}
    R["relative to built-up land"] = {k: fit(diff_series(g, BU, ms), n_boot=nb) for k, g in (("mineral-soil farmland", MIN), ("Hauraki Plains mineral-soil farmland", HAU))}
    (t0, t1), (u0, u1), (v0, v1) = rc["periods"]
    R["periods"] = {k: {"Oct 2018 - Jun 2022": fit(v[(v.t >= u0) & (v.t < u1)], n_boot=nb), "Jul 2022 - Oct 2026": fit(v[(v.t >= v0) & (v.t < v1)], n_boot=nb)} for k, v in S.items()}
    nb6 = rc["n_boot_table6"]   # v98: Table 6 CIs from 20,000 resamples (bog CIs seed-sensitive at 2,000)
    R["segmented"] = {k: segmented(v, rc["break_date"], n_boot=nb6) for k, v in S.items()}
    R["seasonal amplitude"] = {k: amplitude(A[k], rc["n_boot_amplitude"]) for k in ("Organic-Soil farmland", "mineral-soil farmland")}
    R["separate fits, before minus after"] = {k: separate_difference(v, rc["break_date"], n_boot=nb6) for k, v in S.items()}
    # v99: leave-one-pass-out range of the separate-fit change (it rests near zero and on single passes; checker v98: one pass with 12 target segments moved it 2 mm/yr)
    for k, v in S.items():
        jk = [separate_difference(v[v["pass"] != q], rc["break_date"], n_boot=0)["before_minus_after"] for q in v["pass"].unique()]
        R["separate fits, before minus after"][k]["leave_one_pass_out"] = [round(-max(jk), 1), round(-min(jk), 1)]   # as after - before
    LOC3 = MIN[gpd.GeoSeries(gpd.points_from_xy(MIN.X, MIN.Y), crs=2193).within(kp.buffer(3000)).values]
    R["Kopuatai bog, reference within 3 km"] = fit(diff_series(BOG, LOC3, ms), n_boot=0)
    R["Organic-Soil farmland without annual cycle"] = fit(S["Organic-Soil farmland"], n_boot=0, seasonal=False)
    rat = {o: fit(diff_series(F[F.order == o], MIN, ms), n_boot=nb) for o in F.order.unique() if o != "not mapped" and (F.order == o).sum() >= 300}
    rat["Kopuatai bog"] = R["relative to mineral-soil farmland"]["Kopuatai bog"]
    rr = {k: round(v["se"] / v["boot_se"], 2) for k, v in rat.items() if v}
    R["SE ratio formula/bootstrap (soil orders and bog, 2018-2026)"] = dict(ratios=rr, median=round(float(np.median(list(rr.values()))), 2), min=min(rr.values()), max=max(rr.values()))
    seg = pd.read_parquet(f"{work}/segments.parquet", columns=["X", "Y"])
    R["segments after screening in the study area"] = int((seg.X.between(a["xmin"], a["xmax"], inclusive="left") & seg.Y.between(a["ymin"], a["ymax"], inclusive="left")).sum())
    R["survey bias"] = survey_bias(cfg, pd.read_parquet(f"{work}/segments.parquet"))
    json.dump(R, open(f"{out}/tables/icesat2_rates.json", "w"), indent=1)
    with pd.ExcelWriter(f"{out}/plots/fig6_series.xlsx") as xw:
        for k, v in {**S, **{f"abs {k}": A[k] for k in ("Organic-Soil farmland", "mineral-soil farmland")}}.items(): v.to_excel(xw, sheet_name=k[:31], index=False)
    print(json.dumps(R, indent=1), flush=True)
