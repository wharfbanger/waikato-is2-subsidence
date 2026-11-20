"""Step 8 (Table S1): Canterbury and Southland (study areas of Harvey et al. 2026).
SlideRule download (cached in data/raw/icesat2/<Area>_??_atl06sr.parquet, as step 2), screening and NZVD2016 heights as step 2; LiDAR height (bilinear, 4 m)
from each regional DEM in order of precedence (one constant per DEM in the fit); farmland = WorldCover grassland/cropland (pure 3 x 3) inside the
study polygon and outside the reference town; built-up = WorldCover built-up (pure 3 x 3); rates as step 3; published InSAR sampled at the segments.
Output: outputs/tables/south_island.json."""
import glob, json, os, numpy as np, pandas as pd, geopandas as gpd, pyproj, rasterio
from lib.sample import lidar_at, worldcover
from lib.is2 import per_pass, fit, diff_series


def segments(cfg, area, A):
    raw, work = cfg["paths"]["raw"], cfg["paths"]["work"]; S = f"{work}/segments_{area}.parquet"; q = cfg["icesat2"]["screening"]
    if os.path.exists(S): return pd.read_parquet(S)
    from steps.s02_icesat2 import fetch; fetch(cfg, raw, [(area, f"{area}.shp")])   # v99: download in 2-month windows if not cached
    pyproj.datadir.append_data_dir(os.path.abspath(f"{raw}/proj")); TR = pyproj.Transformer.from_crs("EPSG:4979", "EPSG:2193+7839", always_xy=True); fr = []
    for f in sorted(glob.glob(f"{raw}/icesat2/{area}_??_atl06sr.parquet")):
        g = gpd.read_parquet(f); g = g[(g.n_fit_photons >= q["n_photons_min"]) & (g.h_sigma < q["h_sigma_max"]) & (g.dh_fit_dx.abs() < q["dh_fit_dx_max"]) & (g.w_surface_window_final < q["w_surface_window_max"])]
        lon, lat = g.geometry.x.values, g.geometry.y.values; X, Y, Hh = TR.transform(lon, lat, g.h_mean.values)
        fr.append(pd.DataFrame(dict(lon=lon, lat=lat, X=X, Y=Y, H=Hh, t=(g.index.year + (g.index.dayofyear - .5) / 365.25).values, **{"pass": [f"{x:%Y%m%d}_{r}" for x, r in zip(g.index, g.rgt)]})))
    df = pd.concat(fr, ignore_index=True).drop_duplicates(subset=["X", "Y", "t"]); df = df[np.isfinite(df.H)].reset_index(drop=True)
    df["in_area"] = gpd.GeoSeries(gpd.points_from_xy(df.X, df.Y), crs=2193).within(gpd.read_file(f"{raw}/boundaries/{area}.shp").to_crs(2193).unary_union).values
    df["wc"], df["wc_solid"] = worldcover(df.lon.values, df.lat.values, tiles=A["worldcover_tiles"])
    TI = json.load(open(f"{raw}/tile_index.json")); df["L"] = np.nan; df["S"] = np.nan; df["col"] = ""
    for c in A["dems"]:
        m = df.L.isna().values; L, Sl = lidar_at(df.X.values[m], df.Y.values[m], TI[c]["tiles"], TI[c]["base"], factor=cfg["lidar"]["resolutions_m"]["icesat2"])
        idx = np.flatnonzero(m); ok = np.isfinite(L); df.loc[idx[ok], "L"] = L[ok]; df.loc[idx[ok], "S"] = Sl[ok]; df.loc[idx[ok], "col"] = c
    with rasterio.open(f"{raw}/south_island/{A['insar']}") as d:
        v = np.array([x[0] for x in d.sample(zip(df.X, df.Y))], float); v[(v < -1000) | (v > 1000)] = np.nan
        if d.nodata is not None: v[v == d.nodata] = np.nan
    df["insar"] = v; df.to_parquet(S); return df


def run(cfg):
    out = cfg["paths"]["out"]; rc = cfg["icesat2"]["rates"]; ms, mb = rc["min_segments_per_pass"], rc["min_segments_per_pass"]; R = {}
    for area, A in cfg["south_island"].items():
        df = segments(cfg, area, A); lo0, la0, lo1, la1 = A["town_box"]
        d = df[df.L.notna() & (df.S < cfg["icesat2"]["screening"]["terrain_slope_max_deg"])].copy(); d["r"] = d.H - d.L
        d = d[(d.r - d.groupby(["pass", "col"]).r.transform("median")).abs() <= cfg["icesat2"]["screening"]["residual_max_m"]]; town = d.lon.between(lo0, lo1) & d.lat.between(la0, la1)   # v99+: pass and DEM median, as is2rate
        ag = d[d.in_area & d.wc_solid & d.wc.isin(cfg["masks"]["worldcover_farmland"]) & ~town]; bu = d[d.in_area & d.wc_solid & (d.wc == cfg["masks"]["worldcover_builtup"])]   # v99+: built-up land in the study area, as the Waikato (was anywhere in the tiles, +2 km)
        q = ag.insar; fast = ag[q < -20]
        R[area] = dict(segments_used=len(d), farmland_absolute=fit(per_pass(ag, ms), n_boot=rc["n_boot"]), farmland_minus_built_up=fit(diff_series(ag, bu, mb), n_boot=rc["n_boot"]),
                       InSAR_on_farmland=round(float(np.nanmedian(ag.insar)), 1), InSAR_on_built_up=round(float(np.nanmedian(bu.insar)), 1),
                       InSAR_fastest_farmland=dict(InSAR=round(float(np.nanmedian(fast.insar)), 1), minus_built_up=fit(diff_series(fast, bu, mb), n_boot=rc["n_boot"])))
        print(area, json.dumps(R[area]), flush=True)
    json.dump(R, open(f"{out}/tables/south_island.json", "w"), indent=1)
