"""Step 2 (Table S1): ICESat-2 segments and residuals.
(a) SlideRule atl06p request over the study area (2 x 2 tiles of its bounding box + 2 km), cached in data/raw/icesat2/;
(b) quality screening; heights to NZTM 2000 / NZVD2016 (NZGeoid2016, PROJ); one row per segment and time;
(c) labels: study polygon, Organic Soils and the 500 m buffer zone, WorldCover class and 3 x 3 pure window;
(d) LiDAR height (bilinear, 4 m) and terrain slope from the 2024 and 2021 Waikato DEMs.
Output: data/work/segments.parquet (one row per segment; residuals are formed in step 3, where the DEM and the outlier screen are applied)."""
import glob, json, os, numpy as np, pandas as pd, geopandas as gpd, pyproj
from lib.sample import lidar_at, worldcover


def _tile(args):
    """One tile of one area, requested in windows of months_per_request. Each window is retried (up to retries_per_window) while the server reports a
    failure on any granule (alert code < 0, 'Failure on resource'); windows are cached in icesat2/parts/ so a restart resumes; failures left are logged."""
    poly, a_, b_, out, p, area, ij = args
    import json, os, pandas as pd, geopandas as gpd
    from sliderule import sliderule, icesat2
    sliderule.init("slideruleearth.io", verbose=False); parts = []; left = {}
    pdir = os.path.join(os.path.dirname(out), "parts"); os.makedirs(pdir, exist_ok=True)
    ed = [a_] + [d for d in pd.date_range(a_, b_, freq=f"{p['months_per_request']}MS").strftime("%Y-%m-%d") if a_ < d < b_] + [b_]
    for t0, t1 in zip(ed[:-1], ed[1:]):
        pf = f"{pdir}/{area}_{ij}_{t0}.parquet"
        if os.path.exists(pf): parts.append(gpd.read_parquet(pf)); continue
        parms = {"poly": poly, "t0": t0 + "T00:00:00Z", "t1": t1 + "T00:00:00Z", "srt": icesat2.SRT_LAND, "cnf": icesat2.CNF_SURFACE_LOW, "atl08_class": ["atl08_ground"],
                 "len": float(p["len"]), "res": float(p["res"]), "ats": float(p["ats"]), "cnt": p["cnt"], "maxi": p["maxi"], "H_min_win": float(p["H_min_win"]),
                 "sigma_r_max": p["sigma_r_max"], "timeout": p["timeout_s"]}
        best = None; prev = None
        for k in range(p["retries_per_window"]):
            fails = set()
            def alert(rec, session, fails=fails):
                if rec["code"] < 0 and "Failure on resource" in rec["text"]: fails.add(rec["text"].split("Failure on resource ")[1].split(":")[0])
            g = icesat2.atl06p(parms, callbacks={"exceptrec": alert})
            if best is None or len(fails) < len(best[1]): best = (g, set(fails))
            print(area, ij, t0, "try", k + 1, len(g), "failed", len(fails), sorted(fails), flush=True)
            if not fails or fails == prev: break   # complete, or the same granules failed twice: unreadable on the server
            prev = set(fails)
        g, f = best
        if f: left[t0] = sorted(f)
        if not f or f == prev: g.to_parquet(pf)   # cache complete windows, and windows whose failures repeat (persistent)
        if len(g): parts.append(g)
    pd.concat(parts).sort_index().to_parquet(out)
    json.dump(left, open(out.replace(".parquet", "_failures.json"), "w"), indent=1); return out, sum(len(v) for v in left.values())


def _area(jobs):
    return [_tile(j) for j in jobs]


def fetch(cfg, raw, areas=(("Waikato", "Nth_Waikato.shp"),)):
    """SlideRule atl06p over 2 x 2 tiles of each area's bounding box + 2 km; tiles already downloaded are skipped; areas in parallel.
    v99: requests in windows of months_per_request (6). One request for the whole record returned at most ~300 granules (CMR limit),
    so the v98 downloads (17 Nov 2026) stopped in June-July 2026 and missed later passes. Areas run in parallel, tiles and windows within an area in sequence."""
    import concurrent.futures as cf
    from shapely.geometry import box
    p = cfg["icesat2"]["sliderule"]; jobs = []
    for area, shp in areas:
        x0, y0, x1, y1 = gpd.read_file(f"{raw}/boundaries/{shp}").to_crs(2193).total_bounds + [-2000, -2000, 2000, 2000]
        xs, ys = [x0, (x0 + x1) / 2, x1], [y0, (y0 + y1) / 2, y1]
        for i in range(2):
            for j in range(2):
                out = f"{raw}/icesat2/{area}_{i}{j}_atl06sr.parquet"
                if os.path.exists(out): continue
                reg = gpd.GeoSeries([box(xs[i], ys[j], xs[i + 1], ys[j + 1])], crs=2193).to_crs(4326)
                jobs.append(([{"lon": x, "lat": y} for x, y in reg.iloc[0].exterior.coords], p["t0"], p["t1"], out, p, area, f"{i}{j}"))
    if jobs:   # one worker per area, its tiles and windows in sequence (12 parallel tile requests overloaded the server: 10-15 min per request and read failures)
        by = {}
        for j in jobs: by.setdefault(j[5], []).append(j)
        with cf.ProcessPoolExecutor(len(by)) as ex:
            for res in ex.map(_area, by.values()):
                for o, nf in res: print("written", o, "granule-tracks still failing", nf, flush=True)
    # v99: granule-tracks that failed on the server in the final try (persistent: same set on two consecutive tries) are reported, not silently lost
    bad = {f: n for f in glob.glob(f"{raw}/icesat2/*_failures.json") if (n := sum(len(v) for v in json.load(open(f)).values()))}
    print("SlideRule granule-tracks failing after retries (persistent):", bad or "none", flush=True)


def run(cfg):
    raw, work = cfg["paths"]["raw"], cfg["paths"]["work"]; os.makedirs(work, exist_ok=True); q = cfg["icesat2"]["screening"]
    fetch(cfg, raw, [("Waikato", "Nth_Waikato.shp")] + ([(k, f"{k}.shp") for k in cfg["south_island"]] if cfg["icesat2"]["sliderule"].get("fetch_south_island_in_step2", True) else []))   # v99: all three areas, areas in parallel
    pyproj.datadir.append_data_dir(os.path.abspath(f"{raw}/proj"))   # NZGeoid2016 grid (nz_linz_nzgeoid2016.tif, LINZ via PROJ-data)
    TR = pyproj.Transformer.from_crs("EPSG:4979", "EPSG:2193+7839", always_xy=True)
    fr = []
    for f in sorted(glob.glob(f"{raw}/icesat2/Waikato_??_atl06sr.parquet")):
        g = gpd.read_parquet(f)
        g = g[(g.n_fit_photons >= q["n_photons_min"]) & (g.h_sigma < q["h_sigma_max"]) & (g.dh_fit_dx.abs() < q["dh_fit_dx_max"]) & (g.w_surface_window_final < q["w_surface_window_max"])]
        lon, lat = g.geometry.x.values, g.geometry.y.values
        X, Y, Hh = TR.transform(lon, lat, g.h_mean.values)
        fr.append(pd.DataFrame(dict(lon=lon, lat=lat, X=X, Y=Y, H=Hh, t=(g.index.year + (g.index.dayofyear - .5) / 365.25).values,
                                    **{"pass": [f"{x:%Y%m%d}_{r}" for x, r in zip(g.index, g.rgt)]})))
    df = pd.concat(fr, ignore_index=True).drop_duplicates(subset=["X", "Y", "t"]); df = df[np.isfinite(df.H)].reset_index(drop=True)
    print("segments after screening:", len(df), flush=True)
    pts = gpd.GeoSeries(gpd.points_from_xy(df.X, df.Y), crs=2193); P = gpd.GeoDataFrame(geometry=pts.values, crs=2193); P["i"] = np.arange(len(P))
    df["in_area"] = pts.within(gpd.read_file(f"{raw}/boundaries/Nth_Waikato.shp").to_crs(2193).unary_union).values
    org = gpd.read_file(f"{raw}/soils/smap_organic_waikato.shp").to_crs(2193)
    df["peat"] = False; df.loc[gpd.sjoin(P, org[["geometry"]], predicate="within")["i"].unique(), "peat"] = True
    df["near_peat"] = False; df.loc[gpd.sjoin(P, gpd.GeoDataFrame(geometry=org.geometry.buffer(cfg["masks"]["buffer_m"]).values, crs=2193), predicate="within")["i"].unique(), "near_peat"] = True
    df["wc"], df["wc_solid"] = worldcover(df.lon.values, df.lat.values)
    TI = json.load(open(f"{raw}/tile_index.json"))
    for c in cfg["lidar"]["waikato"]:
        print("sampling", c, flush=True); L, S = lidar_at(df.X.values, df.Y.values, TI[c]["tiles"], TI[c]["base"], factor=cfg["lidar"]["resolutions_m"]["icesat2"])
        df["L" + c[-4:]] = L; df["S" + c[-4:]] = S
    df.to_parquet(f"{work}/segments.parquet"); print("written", f"{work}/segments.parquet", len(df), flush=True)
