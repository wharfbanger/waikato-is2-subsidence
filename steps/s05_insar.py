"""Step 5 (Table S1): InSAR vertical velocity.
(a) LiCSBAS in WSL: steps 01-13 (batch_LiCSBAS.sh) unless a time series already exists in the job folder; no noise mask (step 15) by default,
    for ease of reproduction (insar.mask: true applies step 15 with the config thresholds in a copy of the time series); step 16 not used;
(b) vertical = line-of-sight velocity / up component of the unit vector; masked; averaged to the 30 m grid;
(c) each farmland cell referenced locally (as step 4); 240 m map; class values at the ICESat-2 segment locations; InSAR zones.
Outputs: data/work/insar/vel_up.tif, grid30/insar_rel.tif, outputs/maps/insar_rel_240m.tif, outputs/tables/insar.json."""
import json, os, subprocess, numpy as np, pandas as pd, rasterio
from rasterio.transform import from_origin
from rasterio.warp import reproject, Resampling
from lib.maps import read, rel, blocks240
from steps.s04_dem_rates import write


def wsl_path(p): p = os.path.abspath(p).replace("\\", "/"); return f"/mnt/{p[0].lower()}{p[2:]}"


def licsbas(cfg, D):
    """Writes and runs a WSL script; returns the folder holding results/ of the pipeline's step-15 run."""
    i = cfg["insar"]; s15 = i["step15"]; T = f"{i['wsl_pipeline_dir']}/TS"
    step15 = f"""mkdir -p {T}; [ -f {T}/cum.h5 ] || cp $J/TS_GEOCml1clip/cum.h5 {T}/
rm -rf {T}/results {T}/info {T}/network; cp -r $J/TS_GEOCml1clip/results $J/TS_GEOCml1clip/info $J/TS_GEOCml1clip/network {T}/
cd {i['wsl_pipeline_dir']}
LiCSBAS15_mask_ts.py -t TS -c {s15['coh_thre']} -u {s15['n_unw_r_thre']} -v {s15['vstd_thre']} -T {s15['maxTlen_thre']} -g {s15['n_gap_thre']} -s {s15['stc_thre']} -i {s15['n_ifg_noloop_thre']} -l {s15['n_loop_err_thre']} -r {s15['resid_rms_thre']} | sed -n '/Noise index/,/Kept/p'
cp {T}/results/mask {wsl_path(D)}/"""
    sh = f"""#!/bin/bash
set -e
source ~/miniconda3/etc/profile.d/conda.sh; conda activate {i['conda_env']}
export LiCSBAS_HOME={i['licsbas_home']}; export PATH=$LiCSBAS_HOME/bin:$PATH; export PYTHONPATH=$LiCSBAS_HOME/LiCSBAS_lib:$PYTHONPATH
J={i['wsl_job_dir']}
if [ ! -f $J/TS_GEOCml1clip/cum.h5 ]; then
  cd $J; cp $LiCSBAS_HOME/batch_LiCSBAS.sh .
  sed -i 's/^start_step=.*/start_step="01"/; s/^end_step=.*/end_step="13"/; s/^nlook=.*/nlook="1"/' batch_LiCSBAS.sh
  sed -i 's#^p01_frame=.*#p01_frame="{i['frame']}"#; s#^p01_start_date=.*#p01_start_date="{i['start']}"#; s#^p01_end_date=.*#p01_end_date="{i['end']}"#' batch_LiCSBAS.sh
  sed -i 's#^p05_clip_range_geo=.*#p05_clip_range_geo="{'/'.join(str(x) for x in i['clip'])}"#; s#^p12_multi_prime=.*#p12_multi_prime="y"#' batch_LiCSBAS.sh
  ./batch_LiCSBAS.sh
fi
{step15}
echo "interferograms_input $(grep -c '^2' $J/TS_GEOCml1clip/info/11ifg_stats.txt)"; echo "acquisitions_input $(grep '^2' $J/TS_GEOCml1clip/info/11ifg_stats.txt | awk '{{print $1}}' | tr _ '\n' | sort -u | wc -l)"
cp $J/TS_GEOCml1clip/results/vel {wsl_path(D)}/; cp $J/GEOCml1clip/U.geo $J/GEOCml1clip/EQA.dem_par {wsl_path(D)}/
"""
    p = f"{D}/run_licsbas.sh"; open(p, "w", newline="\n").write(sh)
    r = subprocess.run(["wsl.exe", "bash", wsl_path(p)], capture_output=True, text=True); print(r.stdout, r.stderr[-2000:], flush=True)
    if r.returncode: raise RuntimeError("LiCSBAS failed")
    return r.stdout


def run(cfg):
    work, out = cfg["paths"]["work"], cfg["paths"]["out"]; D = f"{work}/insar"; os.makedirs(D, exist_ok=True); W = f"{work}/grid30"; a = cfg["study_area"]
    log = licsbas(cfg, D)
    par = dict(l.split(":")[0:2] for l in open(f"{D}/EQA.dem_par") if ":" in l); w, h = int(par["width"]), int(par["nlines"])
    lat0, lon0, dlat, dlon = (float(par[k].split()[0]) for k in ("corner_lat", "corner_lon", "post_lat", "post_lon"))
    rd = lambda f: np.fromfile(f"{D}/{f}", dtype=np.float32).reshape(h, w)
    vel, U, msk = rd("vel"), rd("U.geo"), rd("mask") > 0; keep = msk if cfg["insar"]["mask"] else np.ones(vel.shape, bool)   # no mask by default (reproducibility)
    vup = np.where(keep, vel / U, np.nan).astype(np.float32); vup[~np.isfinite(vup)] = np.nan
    tf4326 = from_origin(lon0 - dlon / 2, lat0 - dlat / 2, dlon, -dlat)   # EQA corner is a pixel centre
    with rasterio.open(f"{D}/vel_up.tif", "w", driver="GTiff", height=h, width=w, count=1, dtype="float32", crs="EPSG:4326", transform=tf4326, nodata=np.nan) as o: o.write(vup, 1)
    cop, tf = read(f"{W}/cop30_nzvd2016.tif"); g = np.full(cop.shape, np.nan, np.float32)
    reproject(vup, g, src_transform=tf4326, src_crs="EPSG:4326", dst_transform=tf, dst_crs="EPSG:2193", resampling=Resampling.bilinear, src_nodata=np.nan, dst_nodata=np.nan)   # v104: bilinear (upsampling ~90 x 110 m pixels to 30 m; average gave pixel-block steps)
    m = cfg["masks"]; slope, _ = read(f"{W}/slope_deg.tif"); farmf, _ = read(f"{W}/farmfrac.tif"); near = read(f"{W}/near_organic.tif")[0].astype(bool)
    FARM = (farmf >= m["farmland_fraction_min"]) & (slope < m["slope_max_deg"])
    relv, _ = rel(g, FARM, FARM & ~near, cfg); write(f"{W}/insar_up.tif", g, tf); write(f"{W}/insar_rel.tif", relv, tf)
    M = blocks240(relv, cfg); write(f"{out}/maps/insar_rel_240m.tif", M, tf, 240.0)
    # class values: 240 m map at the ICESat-2 segment locations of each class (one value per location)
    S = pd.read_parquet(f"{work}/segments_classified.parquet").drop_duplicates(subset=["X", "Y"])
    F = S[S.land == "farmland"]; MIN = F[(F.order != "Organic") & (F.order != "not mapped") & ~F.near_peat]; hx0, hy0, hx1, hy1 = a["hauraki_box"]
    CL = {"Organic-Soil farmland": F[F.order == "Organic"], "mineral-soil farmland": MIN, "Hauraki Plains mineral-soil farmland": MIN[MIN.X.between(hx0, hx1) & MIN.Y.between(hy0, hy1)]}
    def samp(arr, d, res):
        r_ = ((tf.f - d.Y.values) // res).astype(int); c_ = ((d.X.values - tf.c) // res).astype(int); ok = (r_ >= 0) & (r_ < arr.shape[0]) & (c_ >= 0) & (c_ < arr.shape[1])
        v = np.full(len(d), np.nan); v[ok] = arr[r_[ok], c_[ok]]; return round(float(np.nanmedian(v)), 1)
    YY, XX = np.mgrid[0:g.shape[0], 0:g.shape[1]]; Xc = tf.c + (XX + .5) * 30; Yc = tf.f - (YY + .5) * 30; HBX = (Xc >= hx0) & (Xc < hx1) & (Yc >= hy0) & (Yc < hy1)
    org = read(f"{W}/organic.tif")[0].astype(bool)
    C30 = {"Organic-Soil farmland": FARM & org, "mineral-soil farmland": FARM & ~near, "Hauraki Plains mineral-soil farmland": FARM & ~near & HBX}
    import re as _re
    cnt = {k: int(v) for k, v in _re.findall(r"^(interferograms_input|acquisitions_input) (\d+)", log, _re.M)}
    from lib.maps import block_id, block_boot
    oc = FARM & org & np.isfinite(relv); G = cfg["maps"]["block_bootstrap_km"] * 1000
    org_ci = block_boot(relv[oc], block_id(relv.shape, tf, G, cfg["grids"]["square_origin"])[oc], cfg["maps"]["n_block_boot"])
    M30 = np.full(relv.shape, np.nan, np.float32); M30[:M.shape[0] * 8, :M.shape[1] * 8] = np.repeat(np.repeat(M, 8, 0), 8, 1)   # v96: grid effect alone (same cells)
    def maj(m):   # v97: 240 m cells where most valid 30 m cells are of the class
        ny, nx = M.shape; v = np.isfinite(relv[:ny * 8, :nx * 8]); k = (m[:ny * 8, :nx * 8] & v).reshape(ny, 8, nx, 8).sum((1, 3)); n = v.reshape(ny, 8, nx, 8).sum((1, 3))
        return np.isfinite(M) & (k * 2 > n)
    sites = {}
    if os.path.exists(f"{work}/repeat/sites_30m.tif"):
        sc_ = read(f"{work}/repeat/sites_30m.tif")[0]
        for i, nm in enumerate(("Torehape",), 1): sites[nm] = round(float(np.nanmedian(relv[sc_ == i])), 1)
    R = {"LiCSBAS network": cnt, "Organic-Soil farmland 95% CI (10 km block bootstrap)": org_ci, "repeat LiDAR sites (median of 30 m farmland cells)": sites, "step 15": log.strip(), "classes (median of 30 m farmland cells)": {k: round(float(np.nanmedian(relv[m])), 1) for k, m in C30.items()},
         "classes (240 m map at segment locations)": {k: samp(M, d, 240.0) for k, d in CL.items()},
         "classes (240 m map and 30 m cells, paired: cells with a 240 m value)": {k: {"30 m": round(float(np.nanmedian(relv[m & np.isfinite(M30)])), 1), "240 m": round(float(np.nanmedian(M30[m & np.isfinite(relv)])), 1), "share of cells": round(float((m & np.isfinite(M30) & np.isfinite(relv)).sum() / (m & np.isfinite(relv)).sum()), 3)} for k, m in C30.items()},   # v97: paired (v96 dropped cells of blocks below the fill rule from the 240 m side only)
         "classes (240 m cells, most of their valid 30 m cells in the class)": {k: round(float(np.nanmedian(M[maj(m)])), 1) for k, m in C30.items()},
         "zones_km2": {t: round(float((M <= t).sum() * 0.0576), 1) for t in cfg["insar"]["zones_mm_yr"]}}
    import geopandas as gpd
    from rasterio.features import geometry_mask
    from lib.maps import ref_adapt
    # Section 5.2: unreferenced InSAR at mineral-soil farmland segments minus built-up land within 3 km of Cambridge
    lo, la = cfg["insar"]["cambridge_centre"]; import pyproj
    cx, cy = pyproj.Transformer.from_crs(4326, 2193, always_xy=True).transform(lo, la)
    cam = S[(S.land == "built-up") & (np.hypot(S.X - cx, S.Y - cy) <= cfg["insar"]["cambridge_radius_km"] * 1000)]
    R["mineral-soil farmland relative to Cambridge built-up land"] = round(samp(g, MIN, 30.0) - samp(g, cam, 30.0), 1); R["Cambridge built-up segments"] = int(len(cam))
    # mask sensitivity (Section 5.3): the step-15 mask of Harvey et al. (2026) thresholds, applied to the same velocity
    mn = np.zeros(cop.shape, np.float32)   # mask applied pixel by pixel: a 30 m cell is removed when its InSAR pixel (nearest) is masked
    reproject(msk.astype(np.float32), mn, src_transform=tf4326, src_crs="EPSG:4326", dst_transform=tf, dst_crs="EPSG:2193", resampling=Resampling.nearest)
    gm = np.where(mn > 0.5, g, np.nan).astype(np.float32)
    relm, _ = rel(gm, FARM, FARM & ~near, cfg); Mm = blocks240(relm, cfg); of = FARM & org & np.isfinite(g)
    R["mask sensitivity"] = dict(Organic_farmland_cells_removed_pct=round(100 * float((of & ~np.isfinite(gm)).sum() / of.sum()), 1),
        Organic_change=round(float(np.nanmedian(relm[FARM & org]) - np.nanmedian(relv[FARM & org])), 1),
        zone_area_change_pct={t: round(100 * float(((Mm <= t).sum() / max((M <= t).sum(), 1)) - 1), 1) for t in cfg["insar"]["zones_mm_yr"]})
    farm_area = FARM & np.isfinite(g); R["masked share of study area (%)"] = round(100 * float(1 - np.isfinite(g[np.isfinite(cop)]).mean()), 1)
    json.dump(R, open(f"{out}/tables/insar.json", "w"), indent=1); print(json.dumps(R, indent=1), flush=True)
