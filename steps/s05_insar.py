"""Step 5 (Table S1): InSAR vertical velocity.
(a) LiCSBAS in WSL: steps 01-13 (batch_LiCSBAS.sh) unless a time series already exists in the job folder, then step 15 with the thresholds
    of config.yaml in the pipeline's own copy of the time series (the original job is never modified); step 16 is not used;
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
mkdir -p {T}; [ -f {T}/cum.h5 ] || cp $J/TS_GEOCml1clip/cum.h5 {T}/
rm -rf {T}/results {T}/info {T}/network; cp -r $J/TS_GEOCml1clip/results $J/TS_GEOCml1clip/info $J/TS_GEOCml1clip/network {T}/
cd {i['wsl_pipeline_dir']}
LiCSBAS15_mask_ts.py -t TS -c {s15['coh_thre']} -u {s15['n_unw_r_thre']} -v {s15['vstd_thre']} -T {s15['maxTlen_thre']} -g {s15['n_gap_thre']} -s {s15['stc_thre']} -i {s15['n_ifg_noloop_thre']} -l {s15['n_loop_err_thre']} -r {s15['resid_rms_thre']} | sed -n '/Noise index/,/Kept/p'
cp {T}/results/vel {T}/results/mask {wsl_path(D)}/; cp $J/GEOCml1clip/U.geo $J/GEOCml1clip/EQA.dem_par {wsl_path(D)}/
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
    vel, U, mask = rd("vel"), rd("U.geo"), rd("mask")
    vup = np.where(np.isfinite(mask) & (mask > 0), vel / U, np.nan).astype(np.float32); vup[~np.isfinite(vup)] = np.nan
    tf4326 = from_origin(lon0 - dlon / 2, lat0 - dlat / 2, dlon, -dlat)   # EQA corner is a pixel centre
    with rasterio.open(f"{D}/vel_up.tif", "w", driver="GTiff", height=h, width=w, count=1, dtype="float32", crs="EPSG:4326", transform=tf4326, nodata=np.nan) as o: o.write(vup, 1)
    cop, tf = read(f"{W}/cop30_nzvd2016.tif"); g = np.full(cop.shape, np.nan, np.float32)
    reproject(vup, g, src_transform=tf4326, src_crs="EPSG:4326", dst_transform=tf, dst_crs="EPSG:2193", resampling=Resampling.average, src_nodata=np.nan, dst_nodata=np.nan)
    m = cfg["masks"]; slope, _ = read(f"{W}/slope_deg.tif"); farmf, _ = read(f"{W}/farmfrac.tif"); near = read(f"{W}/near_organic.tif")[0].astype(bool)
    FARM = (farmf >= m["farmland_fraction_min"]) & (slope < m["slope_max_deg"])
    relv, _ = rel(g, FARM, FARM & ~near, cfg); write(f"{W}/insar_up.tif", g, tf); write(f"{W}/insar_rel.tif", relv, tf)
    M = blocks240(relv, cfg); write(f"{out}/maps/insar_rel_240m.tif", M, tf, 240.0)
    # class values: 240 m map at the ICESat-2 segment locations of each class (one value per location)
    S = pd.read_parquet(f"{work}/segments_classified.parquet").drop_duplicates(subset=["X", "Y"])
    F = S[S.land == "farmland"]; MIN = F[(F.order != "Organic") & (F.order != "not mapped") & ~F.near_peat]; hx0, hy0, hx1, hy1 = a["hauraki_box"]
    CL = {"Organic-Soil farmland": F[F.order == "Organic"], "mineral-soil farmland": MIN, "Hauraki Plains mineral-soil farmland": MIN[MIN.X.between(hx0, hx1) & MIN.Y.between(hy0, hy1)],
          "Kopuatai bog": S[(S.land == "natural wetland") & (S.order == "Organic") & (S.body == "Kopuatai")]}
    def samp(arr, d, res):
        r_ = ((tf.f - d.Y.values) // res).astype(int); c_ = ((d.X.values - tf.c) // res).astype(int); ok = (r_ >= 0) & (r_ < arr.shape[0]) & (c_ >= 0) & (c_ < arr.shape[1])
        v = np.full(len(d), np.nan); v[ok] = arr[r_[ok], c_[ok]]; return round(float(np.nanmedian(v)), 1)
    R = {"step 15": log.strip(), "classes (240 m map at segment locations)": {k: samp(M, d, 240.0) for k, d in CL.items()},
         "zones_km2": {t: round(float((M <= t).sum() * 0.0576), 1) for t in cfg["insar"]["zones_mm_yr"]}}
    farm_area = FARM & np.isfinite(g); R["masked share of study area (%)"] = round(100 * float(1 - np.isfinite(g[np.isfinite(cop)]).mean()), 1)
    json.dump(R, open(f"{out}/tables/insar.json", "w"), indent=1); print(json.dumps(R, indent=1), flush=True)
