"""Step 1 (Table S1): download and grid the input data.
(a) LINZ 1 m DEM tiles of each collection intersecting the study area (downloaded once from the public nz-elevation bucket);
(b) each collection averaged to the 30 m analysis grid (exact block means, tile by tile);
(c) Copernicus GLO-30 converted from EGM2008 to NZVD2016 (NZGeoid2016) and reprojected once (bilinear) to the 30 m grid;
(d) 30 m masks: WorldCover farmland, grassland and wetland fractions, slope (Copernicus), Organic Soils and the 500 m buffer zone.
Outputs: data/work/grid30/*.tif (NZTM 2000, 30 m, origin at the study-area corner)."""
import glob, json, os, urllib.request, numpy as np, rasterio, geopandas as gpd, pyproj
from rasterio.transform import from_origin
from rasterio.warp import reproject, Resampling
from rasterio.features import geometry_mask
from rasterio.merge import merge
from scipy import ndimage

R30 = 30.0


def grid(cfg, res=R30):
    a = cfg["study_area"]; return from_origin(a["xmin"], a["ymax"], res, res), (int((a["ymax"] - a["ymin"]) / res), int((a["xmax"] - a["xmin"]) / res))


def write(p, a, tf, nodata=np.nan):
    with rasterio.open(p, "w", driver="GTiff", height=a.shape[0], width=a.shape[1], count=1, dtype=a.dtype, crs="EPSG:2193", transform=tf, nodata=nodata,
                       compress="deflate", tiled=True, BIGTIFF="IF_SAFER") as d: d.write(a, 1)


def download_tiles(cfg, col):
    """Tiles of one LINZ collection that intersect the study area; skipped if present."""
    TI = json.load(open(f"{cfg['paths']['raw']}/tile_index.json")); a = cfg["study_area"]; d = f"{cfg['paths']['linz_dem']}/{col}"; os.makedirs(d, exist_ok=True)
    tr = pyproj.Transformer.from_crs(2193, 4326, always_xy=True); lo0, la0 = tr.transform(a["xmin"], a["ymin"]); lo1, la1 = tr.transform(a["xmax"], a["ymax"])
    for name, b, *_ in TI[col]["tiles"]:
        if b[2] < lo0 or b[0] > lo1 or b[3] < la0 or b[1] > la1: continue
        f = f"{d}/{name}"
        if not os.path.exists(f): urllib.request.urlretrieve(TI[col]["base"] + name, f)
    return sorted(glob.glob(f"{d}/*.tif*"))


def dem30(cfg, col, out):
    tf, sh = grid(cfg); dst = np.full(sh, np.nan, np.float32)
    for f in download_tiles(cfg, col):
        with rasterio.open(f) as d:
            tmp = np.full(sh, np.nan, np.float32)
            reproject(rasterio.band(d, 1), tmp, src_nodata=d.nodata if d.nodata is not None else -9999, dst_nodata=np.nan, dst_transform=tf, dst_crs="EPSG:2193",
                      resampling=Resampling.average, init_dest_nodata=False)
            m = np.isfinite(tmp) & ~np.isfinite(dst); dst[m] = tmp[m]
    write(out, dst, tf); print("wrote", out, f"{np.isfinite(dst).mean():.3f} filled", flush=True)


def copernicus(cfg, out):
    raw = cfg["paths"]["raw"]; tf, sh = grid(cfg)
    srcs = [rasterio.open(f"{raw}/copernicus/Copernicus_DSM_COG_10_{t}_DEM.tif") for t in cfg["copernicus"]["tiles"]]; mos, mtr = merge(srcs); crs = srcs[0].crs
    mos = mos[0].astype(np.float64); h, w = mos.shape; LON, LAT = np.meshgrid(mtr.c + (np.arange(w) + .5) * mtr.a, mtr.f + (np.arange(h) + .5) * mtr.e)
    def samp(p):
        with rasterio.open(p) as g:
            arr = g.read(1).astype(np.float64); t = g.transform
            return ndimage.map_coordinates(arr, [(LAT - t.f) / t.e - .5, (LON - t.c) / t.a - .5], order=1, mode="nearest")
    Hnz = mos + samp(f"{raw}/proj/us_nga_egm08_25.tif") - samp(f"{raw}/proj/nz_linz_nzgeoid2016.tif")   # EGM2008 -> ellipsoid -> NZVD2016
    dst = np.full(sh, np.nan, np.float32)
    reproject(Hnz.astype(np.float32), dst, src_transform=mtr, src_crs=crs, dst_transform=tf, dst_crs="EPSG:2193", resampling=Resampling.bilinear, dst_nodata=np.nan)
    write(out, dst, tf); print("wrote", out, flush=True)


def masks(cfg, W):
    raw = cfg["paths"]["raw"]; tf, sh = grid(cfg); m = cfg["masks"]
    with rasterio.open(f"{raw}/worldcover/ESA_WorldCover_10m_2021_v200_S39E174_Map.tif") as w:
        win = rasterio.windows.from_bounds(174.85, -38.2, 175.8, -36.95, w.transform).round_offsets().round_lengths(); s_ = w.read(1, window=win); wt = w.window_transform(win)
    for name, cls in (("farmfrac", m["worldcover_farmland"]), ("grassfrac", [30]), ("wetfrac", [90]), ("builtfrac", [m["worldcover_builtup"]])):
        a = np.zeros(sh, np.float32)
        reproject(np.isin(s_, cls).astype(np.float32), a, src_transform=wt, src_crs="EPSG:4326", dst_transform=tf, dst_crs="EPSG:2193", resampling=Resampling.average)
        write(f"{W}/{name}.tif", a, tf)
    with rasterio.open(f"{W}/cop30_nzvd2016.tif") as c: dem = c.read(1)
    gy, gx = np.gradient(np.nan_to_num(dem, nan=0), R30); write(f"{W}/slope_deg.tif", np.degrees(np.arctan(np.hypot(gx, gy))).astype(np.float32), tf)
    org = gpd.read_file(f"{raw}/soils/smap_organic_waikato.shp").to_crs(2193)
    write(f"{W}/organic.tif", (~geometry_mask(org.geometry, sh, tf)).astype(np.uint8), tf, nodata=None)
    write(f"{W}/near_organic.tif", (~geometry_mask(org.geometry.buffer(m["buffer_m"]), sh, tf)).astype(np.uint8), tf, nodata=None)
    print("masks written", flush=True)


def run(cfg):
    W = f"{cfg['paths']['work']}/grid30"; os.makedirs(W, exist_ok=True)
    for col in cfg["lidar"]["waikato"]:
        if not os.path.exists(f"{W}/dem30_{col}.tif"): dem30(cfg, col, f"{W}/dem30_{col}.tif")
    if not os.path.exists(f"{W}/cop30_nzvd2016.tif"): copernicus(cfg, f"{W}/cop30_nzvd2016.tif")
    masks(cfg, W)
