"""Sampling helpers: LINZ LiDAR DEM heights and slopes at points (read from the public nz-elevation bucket), and WorldCover class with
the 3 x 3 'pure window' flag."""
import concurrent.futures as cf, numpy as np, pyproj, rasterio
from rasterio.windows import from_bounds
from scipy.ndimage import map_coordinates, minimum_filter

WCURL = "https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_{}_Map.tif"


def lidar_at(X, Y, tiles, base, factor=4, workers=12):
    """Bilinear height and terrain slope (degrees) at NZTM points from the 1 m DEM tiles read at 1/factor resolution (4 m for factor 4).
    tiles: [(name, [lon0, lat0, lon1, lat1]), ...] of one LINZ collection; base: its URL prefix."""
    out = np.full(X.size, np.nan); slope = np.full(X.size, np.nan)
    lon, lat = pyproj.Transformer.from_crs(2193, 4326, always_xy=True).transform(X, Y)
    jobs = [(n, np.flatnonzero((lon >= b[0]) & (lon <= b[2]) & (lat >= b[1]) & (lat <= b[3]))) for n, b, *_ in tiles]
    jobs = [j for j in jobs if j[1].size]

    def work(job):
        name, idx = job
        for k in range(4):
            try:
                with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", GDAL_HTTP_MAX_RETRY="5", GDAL_HTTP_RETRY_DELAY="2"), rasterio.open("/vsicurl/" + base + name) as d:
                    a = d.read(1, out_shape=(d.height // factor, d.width // factor)).astype(float); nod = d.nodata
                    tf = d.transform * d.transform.scale(d.width / a.shape[1], d.height / a.shape[0])
                break
            except Exception:
                if k == 3: return idx, None, None
        if nod is not None: a[a == nod] = np.nan
        c = (X[idx] - tf.c) / tf.a - .5; r = (Y[idx] - tf.f) / tf.e - .5
        v = map_coordinates(np.nan_to_num(a, nan=-9e9), [r, c], order=1, mode="constant", cval=-9e9); v[v < -1e8] = np.nan
        gy, gx = np.gradient(a); sl = np.degrees(np.arctan(np.hypot(gx, gy) / tf.a))
        return idx, v, map_coordinates(np.nan_to_num(sl, nan=90), [r, c], order=0, mode="nearest")
    with cf.ThreadPoolExecutor(workers) as ex:
        for idx, v, s in ex.map(work, jobs):
            if v is None: continue
            ok = np.isfinite(v) & np.isnan(out[idx]); out[idx[ok]] = v[ok]; slope[idx[ok]] = s[ok]
    return out, slope


def worldcover(lon, lat, tiles=("S39E174",), farm=(30, 40), built=50):
    """WorldCover 2021 class at each point, and whether all pixels of its 3 x 3 window are of the same cover (farmland or built-up)."""
    cls = np.zeros(lon.size, np.uint8); solid = np.zeros(lon.size, bool)
    for t in tiles:
        la0, lo0 = -int(t[1:3]), int(t[4:7]); m = (lat >= la0) & (lat < la0 + 3) & (lon >= lo0) & (lon < lo0 + 3)
        if not m.any(): continue
        with rasterio.open("/vsicurl/" + WCURL.format(t)) as d:
            win = from_bounds(lon[m].min() - .002, lat[m].min() - .002, lon[m].max() + .002, lat[m].max() + .002, d.transform).round_offsets().round_lengths()
            a = d.read(1, window=win); tf = d.window_transform(win)
        r = np.clip(np.floor((lat[m] - tf.f) / tf.e).astype(int), 0, a.shape[0] - 1); c = np.clip(np.floor((lon[m] - tf.c) / tf.a).astype(int), 0, a.shape[1] - 1)
        cls[m] = a[r, c]
        pf = minimum_filter(np.isin(a, farm).astype(np.uint8), 3); pb = minimum_filter((a == built).astype(np.uint8), 3)
        solid[m] = np.where(np.isin(cls[m], farm), pf[r, c], np.where(cls[m] == built, pb[r, c], 0)).astype(bool)
    return cls, solid
