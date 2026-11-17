"""Repeat LiDAR helpers: 96 m medians of 8 m differences, and per-cell capture dates from the GPS time of laser returns."""
import concurrent.futures as cf, glob, os, numpy as np, pandas as pd, rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin

GPS0 = np.datetime64("1980-01-06T00:00:00", "s")


def diff96(cfg, colA, colB, F=8, RES=96):
    """2024 minus earlier survey at 8 m, on the 8 m grid from the study-area corner, then the median of each 96 m cell (12 x 12)."""
    a = cfg["study_area"]; tf8 = from_origin(a["xmin"], a["ymax"], F, F); sh8 = (int((a["ymax"] - a["ymin"]) / F), int((a["xmax"] - a["xmin"]) / F))
    d = cfg["paths"]["linz_dem"]; names = sorted(set(os.listdir(f"{d}/{colA}")) & set(os.listdir(f"{d}/{colB}")))
    D = np.full(sh8, np.nan, np.float32)

    def one(n):
        out = []
        for c in (colA, colB):
            with rasterio.open(f"{d}/{c}/{n}") as r:
                x = r.read(1, out_shape=(r.height // F, r.width // F)).astype(np.float32)
                if r.nodata is not None: x[x == r.nodata] = np.nan
                out.append((x, r.transform * r.transform.scale(r.width / x.shape[1], r.height / x.shape[0])))
        if out[0][0].shape != out[1][0].shape: return None
        return out[1][0] - out[0][0], out[0][1]
    with cf.ThreadPoolExecutor(8) as ex:
        for res in ex.map(one, names):
            if res is None: continue
            dd, t = res; tmp = np.full(sh8, np.nan, np.float32)
            reproject(dd, tmp, src_transform=t, src_crs="EPSG:2193", dst_transform=tf8, dst_crs="EPSG:2193", resampling=Resampling.nearest, src_nodata=np.nan, dst_nodata=np.nan)
            m = np.isfinite(tmp) & ~np.isfinite(D); D[m] = tmp[m]
    B = RES // F; ny, nx = sh8[0] // B, sh8[1] // B
    v = D[:ny * B, :nx * B].reshape(ny, B, nx, B).transpose(0, 2, 1, 3).reshape(ny, nx, -1)
    with np.errstate(all="ignore"): med = np.nanmedian(v, axis=-1)
    return med.astype(np.float32), from_origin(a["xmin"], a["ymax"], RES, RES)


def _laz_cells(args):
    f, X0, Y1, RES, step = args
    import laspy
    try:
        with laspy.open(f) as r:
            xs, ys, ts = [], [], []
            for p in r.chunk_iterator(2_000_000): xs.append(np.asarray(p.x[::step])); ys.append(np.asarray(p.y[::step])); ts.append(np.asarray(p.gps_time[::step]))
        x, y, t = np.concatenate(xs), np.concatenate(ys), np.concatenate(ts)
        if t.max() <= 1e6: return None   # GPS week seconds: cannot be dated
        dt = GPS0 + (t + 1e9).astype("timedelta64[s]"); yr = 1970 + dt.astype(np.int64) / (365.2425 * 86400)
        df = pd.DataFrame(dict(ix=((x - X0) // RES).astype(int), iy=((Y1 - y) // RES).astype(int), yr=yr))
        return df.groupby(["ix", "iy"]).yr.agg(["size", "median"]).reset_index()
    except Exception as e:
        print("LAZ error", f, e, flush=True); return None


def capture_dates(cfg, dataset, cache, RES=96, step=20):
    """Median capture date (decimal year) of each 96 m cell from every 20th laser return of each LAZ tile; cached."""
    if os.path.exists(cache): c = pd.read_parquet(cache)
    else:
        a = cfg["study_area"]; fs = sorted(glob.glob(f"{cfg['paths']['laz']}/{dataset}/*.laz")); out = []
        with cf.ProcessPoolExecutor(12) as ex:
            for r in ex.map(_laz_cells, [(f, a["xmin"], a["ymax"], RES, step) for f in fs], chunksize=4):
                if r is not None: out.append(r)
        c = pd.concat(out, ignore_index=True); c = c.groupby(["ix", "iy"]).apply(lambda g: np.average(g["median"], weights=g["size"])).rename("t").reset_index()
        c.to_parquet(cache)
    return c


def cells_to_grid(c, shape):
    g = np.full(shape, np.nan, np.float32); ok = (c.iy >= 0) & (c.iy < shape[0]) & (c.ix >= 0) & (c.ix < shape[1]); g[c.iy[ok], c.ix[ok]] = c.t[ok]; return g
