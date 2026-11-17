"""Map helpers on the 30 m analysis grid: local mineral-soil reference (Section 3.2), 240 m grid, 10 km block bootstrap."""
import numpy as np, pandas as pd, rasterio
from scipy.ndimage import uniform_filter

R = 30.0


def read(p):
    with rasterio.open(p) as d: return d.read(1), d.transform


def local_ref(a, mask, radius_m, cap, min_cells):
    """Mean of a (capped at +-cap) over mask cells in a square window +-radius_m around each cell; NaN where fewer than min_cells."""
    ok = mask & np.isfinite(a); v = np.where(ok, np.clip(a, -cap, cap), 0.0); w = ok.astype(np.float32)
    k = int(2 * radius_m / R) | 1; s = uniform_filter(v, k, mode="constant"); n = uniform_filter(w, k, mode="constant")
    return np.where(n * k * k >= min_cells, s / np.maximum(n, 1e-9), np.nan)


def ref_adapt(a, mask, cfg):
    """+-3 km window, widened to +-6 then +-10 km where it has too few reference cells."""
    m = cfg["maps"]; r = None
    for km in m["ref_window_km"]:
        x = local_ref(a, mask, km * 1000, m["ref_cap_mm_yr"], m["ref_min_cells"]); r = x if r is None else np.where(np.isfinite(r), r, x)
    return r


def rel(a, cells, refmask, cfg):
    """Each cell of `cells` minus its local reference; cells beyond +-cap are dropped."""
    ref = ref_adapt(a, refmask, cfg); cap = cfg["maps"]["ref_cap_mm_yr"]
    return np.where(cells & np.isfinite(a) & (np.abs(a) < cap) & np.isfinite(ref), a - ref, np.nan), ref


def blocks240(a, cfg, B=8):
    """240 m grid: median of the 30 m cells of each 8 x 8 block where at least 25% are valid."""
    ny, nx = (a.shape[0] // B) * B, (a.shape[1] // B) * B
    v = a[:ny, :nx].reshape(ny // B, B, nx // B, B).transpose(0, 2, 1, 3).reshape(ny // B, nx // B, -1); n = np.isfinite(v).sum(-1)
    with np.errstate(all="ignore"): return np.nanmedian(np.where(n[..., None] >= cfg["grids"]["map_cell_min_fill"] * B * B, v, np.nan), axis=-1)


def block_id(shape, tf, size_m, origin=(0.0, 0.0)):
    """Square id of each cell; squares of size_m counted from origin (the study-area corner, Table S1)."""
    YY, XX = np.mgrid[0:shape[0], 0:shape[1]]; xc = tf.c + (XX + .5) * tf.a; yc = tf.f + (YY + .5) * tf.e
    return ((xc - origin[0]) // size_m).astype(np.int64) * 1000 + ((yc - origin[1]) // size_m).astype(np.int64)


def block_boot(vals, blocks, n_boot, stat=np.median, seed=0):
    rng = np.random.default_rng(seed); u = np.unique(blocks); idx = {b: np.flatnonzero(blocks == b) for b in u}
    res = [stat(vals[np.concatenate([idx[b] for b in rng.choice(u, len(u))])]) for _ in range(n_boot)]
    return [round(float(np.percentile(res, 2.5)), 1), round(float(np.percentile(res, 97.5)), 1)]


def bias_grid(csv, regional_m, shape, tf, size_m, origin=(0.0, 0.0)):
    """Survey bias (m) per 10 km square from step 3; squares without a value take the regional value."""
    b = pd.read_csv(csv, index_col=0)["bias_m"]; bid = block_id(shape, tf, size_m, origin)
    return pd.Series(bid.ravel()).map(b).fillna(regional_m).values.reshape(shape).astype(np.float32)
