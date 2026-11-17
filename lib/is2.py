"""ICESat-2 per-pass differencing and trend fitting (paper Section 3.1; Table S1 step 3).

A pass is one date and reference ground track (all beams pooled). Segments are reduced per pass and per LiDAR DEM ("col") to a median residual;
a target is differenced from its reference on the same pass and DEM; the trend is fitted by least squares with one constant per DEM and an
annual sinusoid. SE from the least-squares formula, 95% CI from the t distribution (df = passes - fitted parameters), and a bootstrap of
whole passes as a check."""
import numpy as np, pandas as pd
from scipy.stats import t as tdist

T0 = 2021.5   # time origin of the trend term (does not change the rate)


def design(pp, t_break=None, seasonal=True):
    """Columns: one constant per DEM, trend, sin, cos [, slope change after t_break, step at t_break]. Returns (A, index of the trend column)."""
    G = pd.get_dummies(pp["col"]).astype(float).values; t = pp["t"].values
    cols = [G, t - T0] + ([np.sin(2 * np.pi * t), np.cos(2 * np.pi * t)] if seasonal else [])
    if t_break is not None:
        after = (t >= t_break).astype(float); cols += [after * (t - t_break), after]
    return np.column_stack(cols), G.shape[1]


def per_pass(d, min_seg):
    """Median residual per pass and DEM, kept where at least min_seg segments."""
    pp = d.groupby(["pass", "col"]).agg(t=("t", "mean"), r=("r", "median"), n=("r", "size")).reset_index()
    return pp[pp.n >= min_seg]


def fit(pp, min_passes=8, n_boot=2000, seed=0, seasonal=True):
    """Trend (mm/yr) of per-pass values pp (columns pass, col, t, r in m). None if fewer than min_passes passes."""
    if len(pp) < min_passes or pp["pass"].nunique() < min_passes: return None
    A, k = design(pp, seasonal=seasonal); y = pp["r"].values
    b, *_ = np.linalg.lstsq(A, y, rcond=None); e = y - A @ b
    npass = int(pp["pass"].nunique()); dof = npass - A.shape[1]
    if dof < 1: return None
    se = float(np.sqrt((e @ e) / (len(y) - A.shape[1]) * np.linalg.pinv(A.T @ A)[k, k])) * 1000
    rate = float(b[k]) * 1000; h = float(tdist.ppf(0.975, dof)) * se
    out = dict(rate=round(rate, 1), se=round(se, 2), ci=[round(rate - h, 1), round(rate + h, 1)], dof=dof, passes=npass, rows=len(pp))
    if n_boot:
        passes = pp["pass"].unique(); rng = np.random.default_rng(seed); grp = {p: np.flatnonzero(pp["pass"].values == p) for p in passes}; bs = []
        for _ in range(n_boot):
            i = np.concatenate([grp[p] for p in rng.choice(passes, len(passes))]); Ai, _ = design(pp.iloc[i], seasonal=seasonal)
            if Ai.shape[1] != A.shape[1]: continue
            bs.append(np.linalg.lstsq(Ai, y[i], rcond=None)[0][k])
        bs = np.array(bs) * 1000
        out.update(boot_se=round(float(bs.std(ddof=1)), 2), boot_ci=[round(float(np.percentile(bs, 2.5)), 1), round(float(np.percentile(bs, 97.5)), 1)])
    return out


def diff_series(target, ref, min_seg):
    """Per-pass target minus reference on the same pass and DEM."""
    j = per_pass(target, min_seg).merge(per_pass(ref, min_seg), on=["pass", "col"], suffixes=("", "_ref"))
    j["r"] = j["r"] - j["r_ref"]; j["t"] = j["t"]
    return j[["pass", "col", "t", "r", "n", "n_ref"]]


def trend(d, min_seg=10, **k): return fit(per_pass(d, min_seg), **k)
def diff_trend(target, ref, min_seg=10, **k): return fit(diff_series(target, ref, min_seg), **k)


def period(pp, t0, t1): return pp[(pp["t"] >= t0) & (pp["t"] < t1)]


def segmented(pp, t_break, n_boot=2000, seed=5):
    """Shared DEM constants and annual cycle, separate slopes before and after t_break, and a step at it.
    95% CIs of the two rates, their difference and the step from a bootstrap of whole passes."""
    def one(p):
        A, k = design(p, t_break); b = np.linalg.lstsq(A, p["r"].values, rcond=None)[0]; return b[k] * 1000, (b[k] + b[-2]) * 1000, b[-1] * 1000
    A, k = design(pp, t_break); b0, a0, s0 = one(pp)
    passes = pp["pass"].unique(); rng = np.random.default_rng(seed); grp = {p: np.flatnonzero(pp["pass"].values == p) for p in passes}; bs = []
    for _ in range(n_boot):
        p = pp.iloc[np.concatenate([grp[x] for x in rng.choice(passes, len(passes))])]
        if design(p, t_break)[0].shape[1] != A.shape[1]: continue
        bs.append(one(p))
    bs = np.array(bs); ci = lambda v: [round(float(np.percentile(v, 2.5)), 1), round(float(np.percentile(v, 97.5)), 1)]
    return dict(before=round(float(b0), 1), before_CI=ci(bs[:, 0]), after=round(float(a0), 1), after_CI=ci(bs[:, 1]), after_minus_before=round(float(a0 - b0), 1),
                diff_CI=ci(bs[:, 1] - bs[:, 0]), step_mm=round(float(s0), 1), step_CI=ci(bs[:, 2]), passes=int(len(passes)))


def amplitude(pp, n_boot=500, seed=3):
    """Seasonal amplitude (mm): half the peak-to-trough range of the fitted sinusoid, with a 95% CI from resampled whole passes."""
    def amp(p):
        A, k = design(p); b = np.linalg.lstsq(A, p["r"].values, rcond=None)[0]; return float(np.hypot(b[k + 1], b[k + 2])) * 1000
    a = amp(pp); passes = pp["pass"].unique(); rng = np.random.default_rng(seed); grp = {p: np.flatnonzero(pp["pass"].values == p) for p in passes}
    ncol = design(pp)[0].shape[1]; bs = []
    for _ in range(n_boot):
        p = pp.iloc[np.concatenate([grp[x] for x in rng.choice(passes, len(passes))])]
        if design(p)[0].shape[1] == ncol: bs.append(amp(p))
    return dict(amplitude_mm=round(a, 1), ci=[round(float(np.percentile(bs, 2.5)), 1), round(float(np.percentile(bs, 97.5)), 1)])


def separate_difference(pp, t_break, n_boot=2000, seed=5):
    """Rate before minus rate after t_break, each period fitted separately (own DEM constants and annual cycle); 95% CI from a bootstrap of whole passes."""
    def one(p):
        a, b = p[p.t < t_break], p[p.t >= t_break]
        ra = np.linalg.lstsq(*(lambda A, k: (A, a["r"].values))(*design(a)), rcond=None)[0][design(a)[1]]
        rb = np.linalg.lstsq(*(lambda A, k: (A, b["r"].values))(*design(b)), rcond=None)[0][design(b)[1]]
        return (ra - rb) * 1000
    d0 = one(pp); passes = pp["pass"].unique(); rng = np.random.default_rng(seed); grp = {p: np.flatnonzero(pp["pass"].values == p) for p in passes}; bs = []
    for _ in range(n_boot):
        p = pp.iloc[np.concatenate([grp[x] for x in rng.choice(passes, len(passes))])]
        if (p.t < t_break).sum() > 5 and (p.t >= t_break).sum() > 5: bs.append(one(p))
    return dict(before_minus_after=round(float(d0), 1), CI=[round(float(np.percentile(bs, 2.5)), 1), round(float(np.percentile(bs, 97.5)), 1)])
