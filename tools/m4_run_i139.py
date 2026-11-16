"""I139 M4 (inflection point, Q4) for papers that report a threshold (spec field inflection_point), prereg Section 6 M4.
Usage: python -I tools/m4_run_i139.py <rec_id> [...]  -> extraction/m4/<rec_id>.json
Declared deviation (budget): in the simulations (d), (e) the CI of k is the weighted profile-likelihood interval
(weights rescaled to the Kish effective n, chi-square(1) cut 3.84) on a 50-point grid, not a Muggeo bootstrap B = 200;
real data keep the 100-point grid and Rao-Wu bootstrap B = 500 as registered. PSU random effect SD 0.2 (logit scale).
"""
import os, sys, re, json, hashlib, time, traceback
import numpy as np
from scipy import stats
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m23_common_i139 as C

B_REAL, NSIM, GRID_REAL, GRID_BOOT, GRID_SIM, RE_SD = 500, 100, 100, 25, 25, 0.2
OUT = os.path.join(C.X, "m4")
KSCALE = {"8991f10763": 0.01}   # builder exposure is AISI / 100 (paper's k in raw AISI units)
os.makedirs(OUT, exist_ok=True)


def rng_for(s):
    return np.random.default_rng(int(hashlib.sha256(s.encode()).hexdigest()[:15], 16))


def wq(v, w, q):
    o = np.argsort(v); cw = np.cumsum(w[o]) / w.sum()
    return np.interp(q, cw, v[o])


def loglik(X, y, w, beta):
    eta = np.clip(X @ beta, -30, 30)
    return float(np.sum(w * (y * eta - np.log1p(np.exp(eta)))))


def fit(X, y, w, b0=None):
    try:
        beta, _ = C.irls(X, y, w, beta0=b0, it=30, tol=1e-7)
        return beta
    except np.linalg.LinAlgError:   # fallback added 15/11: ridge 1e-6 IRLS (same as m3_desc)
        b = np.zeros(X.shape[1]) if b0 is None else np.array(b0, float)
        for _ in range(50):
            p = 1 / (1 + np.exp(-np.clip(X @ b, -30, 30)))
            g = X.T @ (w * (y - p)) - 1e-6 * b
            H = (X * (w * p * (1 - p))[:, None]).T @ X + 1e-6 * np.eye(X.shape[1])
            step = np.linalg.lstsq(H, g, rcond=None)[0]
            b = b + step
            if np.max(np.abs(step)) < 1e-7:
                break
        return b


def profile(Z, x, y, w, grid, b0=None):
    """Z: covariate matrix incl. intercept and x; returns loglik per k (weights as given, not normalised)."""
    ll, betas, b = [], [], b0
    wn = w / w[w > 0].mean()
    for k in grid:
        X = np.column_stack([Z, np.maximum(x - k, 0.0)])
        b = fit(X, y, w, None if b is None or len(b) != X.shape[1] else b)
        ll.append(loglik(X, y, wn, b)); betas.append(b)
    return np.array(ll), betas


def run(rid, T, L):
    t0 = time.time()
    spec = T[rid]
    b = C.load_builder(rid)
    if b.get("not_testable"):
        return {"rec_id": rid, "status": "not_testable"}
    df, y = b["df"], b["y"]
    dom = b["domain"].fillna(False).astype(bool)
    terms, tgt = list(b["terms"]), b["target"]
    stem0 = re.split(r"_q\d|__|_t\d|_seg", tgt)[0]
    xcol = next((c for c in ["x", stem0, "pni", "cally", "hrr"] if c in df.columns and df.loc[dom, c].nunique() > 50), None)
    if xcol is None:
        return {"rec_id": rid, "status": "no_continuous_index", "target": tgt}
    stem = re.split(r"_q\d|__|_t\d|_seg", tgt)[0]
    expo = {t for t in terms if t == tgt or t == xcol or t.startswith(stem + "_") or t.startswith(xcol + "_") or t.startswith("x_") or t.startswith("seg")}
    covs = [t for t in terms if t not in expo]
    w, kind = C.correct_weight(rid, b)
    m = dom & df[[y, xcol] + covs].notna().all(axis=1) & (w > 0)
    mm = m.to_numpy()
    x = df.loc[mm, xcol].to_numpy(float); yy = df.loc[mm, y].to_numpy(float); ww = w[mm].to_numpy(float)
    covs = [t for t in covs if df.loc[m, t].nunique() > 1]
    Z = np.column_stack([np.ones(mm.sum()), x, df.loc[mm, covs].to_numpy(float)])
    # drop covariate columns that are linear combinations of earlier ones (added 15/11 after a singular-matrix error)
    keep = [0, 1]
    for j in range(2, Z.shape[1]):
        if np.linalg.matrix_rank(Z[:, keep + [j]]) == len(keep) + 1:
            keep.append(j)
    dropped = [covs[j - 2] for j in range(2, Z.shape[1]) if j not in keep]
    covs = [covs[j - 2] for j in keep[2:]]
    Z = Z[:, keep]
    sdx = float(np.sqrt(np.average((x - np.average(x, weights=ww)) ** 2, weights=ww)))
    lo_g, hi_g = wq(x, ww, 0.01), wq(x, ww, 0.99)
    grid = np.linspace(lo_g, hi_g, GRID_REAL)
    k_rep = float(str(spec["inflection_point"]).split(",")[0]) * KSCALE.get(rid, 1.0)
    res = {"rec_id": rid, "status": "run", "n": int(mm.sum()), "xcol": xcol, "sd_x": sdx, "k_reported": k_rep,
           "grid": [float(lo_g), float(hi_g)], "k_reported_in_range": bool(lo_g <= k_rep <= hi_g), "covariates": covs, "collinear_dropped": dropped}
    # (a) profile k
    ll, betas = profile(Z, x, yy, ww, grid)
    ik = int(np.argmax(ll)); khat = float(grid[ik])
    res["k_hat"] = khat
    # (c) RCS nonlinearity (design-based Wald via svy_logit)
    knots = [wq(x, ww, q) for q in (0.05, 0.35, 0.65, 0.95)]
    def rcs(v, kn):
        t = kn; out = []
        for j in range(len(t) - 2):
            out.append(np.maximum(v - t[j], 0) ** 3 - np.maximum(v - t[-2], 0) ** 3 * (t[-1] - t[j]) / (t[-1] - t[-2])
                       + np.maximum(v - t[-1], 0) ** 3 * (t[-2] - t[j]) / (t[-1] - t[-2]))
        return out
    sc = (knots[-1] - knots[0]) ** 2
    for j, s in enumerate(rcs(df[xcol].to_numpy(float), knots)):
        df[f"_rcs{j}"] = s / sc
    fitn = C.nc.svy_logit(df, y, [xcol, "_rcs0", "_rcs1"] + covs, w.rename("_w") if False else pd_w(df, w), domain=m)
    idx = [fitn["names"].index("_rcs0"), fitn["names"].index("_rcs1")]
    bb = fitn["coef"][idx]; V = fitn["V"][np.ix_(idx, idx)]
    wald = float(bb @ np.linalg.solve(V, bb))
    res["rcs_wald"] = wald; res["rcs_p"] = float(stats.chi2.sf(wald, 2)); res["rcs_F_p"] = float(stats.f.sf(wald / 2, 2, fitn["df"]))
    # (b) Rao-Wu bootstrap of k
    keys, Mx = C.raowu_multipliers(df[df.SDMVSTRA.notna()], B_REAL, f"20261114-i139-m4:{rid}")
    kpos = {k: i for i, k in enumerate(keys)}
    rk = df.SDMVSTRA.fillna(-1).to_numpy(np.int64)[mm] * 1000 + df.SDMVPSU.fillna(-1).to_numpy(np.int64)[mm]
    ki = np.array([kpos[v] for v in rk])
    kb = []
    gb = np.linspace(lo_g, hi_g, GRID_BOOT)
    for r in range(B_REAL):
        llb, _ = profile(Z, x, yy, ww * Mx[r, ki], gb, b0=betas[ik])
        kb.append(float(gb[int(np.argmax(llb))]))
    kb = np.array(kb)
    lo, hi = np.percentile(kb, 2.5), np.percentile(kb, 97.5)
    res["k_ci"] = [float(lo), float(hi)]; res["k_ci_width_sd"] = float((hi - lo) / sdx)
    res["boot_edge_share"] = float(np.mean((kb == gb[0]) | (kb == gb[-1])))
    # simulations: true model = two-segment at k_true (reported k if in range, else k_hat), slope change from fit there
    k_true = k_rep if res["k_reported_in_range"] else khat
    Xt = np.column_stack([Z, np.maximum(x - k_true, 0.0)])
    bt = fit(Xt, yy, ww)
    psu = rk
    upsu, inv = np.unique(psu, return_inverse=True)
    neff = ww.sum() ** 2 / (ww ** 2).sum()
    wsim = ww / ww.sum() * neff
    gs = np.linspace(lo_g, hi_g, GRID_SIM)
    crit = stats.chi2.ppf(0.95, 1) / 2

    def sim_ci(eta, seed):
        g = rng_for(seed)
        re_ = g.normal(0, RE_SD, len(upsu))[inv]
        ys = (g.random(len(eta)) < 1 / (1 + np.exp(-(eta + re_)))).astype(float)
        llp, _ = profile(Z, x, ys, wsim, gs)
        llp = llp * 1.0
        ok = gs[llp >= llp.max() - crit]
        return float(ok.min()), float(ok.max()), float(gs[int(np.argmax(llp))])
    out = {}
    for lev, f in (("paper", 1.0), ("half", 0.5)):
        bt2 = bt.copy(); bt2[-1] = bt[-1] * f
        eta = Xt @ bt2
        cover = loc = 0
        for s in range(NSIM):
            a, c, _ = sim_ci(eta, f"20261114-i139-m4sim:{rid}:{lev}:{s}")
            cover += a <= k_true <= c
            loc += (c - a) < sdx and a <= k_true <= c
        out[lev] = {"coverage": cover / NSIM, "power_localise": loc / NSIM, "slope_change": float(bt2[-1])}
    # (e) negative control: smooth spline model
    Xs = np.column_stack([Z, df.loc[mm, "_rcs0"].to_numpy(float), df.loc[mm, "_rcs1"].to_numpy(float)])
    bs_ = fit(Xs, yy, ww)
    eta_s = Xs @ bs_
    fl = 0
    for s in range(NSIM):
        a, c, _ = sim_ci(eta_s, f"20261114-i139-m4neg:{rid}:{s}")
        fl += (c - a) < sdx
    out["neg_false_localisation"] = fl / NSIM
    res["sim"] = out
    res["k_true_used"] = k_true
    valid = out["paper"]["coverage"] >= 0.90 and out["neg_false_localisation"] <= 0.10
    width = res["k_ci_width_sd"]
    if not valid:
        cl = "classification_invalid"
    elif width < 1.0:
        cl = ("localisable_matching" if lo <= k_rep <= hi else "localisable_different") if res["k_reported_in_range"] else "localisable_reported_k_out_of_range"
    elif out["half"]["power_localise"] >= 0.80:
        cl = "not_localisable"
    else:
        cl = "inconclusive"
    res["class"] = cl
    res["class_0_5sd"] = ("localisable" if (hi - lo) / sdx < 0.5 else "not_localisable_or_inconclusive") if valid else cl
    res["seconds"] = int(time.time() - t0)
    return res


def pd_w(df, w):
    df["_wm4"] = w
    return "_wm4"


def main():
    T, L = C.spec_table(), C.labels()
    for rid in sys.argv[1:]:
        try:
            res = run(rid, T, L)
        except Exception as e:
            traceback.print_exc()
            res = {"rec_id": rid, "status": "error", "error": f"{type(e).__name__}: {e}"[:300]}
        json.dump(res, open(os.path.join(OUT, rid + ".json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1, default=float)
        print(rid, res["status"], res.get("class"), "k_rep", res.get("k_reported"), "k_hat", res.get("k_hat"), "ci", res.get("k_ci"),
              "w/sd", res.get("k_ci_width_sd"), res.get("seconds"), flush=True)


if __name__ == "__main__":
    main()
