"""I139 M3 descriptive analyses (prereg Section 6 M3, "descriptive only"; not used for the main conclusion):
E3b (drop K-C/K-M/K-O/K-X terms), prevalence ratio (survey-weighted Poisson), marginal standardised OR, RoE
(vibration of effects over adjustment subsets), conventional set, adding diabetes or antihypertensive drugs.
Usage: python -I tools/m3_desc_i139.py <rec_id> [...]  -> extraction/m3_desc/<rec_id>.json
Bootstrap here B = 200 (descriptive). Term -> concept mapping by name regex, logged per paper.
"""
import os, sys, re, json, hashlib, time, itertools, traceback
import numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m23_common_i139 as C
from m3_run_i139 import cls, concepts

B = 200
OUT = os.path.join(C.X, "m3_desc")
os.makedirs(OUT, exist_ok=True)
TERM_RX = [("kidney_measure", r"egfr|^scr|creat|cystatin|uacr|^acr|bun|urea|^ucr"),
           ("uric_acid", r"^ua$|^ua_|uric|^sua|hyperua|^hua"),
           ("bp_measured", r"^sbp|^dbp|^bp_mean|systolic|diastolic|^map$"),
           ("glycaemia_measured", r"glucose|hba1c|insulin|homa|^fpg|^glu$|^ghb"),
           ("lipids", r"^tc$|^tc_|hdl|ldl|^tg$|^tg_|trig|chol|nonhdl|lipid_lvl"),
           ("inflammation", r"crp|^wbc|neut|lymph|^plt|platelet|^mono"),
           ("cvd", r"cvd|^chd|^hf$|^hf_|heart|stroke|^chf|angina|^mi$"),
           ("bmi_waist", r"bmi|waist|^wc$|whtr"),
           ("diabetes_dx", r"diab|^dm$|^dm_|^dm3|c_dm|dx_diabetes"),
           ("hypertension_dx", r"htn|hypert|dx_hypertension"),
           ("albumin", r"albumin|^alb$|^alb_|hypoalb"),
           ("hemoglobin", r"hemoglobin|^hgb|^hb$|anemia|anaemia"),
           ("liver", r"^alt|^ast|^ggt|bilirub")]
CONFOUNDER_RX = r"^age|female|^male|^sex|race|^nhw|white|black|hisp|mexican|educ|^edu|^pir|poverty|income|smok|smk|alc|drink|^c_pa|activity|^pa$|^pa_|inactive|marri|cycle"


def term_concept(t):
    base = t.split("__")[0].lower()
    if re.search(CONFOUNDER_RX, base):
        return "confounder_like"
    for c, rx in TERM_RX:
        if re.search(rx, base):
            return c
    return "other"


def irls_pois(Xm, y, w, beta0=None, it=50, tol=1e-9):
    w = w / w[w > 0].mean()
    beta = np.zeros(Xm.shape[1]) if beta0 is None else beta0.copy()
    if beta0 is None:
        beta[0] = np.log(max(np.average(y, weights=w), 1e-6))
    for _ in range(it):
        mu = np.exp(np.clip(Xm @ beta, -30, 5))
        H = Xm.T @ (Xm * (w * mu)[:, None])
        g = Xm.T @ (w * (y - mu))
        step = np.linalg.solve(H + 1e-10 * np.eye(len(beta)), g)
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            break
    return beta


def irls_ridge(Xm, y, w, it=40, tol=1e-7, lam=1e-6):
    """Weighted logistic IRLS with a tiny ridge (RoE only), None when it does not converge."""
    w = w / w[w > 0].mean()
    beta = np.zeros(Xm.shape[1])
    for _ in range(it):
        p = 1 / (1 + np.exp(-np.clip(Xm @ beta, -30, 30)))
        H = Xm.T @ (Xm * (w * p * (1 - p))[:, None]) + lam * len(y) * np.eye(Xm.shape[1])
        step = np.linalg.solve(H, Xm.T @ (w * (y - p)) - lam * len(y) * beta)
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            return beta
    return None


def marginal_or(Xm, beta, w, j, binary):
    X0, X1 = Xm.copy(), Xm.copy()
    if binary:
        X0[:, j], X1[:, j] = 0.0, 1.0
    else:
        X1[:, j] = X1[:, j] + 1.0
    p0 = np.average(1 / (1 + np.exp(-(X0 @ beta))), weights=w)
    p1 = np.average(1 / (1 + np.exp(-(X1 @ beta))), weights=w)
    return float(np.log(p1 / (1 - p1)) - np.log(p0 / (1 - p0)))


def run(rid, T, L):
    t0 = time.time()
    lab = L[rid]
    b = C.load_builder(rid)
    if b.get("not_testable"):
        return {"rec_id": rid, "status": "not_testable"}
    df, y = b["df"], b["y"]
    dom = b["domain"].fillna(False).astype(bool)
    terms, tgt = list(b["terms"]), b["target"]
    w, kind = C.correct_weight(rid, b)
    df["_wc"] = w
    stem = re.split(r"_q\d|__|_t\d|_seg", tgt)[0]
    protected = {t for t in terms if t == tgt or t.startswith(stem + "_") or t.startswith(stem + "__") or t == stem}
    addA, addB = [], []
    for c in concepts(lab["missing_A"]):
        addA += [t for t in C.add_concept(df, c) if t not in terms and t not in addA]
    for c in concepts(lab["missing_B"]):
        addB += [t for t in C.add_concept(df, c) if t not in terms and t not in addB]
    conv = []
    for c in ["age", "sex", "race", "education", "pir", "smoking", "diabetes_dx", "hypertension_dx"]:
        conv += [t for t in C.add_concept(df, c) if t not in conv]
    # drugs
    pills = C.lib.raw(df, C.lib.stems("DIQ"), ["DIQ070", "DID070"])
    df["c_dm_rx"] = np.where((df.DIQ050 == 1) | (pills == 1), 1.0, np.where(df.DIQ010.notna(), 0.0, np.nan))
    df["c_bp_rx"] = np.where(df.BPQ_med_now == 1, 1.0, np.where(df.BPQ020.notna(), 0.0, np.nan))
    kset = set()
    for col in ("adjusts_component", "adjusts_mediator", "conditions_on_outcome", "conditions_on_collider"):
        kset |= set(concepts(lab[col]))
    tmap = {t: term_concept(t) for t in terms if t not in protected}
    # albumin/hemoglobin/liver map to the label concepts used by labels_i139 (component flags use index components)
    kterms = [t for t, c in tmap.items() if c in kset]
    allv = [y] + terms + addA + addB + conv + ["c_dm_rx", "c_bp_rx"]
    common = dom & df[allv].notna().all(axis=1) & (w > 0)
    tz, sd = C.standardise(df, tgt, common, w)
    tzs = lambda ts: [tz if t == tgt else t for t in ts]
    P = tzs(terms)
    binary = df.loc[common, tgt].nunique() <= 2
    models = {"P": P, "E3b": [t for t in P if t not in kterms], "conv": [tz] + conv, "P_dmrx": P + ["c_dm_rx"], "P_bprx": P + ["c_bp_rx"]}
    if addA:
        models["A"] = P + addA
    if addB:
        models["B"] = P + addB
    mm = common.to_numpy()
    ww = w[mm].to_numpy(float)
    yy = df.loc[mm, y].to_numpy(float)
    prep = {}
    res_log = []
    for k, ts in models.items():
        ts = [t for t in ts if t == tz or df.loc[common, t].nunique() > 1]
        Xm = np.column_stack([np.ones(mm.sum()), df.loc[mm, ts].to_numpy(float)])
        j = 1 + ts.index(tz)
        try:
            bl, _ = C.irls(Xm, yy, ww)
        except np.linalg.LinAlgError:
            bl = irls_ridge(Xm, yy, ww)   # singular design (collinear or separated terms): tiny ridge, logged
            res_log.append(f"model {k}: ridge fallback")
        bp = irls_pois(Xm, yy, ww)
        prep[k] = (Xm, j, bl, bp)
    point = {k: {"logor": float(v[2][v[1]]), "logpr": float(v[3][v[1]]),
                 "marg_logor": marginal_or(v[0], v[2], ww, v[1], binary)} for k, v in prep.items()}
    keys, Mx = C.raowu_multipliers(df[df.SDMVSTRA.notna()], B, f"20261114-i139-m3desc:{rid}")
    kpos = {k: i for i, k in enumerate(keys)}
    rk = df.SDMVSTRA.fillna(-1).to_numpy(np.int64)[mm] * 1000 + df.SDMVPSU.fillna(-1).to_numpy(np.int64)[mm]
    ki = np.array([kpos[v] for v in rk])
    boot_keys = [k for k in ("P", "A", "B", "E3b") if k in prep]
    bs = {k: {"logor": np.full(B, np.nan), "logpr": np.full(B, np.nan), "marg": np.full(B, np.nan)} for k in boot_keys}
    for bb in range(B):
        wb = ww * Mx[bb, ki]
        for k in boot_keys:
            Xm, j, bl, bp = prep[k]
            try:
                b1, _ = C.irls(Xm, yy, wb, beta0=bl, it=25, tol=1e-7)
                b2 = irls_pois(Xm, yy, wb, beta0=bp, it=25, tol=1e-7)
                bs[k]["logor"][bb] = b1[j]
                bs[k]["logpr"][bb] = b2[j]
                bs[k]["marg"][bb] = marginal_or(Xm, b1, wb, j, binary)
            except np.linalg.LinAlgError:
                pass

    def dlt(k1, k0, scale):
        pk = {"logor": "logor", "logpr": "logpr", "marg": "marg_logor"}[scale]
        d = point[k1][pk] - point[k0][pk]
        v = bs[k1][scale] - bs[k0][scale]
        v = v[np.isfinite(v)]
        lo, hi = np.percentile(v, 2.5), np.percentile(v, 97.5)
        return {"delta": d, "lo": float(lo), "hi": float(hi), "class": cls(d, lo, hi)}
    res = {"rec_id": rid, "status": "run", "n_common": int(mm.sum()), "weight_kind": kind, "binary_target": bool(binary),
           "term_map": tmap, "log": res_log, "k_concepts": sorted(kset), "k_terms_dropped_E3b": kterms, "point": point}
    for dag in ("A", "B"):
        if dag in prep:
            res[f"E3a_{dag}_PR"] = dlt(dag, "P", "logpr")
            res[f"E3a_{dag}_marginal"] = dlt(dag, "P", "marg")
    if kterms:
        res["E3b"] = {s: dlt("E3b", "P", s) for s in ("logor", "logpr", "marg")}
    for k in ("conv", "P_dmrx", "P_bprx"):
        res[f"shift_{k}"] = point[k]["logor"] - point["P"]["logor"]
    # RoE: subsets of covariate groups
    groups = {}
    for t in terms:
        if t in protected:
            continue
        g = t.split("__")[0]
        groups.setdefault(g, []).append(t)
    for t in addA + addB:
        groups.setdefault(t.split("__")[0], []).append(t)
    gl = sorted(groups)
    rng = np.random.default_rng(int(hashlib.sha256(f"20261114-i139-roe:{rid}".encode()).hexdigest()[:15], 16))
    if len(gl) <= 10:
        subsets = [s for r in range(len(gl) + 1) for s in itertools.combinations(gl, r)]
    else:
        subsets = [tuple(g for g in gl if rng.random() < 0.5) for _ in range(1024)]
    base = [tz] + [t for t in P if t in protected and t != tz]
    ors, n_bad = [], 0
    for s in subsets:
        ts = base + [t for g in s for t in groups[g]]
        ts = [t for t in ts if t == tz or df.loc[common, t].nunique() > 1]
        Xm = np.column_stack([np.ones(mm.sum()), df.loc[mm, ts].to_numpy(float)])
        bl = irls_ridge(Xm, yy, ww)
        v = bl[1 + ts.index(tz)] if bl is not None else np.nan
        if np.isfinite(v) and abs(v) < 10:
            ors.append(float(np.exp(v)))
        else:
            n_bad += 1
    ors = np.array(ors)
    res["RoE"] = {"n_models": int(len(ors)), "n_groups": len(gl), "p2_5": float(np.percentile(ors, 2.5)),
                  "p97_5": float(np.percentile(ors, 97.5)), "roe": float(np.percentile(ors, 97.5) / np.percentile(ors, 2.5)),
                  "share_or_gt1": float((ors > 1).mean()), "n_failed_fits": n_bad}
    res["seconds"] = int(time.time() - t0)
    return res


def main():
    T, L = C.spec_table(), C.labels()
    for rid in sys.argv[1:]:
        try:
            res = run(rid, T, L)
        except Exception as e:
            traceback.print_exc()
            res = {"rec_id": rid, "status": "error", "error": f"{type(e).__name__}: {e}"[:300]}
        json.dump(res, open(os.path.join(OUT, rid + ".json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1, default=float)
        print(rid, res["status"], "E3b", res.get("E3b", {}).get("logor", {}).get("class"), "PR_A", res.get("E3a_A_PR", {}).get("class"),
              "RoE", round(res.get("RoE", {}).get("roe", float("nan")), 3), res.get("seconds"), flush=True)


if __name__ == "__main__":
    main()
