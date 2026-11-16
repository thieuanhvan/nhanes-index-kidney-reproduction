"""I139 M3 (confounder specification, main analysis) per extraction/m2_m3_rules_i139.md.
Usage: python -I tools/m3_run_i139.py <rec_id> [...]   (writes extraction/m3/<rec_id>.json, one file per paper)
Aggregation: tools/m3_summary_i139.py
"""
import os, sys, json, hashlib, time, traceback
import numpy as np, pandas as pd
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m23_common_i139 as C

B = 500
OUT = os.path.join(C.X, os.environ.get("I139_M3_OUT", "m3"))
os.makedirs(OUT, exist_ok=True)


def cls(d, lo, hi):
    if d is None or np.isnan(d):
        return "na"
    if lo > C.LOG110 or hi < -C.LOG110:
        return "firm_positive"
    if abs(d) > C.LOG110:
        return "point_positive"
    return "negative"


def concepts(s):
    return [c for c in (s or "").split("|") if c]


def run(rid, T, L):
    t0 = time.time()
    spec, lab = T[rid], L[rid]
    res = {"rec_id": rid, "family": spec["family"], "dag_family": lab["dag_family"], "flag_outcome_component": C.outcome_flag(spec),
           "missing_A": lab["missing_A"], "missing_B": lab["missing_B"], "label_A": lab["label_A"], "label_B": lab["label_B"]}
    b = C.load_builder(rid)
    if b.get("not_testable"):
        res["status"] = "not_testable"
        res["reason"] = str(b["not_testable"])[:300]
        return res
    df, y = b["df"], b["y"]
    dom = b["domain"].fillna(False).astype(bool)
    terms = list(b["terms"])
    tgt = b["target"]
    w, kind = C.correct_weight(rid, b)
    df["_wc"] = w
    log = []

    def addset(cs):
        out = []
        for c in cs:
            for t in C.add_concept(df, c):
                if t in terms or t in out:
                    log.append(f"{c}: term {t} already in model, not added twice")
                else:
                    out.append(t)
        return out
    addA = addset(concepts(lab["missing_A"]))
    addB = addset(concepts(lab["missing_B"]))
    allv = [y] + terms + addA + addB
    common = dom & df[allv].notna().all(axis=1) & (w > 0)
    psample = dom & df[[y] + terms].notna().all(axis=1) & (w > 0)
    tz, sd = C.standardise(df, tgt, common, w)
    tz_terms = [tz if t == tgt else t for t in terms]
    # negative control
    rng = np.random.default_rng(int(hashlib.sha256(f"20261114-i139-negctl:{rid}".encode()).hexdigest()[:15], 16))
    pir_miss = float(df.loc[dom, "pir"].isna().mean())
    r = rng.standard_normal(len(df))
    r[rng.random(len(df)) < pir_miss] = np.nan
    df["_rand"] = r
    negm = common & df["_rand"].notna()
    at = C.age_terms(tz_terms)
    models = {"P": (common, tz_terms), "Ps": (psample, tz_terms), "Pn": (negm, tz_terms), "N": (negm, tz_terms + ["_rand"])}
    if addA:
        models["A"] = (common, tz_terms + addA)
    if addB and set(addB) != set(addA):
        models["B"] = (common, tz_terms + addB)
    if at:
        models["Q"] = (common, [t for t in tz_terms if t not in at])
    # drop terms constant in the model's sample (singular design), logged
    for k, (m, ts) in list(models.items()):
        keep = [t for t in ts if t == tz or df.loc[m, t].nunique() > 1]
        if len(keep) < len(ts):
            log.append(f"model {k}: constant terms dropped {sorted(set(ts) - set(keep))}")
        models[k] = (m, keep)
    # design-correct full-sample fits (Taylor CI, for reporting)
    full = {}
    for k in models:
        m, ts = models[k]
        fit = C.nc.svy_logit(df, y, ts, "_wc", domain=m)
        row = fit["table"].loc[tz]
        full[k] = {"logor": float(row.coef), "se": float(row.se), "n": fit["n"]}
    # Rao-Wu paired bootstrap
    keys, Mx = C.raowu_multipliers(df[df.SDMVSTRA.notna()], B, f"20261114-i139-m3:{rid}")
    kpos = {k: i for i, k in enumerate(keys)}
    rowkey = (df.SDMVSTRA.fillna(-1).to_numpy(np.int64) * 1000 + df.SDMVPSU.fillna(-1).to_numpy(np.int64))
    prep = {}
    for k, (m, ts) in models.items():
        mm = m.to_numpy()
        Xm = np.column_stack([np.ones(mm.sum()), df.loc[mm, ts].to_numpy(float)])
        yy = df.loc[mm, y].to_numpy(float)
        ww = w[mm].to_numpy(float)
        ki = np.array([kpos[v] for v in rowkey[mm]])
        b0, _ = C.irls(Xm, yy, ww)
        prep[k] = (Xm, yy, ww, ki, b0, 1 + ts.index(tz))
    est = {k: np.full(B, np.nan) for k in models}
    fails = 0
    for bb in range(B):
        for k, (Xm, yy, ww, ki, b0, j) in prep.items():
            try:
                beta, _ = C.irls(Xm, yy, ww * Mx[bb, ki], beta0=b0, it=25, tol=1e-7)
                est[k][bb] = beta[j] if np.all(np.isfinite(beta)) else np.nan
            except np.linalg.LinAlgError:
                fails += 1
    res.update(status="run", weight_kind=kind, target=tgt, target_sd=sd, n_terms=len(terms), added_A=addA, added_B=addB,
               n_common=int(common.sum()), n_paper_sample=int(psample.sum()), n_neg=int(negm.sum()), pir_missing_rate=pir_miss,
               full=full, boot_failed_fits=fails, log=log)

    def delta(k1, k0):
        d = full[k1]["logor"] - full[k0]["logor"]
        v = est[k1] - est[k0]
        v = v[np.isfinite(v)]
        lo, hi = (np.percentile(v, 2.5), np.percentile(v, 97.5)) if len(v) > 20 else (np.nan, np.nan)
        return {"delta": d, "lo": float(lo), "hi": float(hi), "n_boot": int(len(v)), "class": cls(d, lo, hi)}
    res["sample_change"] = full["Ps"]["logor"] - full["P"]["logor"]
    res["E3a_A"] = delta("A", "P") if "A" in models else {"class": "negative", "note": "nothing missing under DAG A (or all present)"}
    if "B" in models:
        res["E3a_B"] = delta("B", "P")
    elif addB:
        res["E3a_B"] = dict(res["E3a_A"], note="DAG B adds the same variables as DAG A")
    else:
        res["E3a_B"] = {"class": "negative", "note": "nothing missing under DAG B"}
    res["neg_control"] = delta("N", "Pn")
    res["pos_control"] = delta("P", "Q") if "Q" in models else {"class": "na", "note": "no age term"}
    res["seconds"] = int(time.time() - t0)
    return res


def main():
    T, L = C.spec_table(), C.labels()
    for rid in sys.argv[1:]:
        try:
            res = run(rid, T, L)
        except Exception as e:
            traceback.print_exc()
            res = {"rec_id": rid, "status": "error", "error": f"{type(e).__name__}: {e}"[:400]}
        json.dump(res, open(os.path.join(OUT, rid + ".json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1, default=float)
        a, bb = res.get("E3a_A", {}), res.get("E3a_B", {})
        print(rid, res["status"], res.get("weight_kind"), "nC", res.get("n_common"), "A", a.get("class"), round(a.get("delta", np.nan), 4) if "delta" in a else "",
              "B", bb.get("class"), "neg", res.get("neg_control", {}).get("class"), "pos", res.get("pos_control", {}).get("class"),
              res.get("seconds"), "s", flush=True)


if __name__ == "__main__":
    main()
