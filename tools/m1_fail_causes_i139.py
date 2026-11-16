"""I139 M1: cause classification of every M1 failure (requested by blind review round 2, point 4(ii); the protocol
prereg_i139_v7_en.md does not define cause categories, so the categories are the reviewer's four (imputation, subsample,
undisclosed exclusion, unexplained) and the rule below is the most literal reading, written 15/11/2026 after M1 results
were seen; descriptive only, it does not change any M1 rank).

For each paper ranked "fail" in extraction/m1_results_i139.tsv, the builder (tools/m1_specs/<rec_id>.py, first
combination, unchanged) is loaded and the sample is decomposed:
  n_cc  = domain & y & all model terms non-missing & weight > 0 (reproduces n_repro of M1);
  n_xy  = domain & y & exposure terms non-missing & weight > 0 (before any covariate is required);
  per covariate: rows of n_xy lost because that covariate is missing, and whether it is a fasting-subsample variable
  (missing in >= 90 % of the n_xy rows whose pooled fasting weight is 0, and in < 50 % of rows with a positive one).
  The full text (fulltext_txt/<rec_id>.txt, before the reference list) is searched for an imputation statement
  (imput, MICE in capitals, chained equation, missForest, na.roughfix; "not imputed" / "without imputation" ignored);
  three papers where the sentence describes imputation only as a sensitivity analysis are set by hand (IMP_SENSITIVITY_ONLY).
Components: OR fails when |dlogOR| > log 1.25 (and the rounding rule did not rescue it); N fails when |dN|/N > 0.15.
Cause of an N failure, first rule that applies:
  1. n_repro > 1.15 x n_reported                      -> "undisclosed exclusion" (the paper removed people it does not
                                                          describe well enough to rebuild);
  2. n_xy < 0.85 x n_reported and the exposure or outcome needs the fasting or dietary subsample, or the weight
     used is a subsample weight that removes the shortfall (n_xy without the weight > 0 condition >= 0.85 x n_reported)
                                                       -> "subsample" (exposure, outcome or weight);
  3. n_xy < 0.85 x n_reported otherwise                -> "unexplained" (population or exposure definition);
  4. covariates cause the shortfall and the paper states imputation -> "imputation";
  5. covariates cause the shortfall and the covariate losing most rows is a fasting-subsample variable -> "subsample";
  6. otherwise                                          -> "unexplained" (covariate missingness, imputation not stated).
Secondary flags (imputation stated, fasting covariate) are recorded even when not the first rule.
An OR failure has no cause category in the reviewer's list; it is labelled "OR differs" with the N cause beside it.
Usage: python -I tools/m1_fail_causes_i139.py [rec_id ...]   -> extraction/m1_fail_causes/<rec_id>.json
       python -I tools/m1_fail_causes_i139.py --summary      -> extraction/m1_fail_causes_i139.tsv
"""
import os, sys, re, csv, json, glob, math, traceback
import numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m23_common_i139 as C

lib = C.lib
OUT = os.path.join(C.X, "m1_fail_causes")
os.makedirs(OUT, exist_ok=True)
FT = os.path.join(C.ROOT, "fulltext_txt")
IMP_RX = re.compile(r"(?i:imput|chained equation|missforest|na\.roughfix)|\bMICE\b")   # "mice" in lower case is the animal
NOT_IMP_RX = re.compile(r"not imputed|without imputation", re.I)
# Manual reading of the quoted sentence (15/11/2026): imputation used only as a sensitivity analysis, main analysis
# complete case or not stated, so it is not the cause under rule 4.
IMP_SENSITIVITY_ONLY = {"9835a67923": "MI is the second of 'two sensitivity analyses'",
                        "76d45e5c64": "MI then 'repeated the primary analyses'; Discussion: 'to enhance the robustness'",
                        "e1d9d038bb": "MI performed 'to evaluate missing-data impact'"}
REF_RX = re.compile(r"^\s*#*\s*(references|bibliography)\s*$", re.I)
LOG125 = math.log(1.25)


def imputation_lines(rid):
    p = os.path.join(FT, rid + ".txt")
    if not os.path.exists(p):
        return None
    out = []
    for i, line in enumerate(open(p, encoding="utf-8", errors="replace"), 1):
        if REF_RX.match(line):
            break
        m = IMP_RX.search(line)
        if m and not NOT_IMP_RX.search(line[max(0, m.start() - 60): m.end() + 60]):
            out.append([i, line[max(0, m.start() - 120): m.end() + 120].strip()])
    return out


def run(rid, m1row):
    b = C.load_builder(rid)
    df, y, terms, tgt = b["df"], b["y"], list(b["terms"]), b["target"]
    dom = b["domain"].fillna(False).astype(bool)
    if b["weight"] is None:
        w = np.ones(len(df))
    else:
        w = df[b["weight"]].fillna(0).to_numpy(float)
    wpos = w > 0
    stem = re.split(r"_q\d|__|_t\d|_seg", tgt)[0]
    expo = [t for t in terms if t == tgt or t == stem or t.startswith(stem + "_")]
    cov = [t for t in terms if t not in expo]
    base = dom.to_numpy() & df[[y]].notna().all(axis=1).to_numpy() & wpos
    xy = base & df[expo].notna().all(axis=1).to_numpy()
    cc = xy & df[cov].notna().all(axis=1).to_numpy() if cov else xy
    xy_now = dom.to_numpy() & df[[y] + expo].notna().all(axis=1).to_numpy()   # before the weight > 0 condition
    wf = lib.pooled_weight(df, "fasting").fillna(0).to_numpy(float)
    f0, f1 = xy & (wf <= 0), xy & (wf > 0)

    def fasting_like(col):
        miss = df[col].isna().to_numpy()
        return bool(f0.sum() and f1.sum() and miss[f0].mean() >= 0.90 and miss[f1].mean() < 0.50)
    loss = []
    for t in cov:
        miss = df[t].isna().to_numpy()
        loss.append([t, int((xy & miss).sum()), fasting_like(t)])
    loss.sort(key=lambda r: -r[1])
    expo_fast = any(fasting_like(t) for t in expo) or fasting_like(y)
    diet = rid in lib.DIET
    nr, npp = float(m1row["n_reported"]), int(m1row["n_repro"])
    return {"rec_id": rid, "n_reported": nr, "n_repro": npp, "n_cc": int(cc.sum()), "n_xy": int(xy.sum()),
            "n_domain_y": int(base.sum()), "n_xy_noweight": int(xy_now.sum()), "weight": b["weight"], "exposure_terms": expo,
            "exposure_or_outcome_fasting": expo_fast, "dietary_index": diet,
            "covariate_loss_top5": loss[:5], "imputation_lines": imputation_lines(rid)}


def classify(d, m1row):
    dl, dn = float(m1row["dlogor"]), float(m1row["dn_rel"])
    nr = d["n_reported"]
    or_fail, n_fail = dl > LOG125, dn > 0.15
    lines = imputation_lines(d["rec_id"])
    imp_any = bool(lines)
    imp = imp_any and d["rec_id"] not in IMP_SENSITIVITY_ONLY
    top = d["covariate_loss_top5"][0] if d["covariate_loss_top5"] else None
    fast_cov = any(r[2] for r in d["covariate_loss_top5"] if r[1] > 0)
    cause = ""
    if n_fail:
        if d["n_repro"] > 1.15 * nr:
            cause = "undisclosed exclusion"
        elif d["n_xy"] < 0.85 * nr:
            sub_w = d.get("n_xy_noweight", 0) >= 0.85 * nr and d["n_xy"] < 0.85 * d.get("n_xy_noweight", 0)
            cause = "subsample" if (d["exposure_or_outcome_fasting"] or d["dietary_index"] or sub_w) else "unexplained"
        elif imp:
            cause = "imputation"
        elif top and top[2]:
            cause = "subsample"
        else:
            cause = "unexplained"
    comp = "OR and N" if or_fail and n_fail else ("OR" if or_fail else ("N" if n_fail else "neither (check)"))
    return {"rec_id": d["rec_id"], "failing_component": comp, "n_cause": cause, "or_differs": or_fail,
            "imputation_stated": imp_any, "imputation_main_analysis": imp, "fasting_covariate_in_top5": fast_cov,
            "n_reported": int(nr), "n_xy_noweight": d.get("n_xy_noweight"), "n_xy": d["n_xy"], "n_repro": d["n_repro"], "n_cc_check": d["n_cc"],
            "dlogor": dl, "dn_rel": dn, "top_covariate_loss": f"{top[0]}:{top[1]}" if top else ""}


def main():
    m1 = {r["rec_id"]: r for r in csv.DictReader(open(os.path.join(C.X, "m1_results_i139.tsv"), encoding="utf-8"), delimiter="\t")}
    if sys.argv[1:] == ["--summary"]:
        rows = []
        for f in sorted(glob.glob(os.path.join(OUT, "*.json"))):
            d = json.load(open(f, encoding="utf-8"))
            if d.get("status") == "error":
                rows.append({"rec_id": d["rec_id"], "failing_component": "error"})
                continue
            rows.append(classify(d, m1[d["rec_id"]]))
        cols = list(rows[0].keys())
        with open(os.path.join(C.X, "m1_fail_causes_i139.tsv"), "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            w.writeheader(); w.writerows(rows)
        for r in rows:
            print("\t".join(str(r.get(c, "")) for c in cols))
        return
    ids = sys.argv[1:] or [r for r, v in m1.items() if v["rank"] == "fail"]
    for rid in ids:
        try:
            d = run(rid, m1[rid])
        except Exception as e:
            traceback.print_exc()
            d = {"rec_id": rid, "status": "error", "error": f"{type(e).__name__}: {e}"[:300]}
        json.dump(d, open(os.path.join(OUT, rid + ".json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1, default=str)
        print(rid, d.get("n_reported"), d.get("n_xy"), d.get("n_repro"), d.get("n_cc"), bool(d.get("imputation_lines")), flush=True)


if __name__ == "__main__":
    main()
