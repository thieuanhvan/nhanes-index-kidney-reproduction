"""I139 M1 sensitivity "any combination passes" (prereg_i139_v7_en.md Section 6, M1: "Outcomes use the first
combination. Sensitivity: 'any combination passes', with the number of passing combinations over the number of
combinations run"). Written 15/11/2026 after the M1 results were seen; descriptive only, the M1 ranks and E1 are not changed.

Combinations come from the reading_1 / reading_2 pairs of the locked specification table (extraction/spec_table_i139.json,
field "ambiguities"). The first combination is the locked builder (tools/m1_specs/<rec_id>.py) unchanged; every other
combination applies the reading_2 of one or more points to that builder through the operations below, then fits and
ranks exactly as m1_run_i139 (same svy_logit, same rank() with the rounding rule). Points are implemented only where
reading_2 is fully defined with variables the builder already uses or with the fixed concept codings of
m23_common_i139.add_concept; the others are listed in NOT_RUN with the reason, so "combinations run" can be smaller than
"combinations defined". Points that do not touch the M1 estimate (inflection point, effect label, an eGFR or component
reading that gives the same computation as reading_1) are not combinations. At most 8 combinations per paper (prereg).
Operations: uacr up/down (people with UACR exactly 30 moved just above / below 30, i.e. >= 30 versus > 30);
egfr_swap (CKD-EPI 2009 and 2021 columns swapped); weight (mec, interview, diet day 1, or the M2/M3 rule-1 correct
weight); add / drop terms; target (another term of the same model); cont (the continuous index instead of its
categories; none, ln or per SD on the domain, unweighted SD as lib.sd_scale); quant (k quantile groups on the domain,
level j versus 1); y_both (CKD or DKD requiring both eGFR < 60 and albuminuria); reported (another reported N).
Usage: python -I tools/m1_anycomb_i139.py <rec_id> [...]   -> extraction/m1_anycomb/<rec_id>.json
       python -I tools/m1_anycomb_i139.py --summary        -> extraction/m1_anycomb_summary_i139.json
"""
import os, sys, json, glob, csv, itertools, traceback, math
import numpy as np
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m23_common_i139 as C
from m1_run_i139 import rank

lib = C.lib
OUT = os.path.join(C.X, "m1_anycomb")
os.makedirs(OUT, exist_ok=True)
STATE = {"uacr": None, "egfr_swap": False}
_orig_table = lib.table


def _kidney(df):
    return df


def table(cycles):
    df = _orig_table(cycles)
    if STATE["egfr_swap"]:
        a, b = df.egfr_ckdepi2009.copy(), df.egfr_ckdepi2021.copy()
        df["egfr_ckdepi2009"], df["egfr_ckdepi2021"] = b, a
    if STATE["uacr"]:
        tie = (df.uacr - 30).abs() < 1e-6
        df.loc[tie, "uacr"] = 30.001 if STATE["uacr"] == "up" else 29.999
    return df


lib.table = table

Q4 = ["x_q2", "x_q3", "x_q4"]
# reading_2 alternatives per point; each point is a list of alternatives, alternative 0 is reading_1 (no operation)
POINTS = {
    "f6e85f9043": [("covariates: Table 3 list adds hypertension", [[], [("add", ["concept:hypertension_dx"])]])],
    "f1ad3fb0ff": [("UACR >= 30", [[], [("uacr", "up")]])],
    "d22fcf0a95": [("estimate: continuous NPAR 1.099 or Q4 vs Q1 1.86", [[], [("cont", "x", "none", Q4, 1.099)], [("target", "x_q4", 1.86)]])],
    "94132d0c1a": [("ACR > 30", [[], [("uacr", "down")]])],
    "fad26f18b1": [("age and BMI continuous", [[], [("drop", ["age2__*", "bmi3__*"]), ("add", ["age", "bmi"])]])],
    "cef539b61e": [("weight: pooled MEC", [[], [("weight", "mec")]])],
    "9cd1df095e": [("covariates: Methods list without serum creatinine", [[], [("drop", ["scr"])]]),
                   ("weight: WTMECPRP", [[], [("weight", "mec")]])],
    "53ca68a923": [("estimate: ln CALLY continuous 0.878", [[], [("cont", "cally", "ln", ["cally_q2", "cally_q3", "cally_q4"], 0.878)]])],
    "c62af9091f": [("eGFR CKD-EPI 2021", [[], [("egfr_swap",)]])],
    "cf0849f17f": [("covariates: Methods list without BMI, Scr, eGFR", [[], [("drop", ["bmi_g__*", "scr", "egfr_ckdepi2009"])]])],
    "4b72270032": [("estimate: Q4 vs Q1 1.83", [[], [("quant", "x", 4, 4, ["x"], 1.83)]])],
    "8991f10763": [("covariates: add blood uric acid", [[], [("add", ["uric_acid"])]])],
    "71eedfeafb": [("eGFR CKD-EPI 2021", [[], [("egfr_swap",)]])],
    "988df86a79": [("reported N 19,371", [[], [("reported", 19371)]])],
    "30eae6044c": [("weighted (correct weight)", [[], [("weight", "correct")]])],
    "aca6f864fb": [("ACR >= 30", [[], [("uacr", "up")]]), ("weight: dietary day 1", [[], [("weight", "diet")]])],
    "78c5287646": [("reported N 14,349 (whole sample)", [[], [("reported", 14349)]])],
    "9835a67923": [("estimate: continuous DOBS 0.98", [[], [("cont", "dobs", "none", ["dobs_q2", "dobs_q3", "dobs_q4"], 0.98)]])],
    "b966d3cfdc": [("estimate: continuous OBS 0.96", [[], [("cont", "x", "none", Q4, 0.96)]]),
                   ("CKD requires both criteria", [[], [("y_both", "ckd_ge30")]])],
    "1c8c318fda": [("covariates: only the nine added variables", [[], [("drop", ["female", "age", "race1__*", "educ3__*", "married", "bmi"])]])],
    "cde66dda6a": [("CKD requires both criteria", [[], [("y_both", "ckd_gt30")]]), ("weight: pooled MEC", [[], [("weight", "mec")]])],
    "ea03123ab4": [("weighted regression (correct weight)", [[], [("weight", "correct")]])],
    "76d45e5c64": [("estimate: continuous UHR 1.10", [[], [("cont", "x", "none", Q4, 1.10)]])],
    "e76d3d1271": [("estimate: T3 vs T1 1.22", [[], [("quant", "x", 3, 3, ["x"], 1.22)]])],
    "450ed6ebb6": [("estimate: RC per SD 1.36", [[], [("cont", "x", "sd", Q4, 1.36)]]),
                   ("DKD requires ACR >= 30 and eGFR < 60", [[], [("y_both", "dkd_ge30")]])],
    "69700c2b09": [("weight: pooled MEC", [[], [("weight", "mec")]])],
    "250e5761d5": [("covariates: Methods list adds PIR", [[], [("add", ["pir"])]]), ("weighted (correct weight)", [[], [("weight", "correct")]]),
                   ("ACR >= 30", [[], [("uacr", "up")]])],
    "ff27565860": [("estimate: Q4 vs Q1 1.22", [[], [("quant", "x", 4, 4, ["x"], 1.22)]]), ("covariates: add HDL-C and TC", [[], [("add", ["hdl", "tc"])]]),
                   ("weighted (correct weight)", [[], [("weight", "correct")]])],
    "259b9d6255": [("estimate: continuous WWI 1.21", [[], [("cont", "x", "none", ["x_q2", "x_q3", "x_q4"], 1.21)]]), ("ACR > 30", [[], [("uacr", "down")]])],
    "ac96396902": [("covariates: add diabetes", [[], [("add", ["concept:diabetes_dx"])]])],
    "c48bec549c": [("ACR >= 30", [[], [("uacr", "up")]])],
    "ba21203b93": [("estimate: continuous FPSA/TPSA x 100, 1.03", [[], [("cont", "x", "none", ["x_t1", "x_t3"], 1.03)]])],
    "4de61a08ed": [("estimate: upper segment 1.24", [[], [("target", "seg_hi", 1.24)]]), ("covariates: add LDL-C, HDL-C, FBG, HbA1c", [[], [("add", ["ldl", "hdl", "glucose", "hba1c"])]])],
    "d84bd11e2e": [("covariates: add BMI", [[], [("add", ["bmi"])]])],
    "a906145f51": [("estimate: Q5 vs Q3 1.79", [[], [("target", "x_q5v3", 1.79)]]), ("weight: interview", [[], [("weight", "interview")]])],
}
NOT_RUN = {
    "79f1e21c4c": ["cycles 1999-2018 (J file): builder coded for the P file (ALQ121, DMDMARTZ branches); needs a new builder"],
    "f1ad3fb0ff": ["model with metabolic syndrome or log HOMA-IR: terms not in the builder", "Table 6 mutually adjusted model: terms not in the builder"],
    "6b7f67827b": ["hospital cohort estimate: not NHANES"],
    "94132d0c1a": ["cycles 1999-2018 (J file): builder coded for the P file"],
    "ea130293e0": ["hospital cohort estimate: not NHANES", "Covariates-section list: needs a new builder"],
    "fad26f18b1": ["KDIGO chronicity: not possible with one NHANES visit"],
    "4b72270032": ["abstract covariate description: list not specified"],
    "78c5287646": ["footnote list with lymphocyte and PNI group: PNI group is a coding of the exposure"],
    "8a9e8d47ec": ["CMI as exposure: units of the CMI estimate not stated"],
    "b966d3cfdc": ["five of six cycles: the omitted cycle is not named"],
    "1c8c318fda": ["CHARLS estimate: not NHANES"],
    "e1d9d038bb": ["Model 4: covariate list not in the builder", "continuous per-SD ln-UHR exists only at Model 4"],
    "450ed6ebb6": ["RCII as exposure: needs a new builder"],
    "ac96396902": ["CVD-CKD comorbidity outcome: excluded by the pre-written rule 4"],
    "c48bec549c": ["exposure scaling not stated (reading_2 undefined)"],
    "ba21203b93": ["TPSA as exposure: needs a new builder", "abstract covariate list: needs a new builder"],
    "3865be04bf": ["cycles 2005-2018 (J file): builder coded for the P file", "marital status: coding not stated"],
    "edaf91172c": ["scaling not stated (reading_2 undefined)", "unspecified 'other key variables'"],
    "aca6f864fb": ["Model 2: gives the same reported OR; its covariate list is not in the builder"],
    "d84bd11e2e": ["type 1 diabetes not separated: no NHANES variable"],
    "a906145f51": ["RCS-segment estimates: model not labelled"],
}


def expand(b, names):
    out = []
    for n in names:
        if n.endswith("*"):
            out += [t for t in b["terms"] if t.startswith(n[:-1])]
        else:
            out.append(n)
    return out


def build(rid, ops):
    STATE["uacr"] = next((o[1] for o in ops if o[0] == "uacr"), None)
    STATE["egfr_swap"] = any(o[0] == "egfr_swap" for o in ops)
    b = C.load_builder(rid)
    df = b["df"]
    dom = b["domain"].fillna(False).astype(bool)
    rep_est, rep_n = None, None
    for o in ops:
        k = o[0]
        if k == "add":
            for t in o[1]:
                new = C.add_concept(df, t.split(":")[1]) if t.startswith("concept:") else [t]
                b["terms"] += [x for x in new if x not in b["terms"]]
        elif k == "drop":
            d = set(expand(b, o[1]))
            missing = [t for t in o[1] if not t.endswith("*") and t not in b["terms"]]
            if missing:
                raise KeyError(f"drop: terms not in model {missing}")
            b["terms"] = [t for t in b["terms"] if t not in d]
        elif k == "weight":
            if o[1] == "correct":
                w, _ = C.correct_weight(rid, b)
            elif o[1] == "diet":
                w = C.diet_weight(df)
            else:
                w = lib.pooled_weight(df, o[1])
            df["_wv"] = w
            b["weight"] = "_wv"
        elif k == "target":
            if o[1] not in b["terms"]:
                raise KeyError(f"target {o[1]} not in model")
            b["target"], rep_est = o[1], o[2]
        elif k == "cont":
            col, trf, dropt, rep_est = o[1], o[2], o[3], o[4]
            v = df[col].astype(float)
            if trf == "ln":
                v = np.log(v.where(v > 0))
            df["_xc"] = v
            if trf == "sd":
                df["_xc"] = lib.sd_scale(df, "_xc", dom)
            b["terms"] = ["_xc"] + [t for t in b["terms"] if t not in dropt]
            b["target"] = "_xc"
        elif k == "quant":
            col, kq, lev, dropt, rep_est = o[1], o[2], o[3], o[4], o[5]
            names, _ = lib.quantile_dummies(df, col, dom & df[col].notna(), k=kq, prefix="_xq")
            b["terms"] = names + [t for t in b["terms"] if t not in dropt]
            b["target"] = f"_xq{lev}"
        elif k == "y_both":
            e, u = df.egfr_ckdepi2009, df.uacr
            alb = (u >= 30) if o[1] in ("ckd_ge30", "dkd_ge30") else (u > 30)
            yy = ((e < 60) & alb)
            ok = e.notna() & u.notna()
            if o[1] == "dkd_ge30":
                yy = yy & (df.dm == 1)
                ok = ok & df.dm.notna()
            df["_yb"] = yy.astype(float).where(ok)
            b["y"] = "_yb"
        elif k == "reported":
            rep_n = o[1]
    return b, rep_est, rep_n


def run(rid, spec):
    pts = POINTS.get(rid, [])
    combos = list(itertools.product(*[alts for _, alts in pts]))[:8] if pts else [()]
    rows = []
    for ci, combo in enumerate(combos):
        ops = [o for alt in combo for o in alt]
        label = "; ".join(pts[i][0] for i, alt in enumerate(combo) if alt) or "first combination"
        try:
            b, rep_est, rep_n = build(rid, ops)
            if b.get("not_testable"):
                rows.append({"combo": ci, "label": label, "status": "not_testable"})
                continue
            df = b["df"]
            w = b["weight"]
            if w is None:
                df["_w1"] = 1.0
                w = "_w1"
            res = lib.svy_logit(df, b["y"], b["terms"], w, domain=b["domain"])
            t = res["table"].loc[b["target"]]
            est = float(rep_est if rep_est is not None else spec["estimate"])
            nr = float(rep_n if rep_n is not None else spec["n_analytic"])
            st, d, dn = rank(est, float(t.OR), nr, res["n"])
            rows.append({"combo": ci, "label": label, "status": "run", "rank": st, "or_reported": est, "or_repro": float("%.4g" % t.OR),
                         "dlogor": math.floor(d * 1e4) / 1e4, "n_reported": nr, "n_repro": int(res["n"]), "dn_rel": math.floor(dn * 1e4) / 1e4})
        except Exception as e:
            traceback.print_exc()
            rows.append({"combo": ci, "label": label, "status": "error", "error": f"{type(e).__name__}: {e}"[:300]})
    n_def = 2 ** (len(pts) + len(NOT_RUN.get(rid, [])))
    return {"rec_id": rid, "points_run": [p for p, _ in pts], "points_not_run": NOT_RUN.get(rid, []),
            "combinations_defined": min(8, n_def), "combinations": rows}


def summary():
    m1 = {r["rec_id"]: r for r in csv.DictReader(open(os.path.join(C.X, "m1_results_i139.tsv"), encoding="utf-8"), delimiter="\t")}
    R = {d["rec_id"]: d for d in (json.load(open(f, encoding="utf-8")) for f in glob.glob(os.path.join(OUT, "*.json")))}
    testable = [r for r, v in m1.items() if v["status"] == "run"]
    out = {"papers_testable": len(testable), "per_paper": {}}
    any_pass = any_strict = n_run = n_pass = n_def = n_err = 0
    changed = []
    for r in testable:
        d = R.get(r)
        rows = d["combinations"] if d else []
        ran = [x for x in rows if x["status"] == "run"]
        p = [x for x in ran if x["rank"] in ("strict_pass", "near_pass")]
        first = m1[r]["rank"]
        first_rep = rows[0].get("rank") if rows else None
        ap = bool(p) or first in ("strict_pass", "near_pass")
        any_pass += ap
        any_strict += any(x["rank"] == "strict_pass" for x in ran) or first == "strict_pass"
        n_run += len(ran); n_pass += len(p); n_err += sum(x["status"] == "error" for x in rows)
        n_def += d["combinations_defined"] if d else 1
        if ap and first == "fail":
            changed.append({r: [x["label"] for x in p]})
        out["per_paper"][r] = {"first_combination_rank": first, "first_reproduced_here": first_rep, "run": len(ran), "pass": len(p),
                               "defined": d["combinations_defined"] if d else 1, "not_run_points": d["points_not_run"] if d else []}
    out.update(any_pass=C_wilson(any_pass, len(testable)), any_strict=any_strict, combinations_run=n_run, combinations_passing=n_pass,
               combinations_defined=n_def, combinations_error=n_err, fail_to_pass=changed,
               first_combination_mismatch=[r for r, v in out["per_paper"].items() if v["first_reproduced_here"] not in (None, v["first_combination_rank"])])
    json.dump(out, open(os.path.join(C.X, "m1_anycomb_summary_i139.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "per_paper"}, ensure_ascii=False, indent=1))


def C_wilson(k, n):
    import m3_summary_i139 as S
    return S.wilson(k, n)


def main():
    if sys.argv[1:] == ["--summary"]:
        return summary()
    T = C.spec_table()
    for rid in sys.argv[1:]:
        res = run(rid, T[rid])
        json.dump(res, open(os.path.join(OUT, rid + ".json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1, default=float)
        print(rid, [(x.get("label"), x.get("rank") or x.get("status")) for x in res["combinations"]], flush=True)


if __name__ == "__main__":
    main()
