"""I139 round-3 review points M1 and M2 (paper6 main v3, section 7.3). Written 16/11/2026 after all main results were
seen; post hoc and descriptive, no model refit, no new extraction.

M1: N_h*-weighted proportions next to the unweighted ones for E1 (strict pass; near pass under denominators i, ii,
    iii), E2 and every reporting-checklist item. Weights and estimator are those of the registered E3a analysis
    (m3_summary_i139.NHS, m3_summary_i139.stratified: fpc, logit CI, stratified bootstrap 2000 when a stratum is 0 or 1).
    Denominator (ii) counts not-testable papers as fail, (iii) as pass. Checklist items exclude "not_applicable".
M2: the same quantities on the 33 primary papers alone (the 14 replacement papers removed), and extreme bounds for
    E1(i) over the 13 paywalled slots (R13 as in sens_sample_i139.py): the 13 replacements are removed and the 13
    unseen originals are counted as all fail (lower) or all pass (upper), each with a Wilson 95 % CI.
Usage: python -I tools/r3_weighted_primary_i139.py   -> extraction/r3_weighted_primary_i139.json
"""
import os, sys, csv, json
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_summary_i139 as S

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
X = os.path.join(ROOT, "extraction")
rd = lambda p: list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"))
PASS = ("strict_pass", "near_pass")


def both(items):
    """items: list of (stratum, 0/1) -> unweighted Wilson and N_h*-weighted estimate."""
    k = sum(v for _, v in items)
    return {"unweighted": S.wilson(k, len(items)), "weighted": S.stratified(items) if items else None}


def block(ids, fam, m1, m2, chk):
    run = [i for i in ids if m1[i]["status"] == "run"]
    nt = [i for i in ids if m1[i]["status"] != "run"]
    out = {"n_papers": len(ids), "n_testable": len(run)}
    out["E1_strict_i"] = both([(fam[i], int(m1[i]["rank"] == "strict_pass")) for i in run])
    out["E1_i_testable"] = both([(fam[i], int(m1[i]["rank"] in PASS)) for i in run])
    out["E1_ii_nt_fail"] = both([(fam[i], int(i in run and m1[i]["rank"] in PASS)) for i in ids])
    out["E1_iii_nt_pass"] = both([(fam[i], int(i in nt or m1[i]["rank"] in PASS)) for i in ids])
    e2 = [i for i in ids if i in m2 and m2[i]["E2_lt_080"] in ("True", "False")]
    out["E2"] = both([(fam[i], int(m2[i]["E2_lt_080"] == "True")) for i in e2])
    items = sorted({int(r["item"]) for r in chk})
    by = {(r["rec_id"], int(r["item"])): r["rating"] for r in chk}
    out["checklist"] = {}
    for it in items:
        rated = [i for i in ids if by.get((i, it)) in ("yes", "no")]
        na = sum(by.get((i, it)) == "not_applicable" for i in ids)
        out["checklist"][it] = dict(both([(fam[i], int(by[(i, it)] == "yes")) for i in rated]), n_not_applicable=na)
    return out


def main():
    sample = rd(os.path.join(X, "sample_47.tsv"))
    fam = {r["rec_id"]: r["family"] for r in sample}
    ids = [r["rec_id"] for r in sample]
    m1 = {r["rec_id"]: r for r in rd(os.path.join(X, "m1_results_i139.tsv"))}
    m2 = {r["rec_id"]: r for r in rd(os.path.join(X, "m2_results_i139.tsv"))}
    chk = [r for r in rd(os.path.join(X, "nchs_checklist_i139.tsv")) if r["rec_id"] in fam]
    acc = json.load(open(os.path.join(X, "sens_sample_i139.json"), encoding="utf-8"))["replacement_accounting"]
    r13, r14 = acc["r13"], acc["r14"]
    assert len(r13) == 13 and len(r14) == 14 and len(ids) == 47
    prim = [i for i in ids if i not in r14]
    assert len(prim) == 33
    out = {"N_h_star": {h: S.tr(v) for h, v in S.NHS.items()},
           "M1_all_47": block(ids, fam, m1, m2, chk),
           "M2_primary_33": block(prim, fam, m1, m2, chk)}
    keep = [i for i in ids if i not in r13]
    run = [i for i in keep if m1[i]["status"] == "run"]
    k = sum(m1[i]["rank"] in PASS for i in run)
    out["M2_E1_i_bounds_R13"] = {"testable_without_R13": len(run), "pass_without_R13": k,
                                 "lower_all_13_fail": S.wilson(k, len(run) + 13),
                                 "upper_all_13_pass": S.wilson(k + 13, len(run) + 13),
                                 "observed_full_sample": S.wilson(sum(m1[i]["rank"] in PASS for i in ids if m1[i]["status"] == "run"),
                                                                  sum(m1[i]["status"] == "run" for i in ids))}
    json.dump(out, open(os.path.join(X, "r3_weighted_primary_i139.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
