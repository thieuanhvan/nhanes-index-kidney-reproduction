"""I139 sample-composition sensitivity analyses on existing outputs (no model refit). Written 15/11/2026 after all main
results were seen; descriptive, not used for the main conclusion.

1. Pilot exclusion (prereg_i139_v7_en.md Section 0: "every outcome is recomputed with the pilot papers excluded").
   Only B1, B3, B7 match the frame query (verdict_vi.md 9.10); B1 = f1ad3fb0ff and B3 = 79f1e21c4c are in the sample,
   B7 = 647b314a7c was a paywalled primary paper and is not. Recomputed: E1 (three denominators), E2, E3a(A) and E3a(B)
   (firm and point estimate). For the weighted E3a, the pilot papers are also removed from the F1 frame count
   (N_h F1 7 -> 4) and from the F1 full-text assessed and eligible counts (4/4 -> 2/2).
2. N_h sensitivity (prereg Section 3.2 step 8, plus verdict 9.9 "N_h include only"): weighted E3a(A) and E3a(B) with
   N_h from (a) reconciled include + unclear (main), (b) reconciled include only, (c) screening pass A alone (include +
   unclear), (d) pass B alone (include + unclear). Full-text eligibility proportions are kept. E1 and E2 are unweighted
   proportions and do not depend on N_h.
3. Replacement-paper bounds (reviewer round 2, point 6): E3a(A) and E3a(B) primary, with replacement papers treated
   (a) as not testable (dropped from the denominator) and (b) as firm positive when they carry a "missing" label.
   Two replacement sets: R13 = the 13 replacements that stand in for paywalled primary papers (10 within their stratum,
   3 filling the F1 shortfall), R14 = all 14 replacements (R13 plus the last-ranked F5 replacement, which stands in for
   the open-access primary 855947c638 that was ineligible at full text). Unweighted proportions with Wilson 95 % CI.
Weighted estimates use m3_summary_i139.stratified (fpc, logit CI; stratified bootstrap 2000 draws when a stratum is 0/1).
Usage: python -I tools/sens_sample_i139.py   -> extraction/sens_sample_i139.json
"""
import os, sys, csv, json, glob
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_summary_i139 as S

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
X = os.path.join(ROOT, "extraction")
rd = lambda p: list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"))
PILOT_IN_SAMPLE = {"f1ad3fb0ff": "B1", "79f1e21c4c": "B3"}
FAMS = ["F1_tyg", "F2_inflam", "F3_diet_life", "F4_lipid_anthro", "F5_other"]


def weighted(items, nh, elig):
    S.NHS.clear()
    S.NHS.update({h: nh[h] * elig[h][0] / elig[h][1] for h in nh})
    return S.stratified(items)


def main():
    sample = rd(os.path.join(X, "sample_47.tsv"))
    queue = rd(os.path.join(ROOT, "screening_full", "sampling_queue.tsv"))
    m1 = {r["rec_id"]: r for r in rd(os.path.join(X, "m1_results_i139.tsv"))}
    m2 = {r["rec_id"]: r for r in rd(os.path.join(X, "m2_results_i139.tsv"))}
    labels = {r["rec_id"]: r for r in rd(os.path.join(X, "labels_i139.tsv"))}
    m3 = {d["rec_id"]: d for d in (json.load(open(f, encoding="utf-8")) for f in glob.glob(os.path.join(X, "m3", "*.json")))}
    ids = [r["rec_id"] for r in sample]
    fam = {r["rec_id"]: r["family"] for r in sample}
    out = {}

    def e1(keep):
        k = [i for i in keep if m1[i]["status"] == "run"]
        nt = [i for i in keep if m1[i]["status"] != "run"]
        npass = sum(m1[i]["rank"] in ("strict_pass", "near_pass") for i in k)
        return {"i_testable": S.wilson(npass, len(k)), "ii_nt_fail": S.wilson(npass, len(keep)),
                "iii_nt_pass": S.wilson(npass + len(nt), len(keep)),
                "strict": sum(m1[i]["rank"] == "strict_pass" for i in k)}

    def e2(keep):
        v = [m2[i] for i in keep if i in m2 and m2[i]["E2_lt_080"] in ("True", "False")]
        return S.wilson(sum(r["E2_lt_080"] == "True" for r in v), len(v))

    def e3(keep, dag, nh, elig, positive_extra=()):
        prim = [i for i in keep if i in m3 and m3[i]["status"] == "run" and not m3[i]["flag_outcome_component"]]
        firm = [(fam[i], int(m3[i][f"E3a_{dag}"]["class"] == "firm_positive" or i in positive_extra)) for i in prim]
        point = [(fam[i], int(m3[i][f"E3a_{dag}"]["class"] in ("firm_positive", "point_positive") or i in positive_extra)) for i in prim]
        return {"n": len(prim), "firm_unweighted": S.wilson(sum(v for _, v in firm), len(prim)),
                "firm_weighted": weighted(firm, nh, elig), "point_weighted": weighted(point, nh, elig)}

    NH_MAIN = dict(S.NH)
    ELIG = dict(S.ELIG)
    # 1. pilot exclusion
    keep = [i for i in ids if i not in PILOT_IN_SAMPLE]
    nh_p = dict(NH_MAIN, F1_tyg=NH_MAIN["F1_tyg"] - 3)
    elig_p = dict(ELIG, F1_tyg=(ELIG["F1_tyg"][0] - 2, ELIG["F1_tyg"][1] - 2))
    out["pilot_exclusion"] = {"dropped": PILOT_IN_SAMPLE, "n_sample": len(keep), "E1": e1(keep), "E1_full": e1(ids),
                              "E2": e2(keep), "E2_full": e2(ids),
                              "E3a_A": e3(keep, "A", nh_p, elig_p), "E3a_B": e3(keep, "B", nh_p, elig_p),
                              "m1_rank_of_dropped": {i: m1[i]["rank"] for i in PILOT_IN_SAMPLE},
                              "N_h_used": nh_p, "elig_used": elig_p}
    # 2. N_h sensitivity
    final = {r["rec_id"]: r for r in rd(os.path.join(ROOT, "screening_full", "screening_final.tsv"))}
    passes = {"A": {r["rec_id"]: r["decision"] for r in rd(os.path.join(ROOT, "screening_full", "screener_A.tsv"))},
              "B": {r["rec_id"]: r["decision"] for r in rd(os.path.join(ROOT, "screening_full", "screener_B.tsv"))}}

    def count(dec, cats):
        return {h: sum(1 for r, v in dec.items() if final[r]["family"] == h and v in cats) for h in FAMS}
    variants = {"main_include_unclear": count({r: v["decision"] for r, v in final.items()}, {"include", "unclear"}),
                "include_only": count({r: v["decision"] for r, v in final.items()}, {"include"}),
                "pass_A_include_unclear": count(passes["A"], {"include", "unclear"}),
                "pass_B_include_unclear": count(passes["B"], {"include", "unclear"})}
    assert variants["main_include_unclear"] == NH_MAIN, (variants["main_include_unclear"], NH_MAIN)
    out["N_h"] = {k: {"N_h": v, "E3a_A": e3(ids, "A", v, ELIG), "E3a_B": e3(ids, "B", v, ELIG)} for k, v in variants.items()}
    # 3. replacement bounds
    repl = [r for r in queue if r["role"] == "replacement" and r["rec_id"] in fam]
    r14 = [r["rec_id"] for r in repl]
    f5 = sorted([r for r in repl if r["family"] == "F5_other"], key=lambda r: int(r["queue_pos"]))
    r13 = [i for i in r14 if i != f5[-1]["rec_id"]]
    rb = {}
    for name, R in (("R13", r13), ("R14", r14)):
        keep = [i for i in ids if i not in R]
        d = {"replacements": R}
        for dag in ("A", "B"):
            prim_all = [i for i in ids if i in m3 and m3[i]["status"] == "run" and not m3[i]["flag_outcome_component"]]
            miss_r = [i for i in R if i in prim_all and labels[i][f"label_{dag}"] == "missing"]
            d[dag] = {"as_not_testable": e3(keep, dag, NH_MAIN, ELIG)["firm_unweighted"],
                      "as_positive_if_missing": S.wilson(sum(m3[i][f"E3a_{dag}"]["class"] == "firm_positive" for i in prim_all) + len(miss_r), len(prim_all)),
                      "replacements_in_primary": [i for i in R if i in prim_all], "replacements_in_primary_missing": miss_r,
                      "all_replacements_missing_label": [i for i in R if labels[i][f"label_{dag}"] == "missing"]}
        rb[name] = d
    out["replacement_bounds"] = rb
    out["replacement_accounting"] = {"r14": r14, "f5_replacements_by_rank": [r["rec_id"] for r in f5], "r13": r13}
    json.dump(out, open(os.path.join(X, "sens_sample_i139.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
