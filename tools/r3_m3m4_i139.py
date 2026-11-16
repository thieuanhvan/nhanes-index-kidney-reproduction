"""I139 round-3 review points M3 and M4 (paper6 main v3, section 7.3). Written 16/11/2026 after all main results were
seen; post hoc and descriptive, no model refit.

M3: secondary reproduction criterion that includes the confidence interval. A paper passes when it is a near pass under
    the registered M1 rule AND |log(width_repro / width_reported)| <= log 1.25, with width = log(upper) - log(lower) of
    the 95 % CI (no sample-size adjustment, since N already has to be within 15 %). Also reported: the CI criterion on
    its own over all 42 testable papers, and the same with the sqrt(N_repro / N_reported) adjustment used for E2.
M4: for each of the 17 M2 papers, the regime whose CI width is closest to the published one (|dlogW| only, ignoring
    the OR, which barely differs across regimes), and E2 restricted to strict-pass papers.
Usage: python -I tools/r3_m3m4_i139.py   -> extraction/r3_m3m4_i139.json
"""
import os, sys, csv, json, math, re
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_summary_i139 as S

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
X = os.path.join(ROOT, "extraction")
rd = lambda p: list(csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"))
L125 = math.log(1.25)
REG = ["i_correct_design", "ii_wrongtype_design", "iii_unweighted_model", "iv_weighted_model", "v_unweighted_design"]


def ci(s):
    m = re.match(r"^\s*([0-9.]+)\s*[-–]\s*([0-9.]+)\s*$", s or "")
    return (float(m.group(1)), float(m.group(2))) if m else None


def logw(c):
    return math.log(c[1]) - math.log(c[0]) if c and c[0] > 0 and c[1] > c[0] else None


def main():
    sample = {r["rec_id"]: r["family"] for r in rd(os.path.join(X, "sample_47.tsv"))}
    m1 = [r for r in rd(os.path.join(X, "m1_results_i139.tsv")) if r["rec_id"] in sample and r["status"] == "run"]
    assert len(m1) == 42
    rows, unparsed = [], []
    for r in m1:
        wr, wp = logw(ci(r["ci_reported"])), logw(ci(r["ci_repro"]))
        if wr is None or wp is None:
            unparsed.append([r["rec_id"], r["ci_reported"], r["ci_repro"]]); continue
        d = math.log(wp / wr)
        dn = d + 0.5 * math.log(float(r["n_repro"]) / float(r["n_reported"]))
        rows.append({"rec_id": r["rec_id"], "rank": r["rank"], "dlogW": S.tr(d), "dlogW_nadj": S.tr(dn),
                     "ci_ok": abs(d) <= L125, "ci_ok_nadj": abs(dn) <= L125})
    near = [x for x in rows if x["rank"] in ("strict_pass", "near_pass")]
    strict = [x for x in rows if x["rank"] == "strict_pass"]
    n_test = 42
    out = {"unparsed": unparsed,
           "M3_near_and_ci": {"k": sum(x["ci_ok"] for x in near), "near_pass_with_ci": len(near),
                               "E1_i_with_ci": S.wilson(sum(x["ci_ok"] for x in near), n_test),
                               "near_failing_ci": [[x["rec_id"], x["dlogW"]] for x in near if not x["ci_ok"]]},
           "M3_strict_and_ci": S.wilson(sum(x["ci_ok"] for x in strict), n_test),
           "M3_ci_only_all": S.wilson(sum(x["ci_ok"] for x in rows), len(rows)),
           "M3_ci_only_all_nadj": S.wilson(sum(x["ci_ok_nadj"] for x in rows), len(rows)),
           "rows": rows}
    m2 = rd(os.path.join(X, "m2_results_i139.tsv"))
    closest = []
    for r in m2:
        mm = {g: re.search(r"dlogW ([0-9.]+)", r[g]) for g in REG}
        if not all(mm.values()):  # published CI of zero width (dropped from E2 by amendment)
            closest.append({"rec_id": r["rec_id"], "m1_rank": r["m1_rank"], "width_ratio_E2": r["width_ratio_E2"],
                            "E2_lt_080": r["E2_lt_080"], "closest_by_width": "not_computable", "dlogW": {}})
            continue
        dw = {g: float(m.group(1)) for g, m in mm.items()}
        best = min(dw, key=dw.get)
        closest.append({"rec_id": r["rec_id"], "m1_rank": r["m1_rank"], "width_ratio_E2": r["width_ratio_E2"],
                        "E2_lt_080": r["E2_lt_080"], "closest_by_width": best, "dlogW": dw})
    tally = {}
    for c in closest:
        key = ("lt080" if c["E2_lt_080"] == "True" else "ge080_or_dropped")
        tally.setdefault(key, {}).setdefault(c["closest_by_width"], 0)
        tally[key][c["closest_by_width"]] += 1
    e2s = [c for c in closest if c["m1_rank"] == "strict_pass" and c["E2_lt_080"] in ("True", "False")]
    out["M4_closest_by_width"] = closest
    out["M4_tally"] = tally
    out["M4_E2_strict_only"] = S.wilson(sum(c["E2_lt_080"] == "True" for c in e2s), len(e2s))
    out["M4_E2_near_only"] = S.wilson(*(lambda v: (sum(c["E2_lt_080"] == "True" for c in v), len(v)))(
        [c for c in closest if c["m1_rank"] == "near_pass" and c["E2_lt_080"] in ("True", "False")]))
    json.dump(out, open(os.path.join(X, "r3_m3m4_i139.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    show = {k: v for k, v in out.items() if k not in ("rows", "M4_closest_by_width")}
    print(json.dumps(show, ensure_ascii=False, indent=1))
    for c in closest:
        print(c["rec_id"], c["m1_rank"], c["width_ratio_E2"], c["closest_by_width"], {k[:4]: v for k, v in c["dlogW"].items()})


if __name__ == "__main__":
    main()
