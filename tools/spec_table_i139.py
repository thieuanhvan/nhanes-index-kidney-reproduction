"""I139 prereg 3.2 step 7: build the locked specification table (first combination) from the reconciled extraction.
Usage: python -I tools/spec_table_i139.py
Reads extraction/sample_47.tsv and extraction/pass_R/<id>.json; writes extraction/spec_table_i139.tsv and .json,
and their sha256 into extraction/spec_table_i139.sha256 (with a UTC time stamp). No model is run before this hash exists.
"""
import os, sys, json, csv, hashlib, datetime
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
X = os.path.join(ROOT, "extraction")


def v(j, *path):
    for k in path:
        j = (j or {}).get(k) if isinstance(j, dict) else None
    return j.get("value") if isinstance(j, dict) and "value" in j else j


COLS = ["family", "rec_id", "doi", "index_name", "index_formula", "index_components", "cycles", "population", "n_analytic",
        "outcome_main", "outcome_definition", "egfr_equation", "uacr_threshold", "weight_variable", "cycle_pooling",
        "model_label", "effect_type", "exposure_scale", "estimate", "ci_low", "ci_high", "p", "source_table",
        "covariates_final", "covariate_sources_differ", "covariates_methods", "inflection_point", "ambiguities", "other_models"]


def main():
    S = list(csv.DictReader(open(os.path.join(X, "sample_47.tsv"), encoding="utf-8"), delimiter="\t"))
    rows = []
    for s in S:
        j = json.load(open(os.path.join(X, "pass_R", s["rec_id"] + ".json"), encoding="utf-8"))
        me = j.get("main_estimate", {})
        rows.append({
            "family": s["family"], "rec_id": j["rec_id"], "doi": j.get("doi"),
            "index_name": v(j, "index_name"), "index_formula": v(j, "index_formula"), "index_components": v(j, "index_components"),
            "cycles": v(j, "cycles"), "population": v(j, "population"), "n_analytic": v(j, "n_analytic"),
            "outcome_main": v(j, "outcome_main"), "outcome_definition": v(j, "outcome_definition"),
            "egfr_equation": v(j, "egfr_equation"), "uacr_threshold": v(j, "uacr_threshold"),
            "weight_variable": v(j, "weight_variable"), "cycle_pooling": v(j, "cycle_pooling"),
            "model_label": v(me, "model_label"), "effect_type": v(me, "effect_type"), "exposure_scale": v(me, "exposure_scale"),
            "estimate": v(me, "estimate"), "ci_low": v(me, "ci_low"), "ci_high": v(me, "ci_high"), "p": v(me, "p"),
            "source_table": v(me, "source_table"), "covariates_final": j.get("covariates_final"),
            "covariate_sources_differ": j.get("covariate_sources_differ"), "covariates_methods": v(j, "covariates_methods"),
            "inflection_point": v(j, "inflection_point"), "ambiguities": j.get("ambiguities") or [], "other_models": j.get("other_models") or [],
        })
    pj = os.path.join(X, "spec_table_i139.json")
    with open(pj, "w", encoding="utf-8") as fj:
        json.dump(rows, fj, ensure_ascii=False, indent=1)
    pt = os.path.join(X, "spec_table_i139.tsv")
    ft = open(pt, "w", encoding="utf-8", newline="")
    w = csv.writer(ft, delimiter="\t", lineterminator="\n")
    w.writerow(COLS)
    for r in rows:
        w.writerow([json.dumps(r[c], ensure_ascii=False) if isinstance(r[c], (list, dict)) else ("" if r[c] is None else r[c]) for c in COLS])
    ft.close()
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    with open(os.path.join(X, "spec_table_i139.sha256"), "w", encoding="utf-8") as f:
        for p in (pt, pj):
            f.write(f"{hashlib.sha256(open(p, 'rb').read()).hexdigest()}  {os.path.basename(p)}  locked {stamp}\n")
    print(open(os.path.join(X, "spec_table_i139.sha256"), encoding="utf-8").read())
    print("papers:", len(rows), "missing estimate:", sum(r["estimate"] is None for r in rows),
          "missing covariates:", sum(not r["covariates_final"] for r in rows), "missing N:", sum(r["n_analytic"] is None for r in rows))


if __name__ == "__main__":
    main()
