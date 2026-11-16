"""I139 prereg 5.3 and 5.4: assign the DAG family from index components, then label each paper's covariate set.
Usage: python -I tools/labels_i139.py
Reads the locked extraction/spec_table_i139.json (hash checked against spec_table_i139.sha256).
Writes extraction/labels_i139.tsv, extraction/unmapped_covariates_i139.tsv, extraction/labels_i139.sha256,
and the author-audit list extraction/author_audit_list_i139.tsv (prereg 3.2 step 10: every paper missing at least one
confounder under DAG A or B, plus 10 others chosen with seed 20261114-i139-audit).
"""
import os, sys, json, csv, hashlib, re
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from concepts_i139 import concepts, concept_set

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
X = os.path.join(ROOT, "extraction")

# component classes for 5.3
C_TG = r"\btg\b|triglycer"
C_GLU = r"glucose|\bfpg\b|\brbg\b|hba1c|insulin"
C_ANTH = r"\bbmi\b|body mass|\bwc\b|waist|height|weight|fat mass|muscle mass"
C_LIPID = r"\btg\b|triglycer|hdl|ldl|\btc\b|cholesterol|lipid"
C_INFL = r"neutrophil|lymphocyte|monocyte|platelet|\bcrp\b|\brdw\b|hemoglobin|haemoglobin|\bhgb\b|albumin|white blood"
C_DIET = r"diet|vitamin|fib(er|re)|fat\b|fatty|pufa|sfa|mufa|carotene|carotenoid|riboflavin|niacin|folate|folic|thiamine|calcium|magnesium|\bmg\b|zinc|copper|selenium|iron|\bfe\b|protein|carbohydrate|caffeine|energy|cholesterol|lycopene|lutein|cryptoxanthin|alcohol|physical activity|cotinine|nicotine|sleep"


def family(components):
    s = [str(c).lower() for c in components]
    has = lambda rx: any(re.search(rx, c) for c in s)
    if has(C_TG) and has(C_GLU):
        return "D1"
    if has(C_GLU) or has(C_ANTH):
        return "D4b"
    if all(re.search(C_LIPID, c) for c in s):
        return "D4a"
    if has(C_INFL):
        return "D2"
    if has(C_DIET):
        return "D3"
    return "D5"


# node table 5.2: concept -> rank per family. G = DAG A and B; B = DAG B only (B+); M = mediator; O = outcome descendant; X = collider
G_ALL = ["age", "sex", "race", "education", "pir", "smoking"]
B_ALL = ["alcohol", "physical_activity"]


def node_ranks(fam):
    g = list(G_ALL)
    b = list(B_ALL)
    m = []
    if fam in ("D1", "D2", "D4a", "D5"):
        b.append("bmi_waist")
    if fam == "D3":
        m.append("bmi_waist")
    if fam == "D4b":
        m.append("bmi_waist")
    if fam in ("D2", "D3"):
        g += ["diabetes_dx", "hypertension_dx"]
        m.append("glycaemia_measured")
    if fam in ("D1", "D4b"):
        m += ["diabetes_dx", "glycaemia_measured", "hypertension_dx"]
    if fam == "D4a":
        b.append("diabetes_dx")
        m.append("hypertension_dx")
    if fam in ("D1", "D4a"):
        b.append("lipid_drugs")
    return g, b, m


O_NODES = ["bp_measured", "uric_acid", "kidney_measure"]
X_NODES = ["cvd"]


def comp_concepts(components):
    cc = concept_set(components)
    for c in components:
        s = str(c).lower()
        if re.search(C_GLU, s):
            cc.add("glycaemia_measured")
        if re.search(C_LIPID, s):
            cc.add("lipids")
        if re.search(C_ANTH, s):
            cc.add("bmi_waist")
        if re.search(r"nicotine|cotinine", s):
            cc.add("smoking")
        if re.search(r"blood pressure|hypertension", s):
            cc.add("bp_measured")
    return cc


def main():
    lock = open(os.path.join(X, "spec_table_i139.sha256"), encoding="utf-8").read()
    pj = os.path.join(X, "spec_table_i139.json")
    h = hashlib.sha256(open(pj, "rb").read()).hexdigest()
    assert h in lock, "spec table changed after lock"
    T = json.load(open(pj, encoding="utf-8"))
    out, unm = [], []
    for r in T:
        comps = r["index_components"] or []
        fam = family(comps)
        g, b, m = node_ranks(fam)
        cc = comp_concepts(comps)
        cov = r["covariates_final"] or []
        if isinstance(cov, dict):  # some reconciled files keep the field form {value, quote, loc}
            cov = cov.get("value") or []
        have = concept_set(cov)
        for x in cov:
            if not concepts(x):
                unm.append((r["rec_id"], x))
        # A node that is constant in the analytic population cannot be "missing" (deviation, logged in verdict 9.11):
        # diabetes-only, hypertension-only or one-sex populations, read by code from the population field.
        pop = str(r["population"] or "").lower()
        const = set()
        if re.search(r"with diabetes|diabetes patients|diabetic participants|diabetes participants|with t2dm|with type 2 diabetes|diagnosed with diabetes|patients with diabetes|non-diabet|non-t2dm|without t2dm|no t2dm|^dm and|restricted to the .* t2dm", pop):
            const.add("diabetes_dx")
        if re.search(r"hypertensive (adults|patients)", pop):
            const.add("hypertension_dx")
        if re.search(r"\bmale patients|\bmen only", pop):
            const.add("sex")
        kc = lambda c: c in cc or c in const  # index component or constant: do not adjust, never "missing"
        miss_A = [c for c in g if c not in have and not kc(c)]
        miss_B = miss_A + [c for c in b if c not in have and not kc(c)]
        adj_comp = sorted(c for c in have & cc if c not in ("cycle",))
        adj_med = sorted(c for c in have if c in m and not kc(c))
        adj_out = sorted(c for c in have if c in O_NODES and not kc(c))
        adj_col = sorted(c for c in have if c in X_NODES)
        out.append({"family": r["family"], "rec_id": r["rec_id"], "dag_family": fam, "index_name": r["index_name"],
                    "n_covariates": len(cov), "concepts_present": "|".join(sorted(have)),
                    "constant_in_population": "|".join(sorted(const)), "label_A": "missing" if miss_A else "sufficient", "missing_A": "|".join(miss_A),
                    "label_B": "missing" if miss_B else "sufficient", "missing_B": "|".join(miss_B),
                    "adjusts_component": "|".join(adj_comp), "adjusts_mediator": "|".join(adj_med),
                    "conditions_on_outcome": "|".join(adj_out), "conditions_on_collider": "|".join(adj_col)})
    pl = os.path.join(X, "labels_i139.tsv")
    fl = open(pl, "w", encoding="utf-8", newline="")
    w = csv.DictWriter(fl, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(out)
    fl.close()
    pu = os.path.join(X, "unmapped_covariates_i139.tsv")
    fu = open(pu, "w", encoding="utf-8", newline="")
    wu = csv.writer(fu, delimiter="\t", lineterminator="\n")
    wu.writerow(["rec_id", "covariate_as_written"])
    wu.writerows(unm)
    fu.close()
    # author audit list (3.2 step 10)
    scope_a = [o["rec_id"] for o in out if o["label_A"] == "missing" or o["label_B"] == "missing"]
    rest = sorted((o["rec_id"] for o in out if o["rec_id"] not in scope_a),
                  key=lambda x: hashlib.sha256(f"20261114-i139-audit:{x}".encode()).hexdigest())[:10]
    pa = os.path.join(X, "author_audit_list_i139.tsv")
    fa = open(pa, "w", encoding="utf-8", newline="")
    wa = csv.writer(fa, delimiter="\t", lineterminator="\n")
    wa.writerow(["rec_id", "scope", "label_A", "missing_A", "label_B", "missing_B"])
    L = {o["rec_id"]: o for o in out}
    for x in scope_a:
        wa.writerow([x, "a_missing", L[x]["label_A"], L[x]["missing_A"], L[x]["label_B"], L[x]["missing_B"]])
    for x in rest:
        wa.writerow([x, "b_random", L[x]["label_A"], L[x]["missing_A"], L[x]["label_B"], L[x]["missing_B"]])
    fa.close()
    with open(os.path.join(X, "labels_i139.sha256"), "w", encoding="utf-8") as f:
        for p in (pl, pu, pa, os.path.join(ROOT, "tools", "concepts_i139.py"), os.path.abspath(__file__)):
            f.write(f"{hashlib.sha256(open(p, 'rb').read()).hexdigest()}  {os.path.basename(p)}\n")
    from collections import Counter
    print("DAG families:", dict(Counter(o["dag_family"] for o in out)))
    print("label A missing:", sum(o["label_A"] == "missing" for o in out), "/", len(out),
          " label B missing:", sum(o["label_B"] == "missing" for o in out))
    print("missing A concepts:", dict(Counter(c for o in out for c in o["missing_A"].split("|") if c)))
    print("missing B concepts:", dict(Counter(c for o in out for c in o["missing_B"].split("|") if c)))
    print("adjusts mediator:", sum(bool(o["adjusts_mediator"]) for o in out), " conditions on outcome:", sum(bool(o["conditions_on_outcome"]) for o in out),
          " collider:", sum(bool(o["conditions_on_collider"]) for o in out), " component:", sum(bool(o["adjusts_component"]) for o in out))
    print("unmapped covariate names:", len(unm), " audit list:", len(scope_a), "+", len(rest))


if __name__ == "__main__":
    main()
