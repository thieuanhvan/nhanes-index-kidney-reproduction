"""I139 M1 builder check WITHOUT fitting: python -I tools/m1_check_i139.py <rec_id> [...]
Runs build(lib) and prints, per term, the non-missing count inside the domain, the number of complete cases,
and the outcome prevalence is NOT shown (the builder author must not tune on results). No model is fitted, no OR printed.
"""
import os, sys, importlib.util
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import m1_lib_i139 as lib

for rid in sys.argv[1:]:
    p = os.path.join(HERE, "m1_specs", rid + ".py")
    sp = importlib.util.spec_from_file_location("b_" + rid, p)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    b = m.build(lib)
    if b.get("not_testable"):
        print(rid, "NOT TESTABLE:", b["not_testable"])
        continue
    df, dom = b["df"], b["domain"].fillna(False).astype(bool)
    print(f"== {rid}: rows {len(df)}, domain {int(dom.sum())}, y defined in domain {int(df.loc[dom, b['y']].notna().sum())}, terms {len(b['terms'])}, target {b['target']}, weight {b['weight']}")
    cc = dom & df[[b["y"]] + b["terms"]].notna().all(axis=1)
    if b["weight"]:
        cc &= df[b["weight"]].fillna(0) > 0
    for t in b["terms"]:
        v = df.loc[dom, t]
        print(f"   {t:40s} non-missing {int(v.notna().sum())}  distinct {v.nunique()}")
    print(f"   complete cases (analysed n): {int(cc.sum())}")
