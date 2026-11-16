"""I139 remaining descriptive sensitivity analyses (prereg Section 6 M3, "other sensitivity analyses"), 15/11/2026.
Re-runs the M3 comparison (m3_run_i139.run, B = 100, descriptive) under a modified NHANES table or label:
  egfr  : CKD-EPI 2009 and 2021 columns swapped (a paper using 2009 gets 2021 and vice versa)
  race  : binary race (one race term in the paper's model) treated as missing under DAG A and B (RIDRETH groups added)
  conf2 : albuminuria confirmed by the second urine sample (UACR = min of the two samples), 2009-2010 only;
          compared with conf1 = the same people, first sample only
uacr_down / uacr_up: UACR > 30 vs >= 30 differs only for people with UACR exactly 30 (378 person-papers, 41/42 papers);
  they are moved below or above 30 (B from env I139_OUTSENS_B, 20 used).
Usage: python -I tools/m3_outsens_i139.py <variant> <rec_id> [...]  -> extraction/m3_outsens/<variant>/<rec_id>.json
"""
import os, sys, re, json, traceback
import numpy as np, pandas as pd
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m23_common_i139 as C
import m3_run_i139 as M

VAR = sys.argv[1]
M.B = int(os.environ.get("I139_OUTSENS_B", "100"))
M.OUT = os.path.join(C.X, "m3_outsens", VAR)
os.makedirs(M.OUT, exist_ok=True)
lib = C.lib
_orig_table = lib.table
RACE_RX = re.compile(r"race|ethn|nhw|white|black|hispanic", re.I)


def _recode_kidney(df):
    ok = df.uacr.notna()
    df["albuminuria"] = (df.uacr >= 30).astype(float).where(ok)
    for eq in ("2009", "2021"):
        e = df[f"egfr_ckdepi{eq}"]
        df[f"ckd_{eq}"] = ((df.uacr >= 30) | (e < 60)).astype(float).where(ok & e.notna())
    return df


def _second_sample(df):
    x = lib._read("2009-2010", "ALB_CR_F")
    x = x.drop_duplicates("SEQN").set_index("SEQN")
    u2 = np.where(x.URXUCR2 > 0, x.URXUMA2 / x.URXUCR2 * 100, np.nan)
    return pd.Series(u2, index=x.index).reindex(df.SEQN).to_numpy(float)


def table(cycles):
    df = _orig_table(cycles)
    if VAR == "egfr":
        a, b = df.egfr_ckdepi2009.copy(), df.egfr_ckdepi2021.copy()
        df["egfr_ckdepi2009"], df["egfr_ckdepi2021"] = b, a
        df = _recode_kidney(df)
    elif VAR in ("uacr_down", "uacr_up"):   # people with UACR exactly 30 moved just below / just above 30
        tie = (df.uacr - 30).abs() < 1e-6
        df.loc[tie, "uacr"] = 29.999 if VAR == "uacr_down" else 30.001
        df = _recode_kidney(df)
    elif VAR in ("conf1", "conf2"):
        u2 = _second_sample(df)
        keep = (df.cycle == "2009-2010").to_numpy() & np.isfinite(u2) & df.uacr.notna().to_numpy()
        if VAR == "conf2":
            df.loc[keep, "uacr"] = np.minimum(df.loc[keep, "uacr"].to_numpy(float), u2[keep])
        df.loc[~keep, "uacr"] = np.nan        # only people with two samples in 2009-2010
        df = _recode_kidney(df)
    return df


lib.table = table


def main():
    T, L = C.spec_table(), C.labels()
    for rid in sys.argv[2:]:
        try:
            if VAR == "race":
                b = C.load_builder(rid)
                if b.get("not_testable"):
                    continue
                rt = [t for t in b["terms"] if RACE_RX.search(t)]
                if len(rt) != 1:
                    print(rid, "race terms", len(rt), "-> not binary, skipped", flush=True)
                    continue
                L = dict(L)
                lab = dict(L[rid])
                for k in ("missing_A", "missing_B"):
                    lab[k] = "|".join(sorted(set(filter(None, (lab[k] or "").split("|"))) | {"race"}))
                L[rid] = lab
                print(rid, "binary race term", rt, flush=True)
            res = M.run(rid, T, L)
        except Exception as e:
            traceback.print_exc()
            res = {"rec_id": rid, "status": "error", "error": f"{type(e).__name__}: {e}"[:400]}
        json.dump(res, open(os.path.join(M.OUT, rid + ".json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1, default=float)
        a = res.get("E3a_A", {})
        print(rid, VAR, res.get("status"), "nC", res.get("n_common"), "A", a.get("class"), a.get("delta"), flush=True)


if __name__ == "__main__":
    main()
