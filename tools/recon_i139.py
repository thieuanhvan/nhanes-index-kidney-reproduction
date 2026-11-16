"""I139 prereg section 4 and KD1: reconstruct each sampled paper's index on the paper's cycles and broad domain,
before any model is run. Usage: python -I tools/recon_i139.py
Reads the locked spec table (hash checked). Writes extraction/recon_i139.tsv (one row per paper) and prints the KD1 count.

Status per paper:
  computed            index computed from public CDC files with the paper's formula
  pending_diet        dietary score (DII, OBS, DOBS) not yet implemented here; components available in public dietary files
  not_reconstructable reason given (formula not in the paper, input not in public CDC files, ...)
Broad domain only (age limits, pregnancy, population restriction such as diabetes or hypertension); the paper's full
exclusion list and covariate completeness are applied in M1. N columns are therefore upper bounds, reported for
orientation, not used for KD1. No outcome model is fitted; the outcome inputs (creatinine, UACR) are only counted.
"""
import os, sys, json, csv, hashlib
import numpy as np, pandas as pd
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nhanes_core as nc
import nhanes_indices as ni
import nhanes_extra_i139 as nx

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
X = os.path.join(ROOT, "extraction")
TG, CH = ni.TG_MMOL, ni.CHOL_MMOL
UKAT = 1 / 60  # U/L to ukat/L


def s(d, c):
    return d[c].astype(float)


def div(a, b):
    return a / b.where(b > 0)


def tyg(d):
    return np.log(s(d, "tg") * s(d, "glucose") / 2)


def mets_vf(coef):
    def f(d):
        m = ni.mets_ir(d)
        return 4.466 + 0.011 * np.log(m) ** 3 + coef * np.log(ni.whtr(d)) ** 3 + 0.319 * (d.female == 0) + 0.594 * np.log(s(d, "age"))
    return f


def fsi(d):
    ratio = div(s(d, "alt"), s(d, "ast"))
    return (-7.981 + 0.011 * s(d, "age") - 0.146 * s(d, "female") + 0.173 * s(d, "bmi") + 0.007 * s(d, "tg")
            + 0.593 * s(d, "hypertension_any_140") + 0.789 * s(d, "diabetes_any") + 1.1 * (ratio >= 1.33).where(ratio.notna()))


INDEX = {  # key: (function, weight kind)
    "CTI_mgdl": (lambda d: 0.412 * np.log(s(d, "crp_mgL") / 10) + tyg(d), "fasting"),
    "TyG": (tyg, "fasting"),
    "NTI": (lambda d: ni.nlr(d) * tyg(d), "fasting"),
    "NPAR": (ni.npar, "mec"),
    "SIRI": (ni.siri, "mec"),
    "RAR": (lambda d: div(s(d, "rdw"), s(d, "albumin_serum")), "mec"),
    "HALP": (ni.halp, "mec"),
    "SII": (ni.sii, "mec"),
    "NHR": (lambda d: div(s(d, "neut_n"), s(d, "hdl") / CH), "mec"),
    "CALLY_gL": (lambda d: div(s(d, "albumin_serum") * 10 * s(d, "lymph_n"), s(d, "crp_mgL") * 10), "mec"),
    "CALLY_gdL": (lambda d: div(s(d, "lymph_n") * s(d, "albumin_serum"), s(d, "crp_mgL") / 10 * 10), "mec"),
    "AISI": (ni.aisi, "mec"),
    "ABSI_x100": (lambda d: ni.absi(d) * 100, "mec"),
    "PNI": (lambda d: 5 * s(d, "lymph_n") + s(d, "albumin_serum") * 10, "mec"),
    "WWI": (ni.wwi, "mec"),
    "UHR": (lambda d: div(s(d, "uric_acid"), s(d, "hdl")), "mec"),
    "UHR_x100": (lambda d: div(s(d, "uric_acid"), s(d, "hdl")) * 100, "mec"),
    "METS-VF_3.329": (mets_vf(3.329), "fasting"),
    "METS-VF_3.239": (mets_vf(3.239), "fasting"),
    "CMI": (ni.cmi, "fasting"),
    "RC": (lambda d: s(d, "tc") - s(d, "hdl") - s(d, "ldl"), "fasting"),
    "NHHR": (ni.non_hdl_hdl, "mec"),
    "ALBI": (lambda d: np.log10(s(d, "bilirubin_umol").where(lambda v: v > 0)) * 0.66 - 0.085 * s(d, "albumin_serum") * 10, "mec"),
    "FSI": (fsi, "fasting"),
    "FPSA_TPSA": (lambda d: div(s(d, "psa_free"), s(d, "psa_total")) * 100, "mec"),
    "CuZn": (lambda d: div(s(d, "copper"), s(d, "zinc")), "mec"),
    "FIB4_ukat": (lambda d: div(s(d, "age") * s(d, "ast") * UKAT, s(d, "platelets") * np.sqrt(s(d, "alt").where(lambda v: v > 0) * UKAT)), "mec"),
    "HRR": (lambda d: div(s(d, "hgb"), s(d, "rdw")), "mec"),
    # AVAI: formula in a MathML display formula of the JATS file (dropped by the text converter); HDL and TG read in mmol/L
    "AVAI_mmol": (lambda d: -16.186 - 1.369 * s(d, "hdl") / CH + 0.38 * s(d, "waist") + 0.144 * s(d, "age") - 0.013 * s(d, "bmi") - 0.151 * s(d, "tg") / TG, "fasting"),
    "TrunkFMR": (lambda d: div(s(d, "trunk_fat_g"), s(d, "trunk_lean_g")), "mec"),
    "SHR": (lambda d: div(s(d, "glucose_serum") / 18.0, 1.59 * s(d, "hba1c") - 2.59), "mec"),
}

# rec_id: (index key | status string, min age, max age, restriction list, outcome, eGFR equation)
DIAB, HTN = "diabetes", "hypertension"
CFG = {
    "f6e85f9043": ("CTI_mgdl", 20, None, [], "albuminuria", "2009"),
    "79f1e21c4c": ("TyG", 20, None, ["nopreg"], "ckd", "2009"),
    "f1ad3fb0ff": ("TyG", 20, None, [], "albuminuria", "2009"),
    "6b7f67827b": ("NTI", 20, None, [HTN], "ckd", "2009"),
    "d22fcf0a95": ("NPAR", 20, None, ["nopreg"], "ckd", "2009"),
    "94132d0c1a": ("SIRI", 20, None, ["nopreg"], "ckd", "2009"),
    "ea130293e0": ("RAR", 18, None, [DIAB], "ckd", "2009"),
    "fad26f18b1": ("HALP", 20, None, ["nopreg"], "ckd", "2009"),
    "cef539b61e": ("SII", 18, None, ["nopreg"], "ckd", "2009"),
    "9cd1df095e": ("NHR", 20, None, [], "albuminuria", "2009"),
    "53ca68a923": ("CALLY_gL", 20, None, ["nopreg"], "ckd", "2009"),
    "c62af9091f": ("CALLY_gdL", 20, None, ["nopreg"], "ckd", "2009"),
    "cf0849f17f": ("SII", 18, None, [DIAB, "nopreg"], "ckd", "2009"),
    "4b72270032": ("RAR", 20, None, [], "ckd", "2009"),
    "8991f10763": ("AISI", 20, None, [], "ckd", "2009"),
    "71eedfeafb": ("pending_diet:OBS", 20, None, ["bmi25", "nopreg"], "ckd", "2009"),
    "9e2e6087fa": ("pending_diet:DII", 20, None, ["nopreg"], "albuminuria", "2009"),
    "988df86a79": ("pending_diet:DII", 18, None, [], "ckd", "2009"),
    "30eae6044c": ("ABSI_x100", 20, 80, [DIAB, "nopreg"], "ckd", "2021"),
    "fe82b9da6a": ("pending_diet:OBS", 20, None, [DIAB], "ckd", "2009"),
    "aca6f864fb": ("pending_diet:OBS", 18, None, [], "ckd", "2009"),
    "78c5287646": ("PNI", 18, None, [DIAB, "nopreg"], "ckd", "2009"),
    "9835a67923": ("pending_diet:DOBS", 20, None, [DIAB, "nopreg"], "ckd", "2009"),
    "8a9e8d47ec": ("WWI", 20, None, [HTN], "albuminuria", "2009"),
    "b966d3cfdc": ("pending_diet:OBS", 18, None, [DIAB], "ckd", "2009"),
    "1c8c318fda": ("UHR", 45, None, [], "ckd", "2021"),
    "cde66dda6a": ("METS-VF_3.329", 20, None, [], "ckd", "2009"),
    "ea03123ab4": ("METS-VF_3.239", 20, None, [], "albuminuria", "2009"),
    "76d45e5c64": ("UHR_x100", 20, None, [DIAB], "ckd", "2009"),
    "e1d9d038bb": ("UHR_x100", 20, None, [DIAB, "nopreg"], "ckd", "2009"),
    "e76d3d1271": ("CMI", 20, None, [], "ckd", "2009"),
    "450ed6ebb6": ("RC", 20, None, [DIAB, "nopreg"], "ckd", "2009"),
    "69700c2b09": ("NHHR", 18, None, [DIAB, "nopreg"], "ckd", "2009"),
    "250e5761d5": ("AVAI_mmol", 20, None, [], "ckd", "2009"),
    "ff27565860": ("NHHR", 20, None, [DIAB], "ckd", "2009"),
    "259b9d6255": ("WWI", 20, None, ["hyperuricemia"], "ckd", "2009"),
    "fea693ee81": ("not_reconstructable:LC9 diet metric needs HEI-2015 from the USDA FPED files (third-party processed, prereg section 4)", 18, None, [], "ckd", "2009"),
    "ac96396902": ("ALBI", 20, None, [], "ckd", "2009"),
    "c48bec549c": ("FSI", 18, None, ["nopreg"], "albuminuria", "2009"),
    "ba21203b93": ("FPSA_TPSA", 20, None, [DIAB, "male"], "ckd", "2009"),
    "0717542568": ("CuZn", 20, None, [], "ckd", "2009"),
    "3865be04bf": ("FIB4_ukat", 18, None, [], "ckd", "2009"),
    "edaf91172c": ("HRR", 20, None, [], "ckd", "2009"),
    "2ecb3dff61": ("TrunkFMR", 20, None, [DIAB], "ckd", "2009"),
    "4de61a08ed": ("HRR", 20, None, [], "albuminuria", "2009"),
    "d84bd11e2e": ("CMI", 20, None, [DIAB, "nopreg"], "ckd", "2009"),
    "a906145f51": ("SHR", 20, None, [DIAB], "ckd", "2009"),
}


def cyc(c):
    return "2017-2020" if str(c).startswith("2017-2020") else str(c)[:9]


_tables = {}


def table(cycles):
    k = tuple(cycles)
    if k not in _tables:
        _tables[k] = nx.add_extra(ni.add_anthropometry(nc.build(list(cycles))))
    return _tables[k]


def domain(d, amin, amax, restr):
    m = s(d, "age") >= amin
    if amax:
        m &= s(d, "age") <= amax
    if "nopreg" in restr:
        m &= d.pregnant != 1
    if DIAB in restr:
        m &= d.diabetes_any == 1
    if HTN in restr:
        m &= d.hypertension_any_140 == 1
    if "bmi25" in restr:
        m &= s(d, "bmi") >= 25
    if "male" in restr:
        m &= d.female == 0
    if "hyperuricemia" in restr:
        m &= ((d.female == 0) & (s(d, "uric_acid") >= 7)) | ((d.female == 1) & (s(d, "uric_acid") >= 6))
    return m.fillna(False)


def main():
    lock = open(os.path.join(X, "spec_table_i139.sha256"), encoding="utf-8").read()
    pj = os.path.join(X, "spec_table_i139.json")
    assert hashlib.sha256(open(pj, "rb").read()).hexdigest() in lock, "spec table changed after lock"
    T = json.load(open(pj, encoding="utf-8"))
    rows = []
    for r in T:
        rid = r["rec_id"]
        key, amin, amax, restr, outc, eq = CFG[rid]
        cycles = [cyc(c) for c in r["cycles"]]
        row = {"family": r["family"], "rec_id": rid, "index_name": r["index_name"], "index_key": key, "cycles": " ".join(cycles),
               "n_reported": r["n_analytic"], "status": "", "reason": "", "n_domain": "", "n_index": "", "n_index_outcome": "",
               "ratio_to_reported": "", "cycles_without_index": "", "index_mean": "", "index_sd": ""}
        if ":" in key:
            row["status"], row["reason"] = key.split(":", 1)
        try:
            d = table(cycles)
        except Exception as e:
            row["status"], row["reason"] = "not_reconstructable", f"data load failed: {type(e).__name__} {e}"[:200]
            rows.append(row)
            continue
        dm = domain(d, amin, amax, restr)
        row["n_domain"] = int(dm.sum())
        if key in INDEX:
            f, _ = INDEX[key]
            v = f(d).replace([np.inf, -np.inf], np.nan)
            okx = dm & v.notna()
            e = d[f"egfr_ckdepi{eq}"]
            oko = d.uacr.notna() & (e.notna() if outc == "ckd" else True)
            row["n_index"] = int(okx.sum())
            row["n_index_outcome"] = int((okx & oko).sum())
            if row["n_reported"]:
                row["ratio_to_reported"] = int(10000 * row["n_index_outcome"] / row["n_reported"]) / 10000
            miss = [c for c in cycles if not ((d.cycle == c) & okx).any()]
            row["cycles_without_index"] = " ".join(miss)
            if okx.sum():
                row["index_mean"] = float("%.4g" % v[okx].mean())
                row["index_sd"] = float("%.4g" % v[okx].std())
            if okx.sum() == 0:
                row["status"], row["reason"] = "not_reconstructable", "index inputs absent in the public files of these cycles"
            else:
                row["status"] = "computed"
                if miss:
                    row["reason"] = "index absent in some cycles: " + " ".join(miss)
        elif row["status"] == "pending_diet":
            has = d.loc[dm, "has_diet1"].mean() if dm.any() else 0
            row["reason"] = f"day-1 dietary recall present for {int(1000 * has) / 10} % of the domain"
        rows.append(row)
        print(rid, row["status"], row["n_index_outcome"], row["n_reported"], row["ratio_to_reported"], row["reason"][:80], flush=True)
    p = os.path.join(X, "recon_i139.tsv")
    with open(p, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    from collections import Counter
    c = Counter(r["status"] for r in rows)
    print("STATUS", dict(c))
    print("KD1 lower bound (computed only):", c["computed"], " upper bound (+ pending diet):", c["computed"] + c["pending_diet"], " threshold 30")
    print("sha256", hashlib.sha256(open(p, "rb").read()).hexdigest())


if __name__ == "__main__":
    main()
