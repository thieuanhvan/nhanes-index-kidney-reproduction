"""M1 builder, 9835a67923 (DOBS and diabetic nephropathy; source withheld).
Spec: DOBS Q4 vs Q1, DN = diabetes with UACR > 30 mg/g, NHANES 2001-2010, Table 2 Model 3 before PSM,
dietary weights WTDR2D (two complete recalls) else WTDRD1, times 1/5.
Sources in fulltext_txt/9835a67923.txt: population (Study population, lines 22-23), DOBS (line 26), diabetes and DN
(line 29), covariates (line 32), inflammatory factors (line 35), weights (Statistical analysis, line 39),
Model 3 (Table 2 footnote).
The DOBS scoring rule is in an unread supplement; index = diet_indices_i139 dobs_9835a67923 reading 1
(tertile scores 0-2, not sex-specific; carotene alpha + beta), as the author chose.
"""
import numpy as np
import pandas as pd

REC_ID = "9835a67923"
CYCLES = ["2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010"]


def build(lib):
    import diet_indices_i139 as di   # shared module, read only (energy intake from the same recall rule)
    df = lib.table(CYCLES)
    s = lib.stems
    df["dobs"] = lib.DIET[REC_ID](df, reading=1)
    kcal = di.diet(df, "mean").kcal   # mean of two recalls when the second is complete, else day 1 (paper's rule)
    # diabetes: FBG >= 126, HbA1c >= 6.5, doctor told, 2-h OGTT >= 200, oral hypoglycaemic drugs or insulin
    ogtt = lib.raw(df, s("OGTT"), "LBXGLT")
    diq070 = lib.raw(df, s("DIQ"), ["DIQ070", "DID070"])
    df["diab"] = ((df.glucose >= 126) | (df.hba1c >= 6.5) | (df.DIQ010 == 1) | (ogtt >= 200)
                  | (diq070 == 1) | (df.DIQ050 == 1)).astype(float)
    df["y"] = (df.uacr > 30).astype(float).where(df.uacr.notna())
    # hypertension: SBP >= 130 or DBP >= 80 or on antihypertensive medication (BPQ050A now taking prescribed medicine)
    med = lib.raw(df, s("BPQ"), "BPQ050A")
    known = (df.sbp.notna() & df.dbp.notna()) | (med == 1)
    df["htn"] = ((df.sbp >= 130) | (df.dbp >= 80) | (med == 1)).astype(float).where(known)
    # marital status (DMDMARTL): married or living with partner / widowed, divorced or separated / never married
    m = lib.raw(df, s("DEMO"), "DMDMARTL")
    df["marital"] = pd.Series(np.select([m.isin([1, 6]), m.isin([2, 3, 4]), m == 5], ["married_partner", "wid_div_sep", "never"], default=""),
                              index=df.index).where(m.isin([1, 2, 3, 4, 5, 6]))
    dom = ((df.pregnant != 1) & (df.diab == 1) & df.uacr.notna() & df.egfr_ckdepi2009.notna()
           & df.dobs.notna() & kcal.between(500, 5000)
           & df[["wbc", "lymph_n", "mono_n", "neut_n", "crp_mgL"]].notna().all(axis=1))
    qnames, cuts = lib.quantile_dummies(df, "dobs", dom, k=4, prefix="dobs_q")
    terms = qnames + ["age", "female"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += lib.dummies(df, "marital", "married_partner")
    terms += ["bmi", "htn"]
    terms += lib.dummies(df, "smoking", "never")
    # dietary weight: WTDR2D when both recalls are complete (WTDR2D > 0), else WTDRD1; 2001-2002 has day 1 only
    dst = {c: ("DRXTOT_B" if c == "2001-2002" else ("DR1TOT_C" if c == "2003-2004" else "DR1TOT" + lib.SUF[c])) for c in CYCLES}
    w1 = lib.raw(df, dst, "WTDRD1")
    w2 = lib.raw(df, dst, "WTDR2D")
    df["w"] = (w2.where(w2 > 0, w1) / 5).fillna(0)
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "dobs_q4", "weight": "w",
            "readings": {
                "exposure scale of the main estimate": "reading_1: Q4 vs Q1, Model 3, before PSM, OR 0.78",
                "DOBS scoring": "diet_indices_i139 reading_1: tertile scores 0-2 (reversed for fat and iron), not sex-specific, carotene alpha + beta, mean of two recalls (supplement unread; author's choice)",
                "quartile cut": "unweighted quartiles of DOBS in the domain",
                "age": "no age limit stated; none applied (marital status, asked of adults, limits complete cases)",
                "diabetes": "paper's five criteria; oral drugs DIQ070 (DID070 in 2005-2008); OGTT exists 2005-2010 only",
                "outcome": "UACR > 30 mg/g as written",
                "eGFR": "CKD-EPI 2009 with race; used only for the 'missing eGFR' exclusion",
                "energy": "mean daily energy of the recalls used for the DOBS, 500-5000 kcal kept",
                "inflammatory factors": "exclusion for missing WBC, lymphocyte, monocyte, neutrophil counts or CRP",
                "race": "RIDRETH1 5 groups (coding not stated)",
                "marital status": "3 groups: married or living with partner; widowed, divorced or separated; never married (coding not stated, common NHANES grouping)",
                "BMI, age": "continuous",
                "hypertension": "SBP >= 130 or DBP >= 80 (mean of readings) or current prescribed antihypertensive medication (BPQ050A); no self-reported diagnosis, as written",
                "smoking": "never, former, current",
                "weight": "WTDR2D if > 0 (two complete recalls) else WTDRD1, divided by 5 (2001-2002: WTDRD1)",
            },
            "notes": f"DOBS quartile cut points in the domain: {[int(c * 1000) / 1000 for c in cuts]}"}
