"""M1 builder, a906145f51 (stress hyperglycemia ratio and CKD in diabetes, Front Endocrinol 2024, doi [DOI withheld]).
Spec: SHR quintile 1 vs quintile 3 (reference), CKD (eGFR <= 60 or uACR >= 30), Model 3 (15 covariates),
NHANES 1999-March 2020, MEC weights pooled over 21.2 years.
Sources in the paper (fulltext_txt/a906145f51.txt): population line 14, SHR, DM and CKD line 17, covariates line 20,
weights lines 24-28, quintiles and final model lines 29-30, quintile limits line 33.
"""
import re
import numpy as np
import pandas as pd
import diet_indices_i139 as di

REC_ID = "a906145f51"
CYCLES = ["1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014",
          "2015-2016", "2017-2020"]
KIQ = {"2001-2002": "KIQ_U_B", "2003-2004": "KIQ_U_C", "2005-2006": "KIQ_U_D", "2007-2008": "KIQ_U_E", "2009-2010": "KIQ_U_F",
       "2011-2012": "KIQ_U_G", "2013-2014": "KIQ_U_H", "2015-2016": "KIQ_U_I", "2017-2020": "P_KIQ_U"}
ANTIDIAB = (r"metformin|glipizide|glyburide|glimepiride|glibenclamide|chlorpropamide|tolbutamide|tolazamide|acetohexamide|"
            r"pioglitazone|rosiglitazone|troglitazone|sitagliptin|saxagliptin|linagliptin|alogliptin|exenatide|liraglutide|"
            r"dulaglutide|semaglutide|lixisenatide|albiglutide|canagliflozin|dapagliflozin|empagliflozin|ertugliflozin|"
            r"repaglinide|nateglinide|acarbose|miglitol|pramlintide|insulin")
RASI = (r"benazepril|captopril|enalapril|fosinopril|lisinopril|moexipril|perindopril|quinapril|ramipril|trandolapril|"
        r"losartan|valsartan|irbesartan|candesartan|olmesartan|telmisartan|eprosartan|azilsartan")
STATIN = r"atorvastatin|simvastatin|rosuvastatin|pravastatin|lovastatin|fluvastatin|pitavastatin|cerivastatin"


def _miss(v, codes=(7, 9, 77, 99)):
    return v.where(~v.isin(codes))


def _rx_any(lib, df, pattern):
    """lib.rx_any with a fallback to RXD240B, the drug-name column of RXQ_RX in 1999-2002 (lib.rx_any reads only
    RXDDRUG and returns all-missing for 1999-2002). 1 = any matching prescription, 0 = person in RXQ_RX without one."""
    out = pd.Series(np.nan, index=df.index, dtype=float)
    rx = re.compile(pattern, re.I)
    for c in df.cycle.unique():
        st = ("RXQ_RX" + lib.SUF[c]) if c in lib.SUF else "P_RXQ_RX"
        x = lib._read(c, st)
        if x is None:
            continue
        col = "RXDDRUG" if "RXDDRUG" in x.columns else ("RXD240B" if "RXD240B" in x.columns else None)
        if col is None:
            continue
        names = x[col].map(lambda v: v.decode("latin-1") if isinstance(v, bytes) else str(v))
        hit = x.assign(h=names.str.contains(rx)).groupby("SEQN").h.any().astype(float)
        m = (df.cycle == c).to_numpy()
        out[m] = hit.reindex(df.loc[m, "SEQN"]).to_numpy(float)
    return out


def build(lib):
    df = lib.table(CYCLES)
    s = lib.stems
    df["x"] = lib.INDEX["SHR"][0](df).replace([np.inf, -np.inf], np.nan)
    e = df.egfr_ckdepi2009
    df["y"] = ((e <= 60) | (df.uacr >= 30)).astype(float).where(df.uacr.notna() & e.notna())

    # DM: HbA1c >= 6.5 %, 'doctor told you have diabetes' = yes, or 'taking insulin now' = yes
    df["dm"] = ((df.hba1c >= 6.5) | (df.DIQ010 == 1) | (df.DIQ050 == 1)).astype(float)
    dialysis = lib.raw(df, KIQ, "KIQ025") == 1
    liver = lib.raw(df, s("MCQ"), "MCQ160L") == 1
    cancer = lib.raw(df, s("MCQ"), "MCQ220") == 1

    # alcohol g/day from the day-1 24-h recall: none 0, moderate 0.1-27.9 men / 0.1-13.9 women, heavy >= 28 / >= 14
    alco = di.diet(df, "day1")["alco"]
    heavy = np.where(df.female == 1, alco >= 14, alco >= 28)
    df["alcohol"] = pd.Series(np.select([alco <= 0, heavy], ["none", "heavy"], default="moderate"),
                              index=df.index).astype(object).where(alco.notna(), np.nan)
    # DM duration = age - age at diagnosis (1999-2000 DIQ040G/DIQ040Q, 2001-2004 DID040G/DID040Q, 2005+ DID040);
    # under 1 year old = 0; 0 for diabetes without a diagnosis
    g = lib.raw(df, s("DIQ"), ["DIQ040G", "DID040G"])
    q = lib.raw(df, s("DIQ"), ["DIQ040Q", "DID040Q"])
    q = q.where(~q.isin([77777, 99999]))
    d5 = lib.raw(df, s("DIQ"), "DID040")
    d5 = d5.where(~d5.isin([777, 999])).replace(666, 0)
    old = df.cycle.isin(["1999-2000", "2001-2002", "2003-2004"])
    age_dx = pd.Series(np.where(old, np.where(g == 2, 0.0, np.where(g == 1, q, np.nan)), d5), index=df.index)
    dur = (df.age - age_dx).clip(lower=0)
    df["dm_dur"] = np.where(df.DIQ010 == 1, dur, np.where(df.DIQ010.isin([2, 3]), 0.0, np.nan)).astype(float)
    # anemia: Hb < 13 g/dL men, < 12 g/dL women (WHO)
    df["anemia"] = pd.Series(np.where(df.female == 1, df.hgb < 12, df.hgb < 13), index=df.index).astype(float).where(df.hgb.notna())
    df["antidiab"] = _rx_any(lib, df, ANTIDIAB)
    df["rasi"] = _rx_any(lib, df, RASI)
    df["statin"] = _rx_any(lib, df, STATIN)

    dom = ((df.dm == 1) & (df.age > 20) & ~dialysis & ~liver & ~cancer & df.x.notna() & df.uacr.notna() & e.notna())
    _, cuts = lib.quantile_dummies(df, "x", dom, k=5)
    grp = df["x_qgroup"]
    qn = []
    for j in (1, 2, 4, 5):
        df[f"x_q{j}v3"] = (grp == j).astype(float).where(grp.notna())
        qn.append(f"x_q{j}v3")

    terms = qn + ["age", "female", "bmi"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += lib.dummies(df, "smoking", "never")
    terms += lib.dummies(df, "alcohol", "none")
    terms += ["dm_dur", "hypertension_any_140", "cvd_history", "anemia", "uric_acid", "bun", "antidiab", "rasi", "statin"]

    df["w"] = lib.pooled_weight(df, "mec")
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x_q1v3", "weight": "w",
            "readings": {
                "main estimate": "reading_1 quintile 1 vs quintile 3, Model 3 OR 1.50",
                "continuous estimate": "reading_1 keep the category estimate",
                "weight type": "reading_1 MEC weights WTMEC4YR (4/21.2), WTMEC2YR (2/21.2), WTMECPRP (3.2/21.2) = lib pooled 'mec'",
                "eGFR equation": "not stated: CKD-EPI 2009 with race",
                "CKD": "eGFR <= 60 ('less than or equal' as written) or uACR >= 30",
                "SHR": "RBG = serum biochemistry glucose LBXSGL (mg/dL / 18), HbA1c LBXGH (recon INDEX SHR)",
                "quintiles": "unweighted quintiles of SHR cut on the domain (paper prints limits 0.707, 0.815, 0.926, 1.072); value at a cut point goes to the lower group",
                "age": "'age over 20 years' read literally as age > 20",
                "DM": "HbA1c >= 6.5, DIQ010 = 1, or DIQ050 insulin now = 1 (as written; pills not in the definition)",
                "dialysis": "KIQ025 = 1 excluded; KIQ_U absent in 1999-2000, so no dialysis exclusion there",
                "hepatic insufficiency": "MCQ160L ever told any liver condition = 1 excluded (closest public item)",
                "malignant tumor": "MCQ220 ever told cancer or malignancy = 1 excluded",
                "race": "RIDRETH1 5 groups",
                "smoking": "never / former / current (nhanes_core smoking), as defined",
                "alcohol": "g/day from the day-1 dietary recall (DR1TALCO / DRXTALCO, reliable recalls only); none / moderate / heavy by the stated sex cut points",
                "DM duration": "age minus self-reported age at diagnosis; 0 for diabetes without diagnosis",
                "hypertension": "not defined in the paper: diagnosis, current medication, or mean BP >= 140/90 (hypertension_any_140)",
                "CVD": "heart attack, CHF, angina, CHD, stroke (cvd_history, MCQ160B-F)",
                "anemia": "WHO: Hb < 13 g/dL men, < 12 g/dL women",
                "antidiabetic drugs, RASi, statins": "RXQ_RX generic-name regexes (antidiabetic incl. insulin; RASi = ACE inhibitors + ARBs)",
                "age, BMI, uric acid, BUN": "continuous",
                "pregnancy": "not excluded (not stated)",
            },
            "notes": "builder written before first execution. Post-run technical fix 1 (after first m1_check): medication terms were all missing in 1999-2000 and 2001-2002 because lib.rx_any reads only RXDDRUG while those releases name the drug in RXD240B; replaced by the local _rx_any with the RXD240B fallback. No reading changed."}
