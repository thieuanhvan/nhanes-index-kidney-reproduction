"""M1 builder, 71eedfeafb (OBS and CKD in overweight/obese adults; source withheld).
Spec: OBS per 1 unit, CKD (eGFR < 60 or ACR > 30), Model 3 (11 covariates), NHANES 1999-2018, survey weighted.
Sources in the paper (fulltext_txt/71eedfeafb.txt): population (Methods, Study population), OBS (Assessment of OBS),
outcome (Definition of CKD), covariates and medication classes (Covariates), models (Statistical analysis).
OBS from diet_indices_i139.DIET['71eedfeafb'] reading 1 (extraction/diet_recon_notes_i139.md).
"""
import re
import numpy as np
import pandas as pd

REC_ID = "71eedfeafb"
CYCLES = ["1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008", "2009-2010", "2011-2012", "2013-2014",
          "2015-2016", "2017-2018"]
OGTT = {"2005-2006": "OGTT_D", "2007-2008": "OGTT_E", "2009-2010": "OGTT_F", "2011-2012": "OGTT_G", "2013-2014": "OGTT_H",
        "2015-2016": "OGTT_I"}
RASI = (r"benazepril|captopril|enalapril|fosinopril|lisinopril|moexipril|perindopril|quinapril|ramipril|trandolapril|"
        r"losartan|valsartan|irbesartan|candesartan|olmesartan|telmisartan|eprosartan|azilsartan")
NSAID = (r"ibuprofen|naproxen|diclofenac|celecoxib|meloxicam|indomethacin|ketorolac|etodolac|nabumetone|piroxicam|"
         r"sulindac|oxaprozin|ketoprofen|flurbiprofen|diflunisal|salsalate|tolmetin|fenoprofen|mefenamic|meclofenamate|"
         r"rofecoxib|valdecoxib")


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
    df["x"] = lib.DIET[REC_ID](df, reading=1)
    e = df.egfr_ckdepi2009
    df["y"] = ((df.uacr > 30) | (e < 60)).astype(float).where(df.uacr.notna() & e.notna())

    # education: below high school vs high school and above
    df["educ_hs_plus"] = df.educ3.isin(["hs_ged", "gt_hs"]).astype(float).where(df.educ3.notna())
    # marital: married/living with partner vs not
    mar = _miss(lib.raw(df, lib.stems("DEMO"), "DMDMARTL"))
    df["married"] = mar.isin([1, 6]).astype(float).where(mar.isin([1, 2, 3, 4, 5, 6]))
    # hypertension: mean SBP >= 140, DBP >= 90, self-reported diagnosis, or current antihypertensive medication
    df["htn"] = df.hypertension_any_140
    # hyperlipidemia: TG >= 150, TC >= 200, LDL >= 130, HDL <= 40 men / <= 50 women, or lipid-lowering therapy
    bpq100d = _miss(lib.raw(df, lib.stems("BPQ"), "BPQ100D"))
    hdl_low = np.where(df.female == 1, df.hdl <= 50, df.hdl <= 40)
    pos = (df.tg >= 150) | (df.tc >= 200) | (df.ldl >= 130) | hdl_low | (bpq100d == 1)
    df["hyperlip"] = pd.Series(np.where(pos, 1.0, np.where(df.tc.notna() & df.hdl.notna(), 0.0, np.nan)), index=df.index)
    # diabetes: diagnosis, HbA1c >= 6.5, FPG >= 7.0 mmol/L, random glucose >= 11.1, OGTT >= 11.1, antidiabetic medication
    ogtt = lib.raw(df, OGTT, "LBXGLT")
    pos = ((df.DIQ010 == 1) | (df.hba1c >= 6.5) | (df.glucose >= 126) | (df.glucose_serum >= 200) | (ogtt >= 200)
           | (df.DIQ050 == 1) | (df.DIQ070 == 1))
    known = df.DIQ010.isin([1, 2, 3]) | df.hba1c.notna() | df.glucose.notna() | df.glucose_serum.notna()
    df["dm"] = pd.Series(np.where(pos, 1.0, np.where(known, 0.0, np.nan)), index=df.index)
    # ASCVD: CHD, angina, MI or stroke (MCQ160C, D, E, F)
    mc = df[["MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F"]]
    df["ascvd"] = pd.Series(np.where((mc == 1).any(axis=1), 1.0, np.where(mc.isin([1, 2]).all(axis=1), 0.0, np.nan)), index=df.index)
    df["rasi"] = _rx_any(lib, df, RASI)
    df["nsaid"] = _rx_any(lib, df, NSAID)

    terms = ["x"]
    terms += lib.dummies(df, "race1", "nh_white")
    terms += ["age", "female", "educ_hs_plus", "married", "htn", "hyperlip", "ascvd", "dm", "rasi", "nsaid"]

    df["w"] = lib.pooled_weight(df, "mec")
    dom = (df.age >= 20) & (df.bmi >= 25) & (df.pregnant != 1) & df.x.notna() & df.uacr.notna() & e.notna()
    return {"df": df, "domain": dom, "y": "y", "terms": terms, "target": "x", "weight": "w",
            "readings": {
                "effect measure label": "reading_1 OR (logistic regression), Table 2 header 'HR' ignored",
                "eGFR equation": "reading_1 CKD-EPI with race term, 2009 form ('updated' CKD-EPI with race)",
                "CKD threshold": "ACR > 30 (strict, as written) or eGFR < 60",
                "OBS": "DIET['71eedfeafb'] reading 1 (tertiles in age >= 20, BMI >= 25, not pregnant; 5 carotenoids; alcohol tertiles; >= 19 of 20 components)",
                "survey weighting": "weights not named; dietary weight is not a lib kind; smallest lib kind covering the sample = pooled MEC weight (WTMEC2YR, WTMEC4YR 1999-2002)",
                "race": "RIDRETH1 5 groups",
                "education": "binary below high school (DMDEDUC2 1-2) vs high school and above (3-5)",
                "marital status": "married/living with partner (DMDMARTL 1, 6) vs not (2-5)",
                "hypertension": "hypertension_any_140 (diagnosis, BPQ050A, mean SBP >= 140 or DBP >= 90), matches the paper",
                "hyperlipidemia": "as written; lipid-lowering therapy = BPQ100D now taking prescribed cholesterol medicine; 0 when TC and HDL known and no criterion met",
                "diabetes": "as written; random glucose = serum glucose LBXSGL >= 200 mg/dL; OGTT 2005-2016 only; medication = DIQ050 insulin or DIQ070 pills (DIQ070 absent 2005-2008)",
                "ASCVD": "MCQ160C, MCQ160D, MCQ160E, MCQ160F (no heart failure)",
                "RASI": "RXQ_RX generic names of ACE inhibitors and ARBs (combinations included)",
                "NSAIDs": "RXQ_RX generic names of non-salicylate NSAIDs; aspirin not counted",
                "age, sex": "age continuous, sex binary",
            },
            "notes": "builder written before first execution. Post-run technical fix 1 (after first m1_check): medication terms were all missing in 1999-2000 and 2001-2002 because lib.rx_any reads only RXDDRUG while those releases name the drug in RXD240B; replaced by the local _rx_any with the RXD240B fallback. No reading changed."}
