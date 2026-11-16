"""nhanes_indices: outcome-blind reconstruction of exposure indices for study I139 (prereg_i139_v7_en.md).

One function per index. Each takes the harmonised person table of nhanes_core (plus `height` and `weight`, added by
add_anthropometry()) and returns a pandas Series aligned to the table index. Missing inputs give NaN; a zero or
negative denominator gives NaN. No function reads, derives or models a kidney outcome.

Units of the harmonised table (NHANES as published): glucose, tg, hdl, tc mg/dL; insulin uU/mL; hba1c %;
waist, height cm; weight kg; bmi kg/m2; blood counts 1000 cells/uL; neut_pct %; hgb g/dL; albumin_serum g/dL.
Conversions used: TG mg/dL / 88.57 = mmol/L; cholesterol (HDL, TC) mg/dL / 38.67 = mmol/L; g/dL * 10 = g/L.

OUTCOME_COMPONENT flags (prereg "outcome-determined component"): True when the index contains creatinine,
cystatin C, BUN, urine albumin, uric acid, haemoglobin, serum albumin, phosphorus or potassium. Only NPAR and HALP
are True in this library. Each docstring repeats the flag.

Families (prereg): D1 TyG family, D2 inflammatory/haematological, D3 dietary/lifestyle (stubs), D4a lipid and
anthropometric, D4b insulin-based. DOIs checked on CrossRef on 2026-11-14 (all 28 resolved).

Outcome-blind loader: load_blind(cycles) builds the harmonised table without reading any kidney column
(creatinine, urine albumin/creatinine, BUN, uric acid); the derived eGFR/UACR/CKD columns are dropped.
"""
import os, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nhanes_core as nc

TG_MMOL = 88.57     # mg/dL per mmol/L, triglycerides
CHOL_MMOL = 38.67   # mg/dL per mmol/L, cholesterol

KIDNEY_VARS = ("scr_raw", "uric_acid", "bun", "u_alb", "u_cr", "uacr_cdc")
KIDNEY_DERIVED = ("scr", "scr_calibrated", "egfr_ckdepi2009", "egfr_ckdepi2021", "uacr", "albuminuria",
                  "egfr_lt60_2009", "egfr_lt60_2021", "ckd_2009", "ckd_2021")


# ---------------------------------------------------------------- loading (outcome-blind)
def add_anthropometry(df):
    """Add `height` (BMXHT, cm) and `weight` (BMXWT, kg) from the BMX file of each cycle, matched on cycle and SEQN."""
    df = df.copy(); df["height"] = np.nan; df["weight"] = np.nan
    for c in df.cycle.unique():
        x = nc.read(c, "BMX")
        if x is None: continue
        x = x.drop_duplicates("SEQN").set_index("SEQN")
        m = (df.cycle == c).to_numpy()
        for col, src in (("height", "BMXHT"), ("weight", "BMXWT")):
            if src in x.columns:
                df.loc[m, col] = x[src].reindex(df.loc[m, "SEQN"]).to_numpy(float)
    return df


def load_blind(cycles):
    """Harmonised table for `cycles` with the kidney columns never read from disk, plus height and weight."""
    saved = {k: nc.VARS[k] for k in KIDNEY_VARS}
    try:
        for k in KIDNEY_VARS: nc.VARS[k] = ("NONE", [], ())   # component absent from every file map: column is NaN
        df = nc.build(list(cycles))
    finally:
        nc.VARS.update(saved)
    df = df.drop(columns=[c for c in KIDNEY_VARS + KIDNEY_DERIVED if c in df.columns])
    return add_anthropometry(df)


# ---------------------------------------------------------------- helpers
def _s(df, col):
    if col not in df.columns:
        raise KeyError(f"column '{col}' missing; build the table with load_blind() or add_anthropometry()")
    return pd.to_numeric(df[col], errors="coerce").astype(float)

def _div(a, b):
    return a / b.where(b > 0)

def _female(df):
    return _s(df, "female")

def _log(x):
    return np.log(x.where(x > 0))

def _log10(x):
    return np.log10(x.where(x > 0))


# ================================================================ D1 TyG family
def tyg(df):
    """TyG = ln(TG [mg/dL] * fasting glucose [mg/dL] / 2). Unitless.
    Columns: tg, glucose (fasting subsample). OUTCOME_COMPONENT: False.
    Simental-Mendia LE, Rodriguez-Moran M, Guerrero-Romero F. Metab Syndr Relat Disord 2008;6(4):299-304.
    https://doi.org/[DOI withheld]"""
    return _log(_s(df, "tg") * _s(df, "glucose") / 2)

def tyg_bmi(df):
    """TyG-BMI = TyG * BMI [kg/m2]. OUTCOME_COMPONENT: False.
    Er LK, Wu S, Chou HH, et al. PLoS One 2016;11(3):e0149731. https://doi.org/[DOI withheld]"""
    return tyg(df) * _s(df, "bmi")

def tyg_wc(df):
    """TyG-WC = TyG * waist circumference [cm]. OUTCOME_COMPONENT: False.
    Zheng S, Shi S, Ren X, et al. J Transl Med 2016;14:260. https://doi.org/[DOI withheld]"""
    return tyg(df) * _s(df, "waist")

def tyg_whtr(df):
    """TyG-WHtR = TyG * waist [cm] / height [cm]. OUTCOME_COMPONENT: False.
    Lim J, Kim J, Koo SH, Kwon GC. PLoS One 2019;14(3):e0212963. https://doi.org/[DOI withheld]"""
    return tyg(df) * whtr(df)


# ================================================================ D2 inflammatory / haematological
def nlr(df):
    """NLR = neutrophils / lymphocytes (both 1000 cells/uL). Unitless. OUTCOME_COMPONENT: False.
    Forget P, Khalifa C, Defour JP, et al. BMC Res Notes 2017;10:12. https://doi.org/[DOI withheld]"""
    return _div(_s(df, "neut_n"), _s(df, "lymph_n"))

def plr(df):
    """PLR = platelets / lymphocytes (both 1000 cells/uL). Unitless. OUTCOME_COMPONENT: False.
    Smith RA, Bosonnet L, Raraty M, et al. Am J Surg 2009;197(4):466-472. https://doi.org/[DOI withheld]"""
    return _div(_s(df, "platelets"), _s(df, "lymph_n"))

def mlr(df):
    """MLR = monocytes / lymphocytes (both 1000 cells/uL); the reciprocal of the lymphocyte-to-monocyte ratio.
    OUTCOME_COMPONENT: False.
    Stotz M, Pichler M, Absenger G, et al. Br J Cancer 2014;110(2):435-440. https://doi.org/[DOI withheld]"""
    return _div(_s(df, "mono_n"), _s(df, "lymph_n"))

def sii(df):
    """SII = platelets * neutrophils / lymphocytes, counts in 1000 cells/uL (= 10^9/L). Unit 10^9/L.
    OUTCOME_COMPONENT: False.
    Hu B, Yang XR, Xu Y, et al. Clin Cancer Res 2014;20(23):6212-6222. https://doi.org/[DOI withheld]"""
    return _div(_s(df, "platelets") * _s(df, "neut_n"), _s(df, "lymph_n"))

def siri(df):
    """SIRI = neutrophils * monocytes / lymphocytes, counts in 10^9/L. Unit 10^9/L. OUTCOME_COMPONENT: False.
    Qi Q, Zhuang L, Shen Y, et al. Cancer 2016;122(14):2158-2167. https://doi.org/[DOI withheld]"""
    return _div(_s(df, "neut_n") * _s(df, "mono_n"), _s(df, "lymph_n"))

def aisi(df):
    """AISI = neutrophils * monocytes * platelets / lymphocytes, counts in 10^9/L. Unit (10^9/L)^2.
    OUTCOME_COMPONENT: False.
    Zinellu A, Collu C, Nasser M, et al. J Clin Med 2021;10(18):4134. https://doi.org/[DOI withheld]"""
    return _div(_s(df, "neut_n") * _s(df, "mono_n") * _s(df, "platelets"), _s(df, "lymph_n"))

def npar(df):
    """NPAR = neutrophil percentage [%] / serum albumin [g/dL]. Unitless as reported (typical adult value 12 to 16).
    Ambiguity: some papers use albumin in g/L and multiply by 100 (value = 10 * this one); this library uses g/dL,
    the scale of Wang et al. 2020 and of most NHANES reports.
    OUTCOME_COMPONENT: True (serum albumin).
    Wang B, Li D, Cheng B, Ying B, Gong Y. Biomed Res Int 2020;2020:5687672. https://doi.org/[DOI withheld]"""
    return _div(_s(df, "neut_pct"), _s(df, "albumin_serum"))

def halp(df):
    """HALP = haemoglobin [g/L] * albumin [g/L] * lymphocytes [10^9/L] / platelets [10^9/L].
    NHANES g/dL values are multiplied by 10 (original definition in g/L; papers that keep g/dL get 1/100 of this).
    OUTCOME_COMPONENT: True (haemoglobin, serum albumin).
    Chen XL, Xue L, Wang W, et al. Oncotarget 2015;6(38):41370-41382. https://doi.org/[DOI withheld]"""
    return _div(_s(df, "hgb") * 10 * _s(df, "albumin_serum") * 10 * _s(df, "lymph_n"), _s(df, "platelets"))


# ================================================================ D4a lipid / anthropometric
def aip(df):
    """AIP = log10(TG [mmol/L] / HDL-C [mmol/L]); TG mg/dL / 88.57, HDL mg/dL / 38.67. OUTCOME_COMPONENT: False.
    Dobiasova M, Frohlich J. Clin Biochem 2001;34(7):583-588. https://doi.org/[DOI withheld])00263-6"""
    return _log10(_div(_s(df, "tg") / TG_MMOL, _s(df, "hdl") / CHOL_MMOL))

def non_hdl_hdl(df):
    """Non-HDL-C/HDL-C ratio (NHHR) = (TC - HDL-C) / HDL-C, same units (unit-free). OUTCOME_COMPONENT: False.
    Kim SW, Jee JH, Kim HJ, et al. Int J Cardiol 2013;168(3):2678-2683. https://doi.org/[DOI withheld]"""
    h = _s(df, "hdl")
    return _div(_s(df, "tc") - h, h)

def vai(df):
    """VAI, sex-specific, TG and HDL-C in mmol/L, waist cm, BMI kg/m2:
    men   = WC / (39.68 + 1.88 * BMI) * (TG / 1.03) * (1.31 / HDL)
    women = WC / (36.58 + 1.89 * BMI) * (TG / 0.81) * (1.52 / HDL). OUTCOME_COMPONENT: False.
    Amato MC, Giordano C, Galia M, et al. Diabetes Care 2010;33(4):920-922. https://doi.org/[DOI withheld]"""
    wc, bmi, f = _s(df, "waist"), _s(df, "bmi"), _female(df)
    tg, h = _s(df, "tg") / TG_MMOL, _s(df, "hdl") / CHOL_MMOL
    men = wc / (39.68 + 1.88 * bmi) * (tg / 1.03) * _div(1.31 + 0 * h, h)
    wom = wc / (36.58 + 1.89 * bmi) * (tg / 0.81) * _div(1.52 + 0 * h, h)
    return pd.Series(np.where(f == 1, wom, np.where(f == 0, men, np.nan)), index=df.index)

def lap(df):
    """LAP = (WC [cm] - 65) * TG [mmol/L] in men, (WC - 58) * TG in women. No floor is applied, so a waist below
    the sex threshold gives a value of 0 or less (kept as computed). OUTCOME_COMPONENT: False.
    Kahn HS. BMC Cardiovasc Disord 2005;5:26. https://doi.org/[DOI withheld]"""
    f = _female(df)
    k = pd.Series(np.where(f == 1, 58.0, np.where(f == 0, 65.0, np.nan)), index=df.index)
    return (_s(df, "waist") - k) * _s(df, "tg") / TG_MMOL

def cmi(df):
    """CMI = (TG [mmol/L] / HDL-C [mmol/L]) * (waist / height). Ambiguity: the TG/HDL ratio depends on units
    (mg/dL ratio = 2.290 * mmol/L ratio); this library uses mmol/L, as most NHANES reports. OUTCOME_COMPONENT: False.
    Wakabayashi I, Daimon T. Clin Chim Acta 2015;438:274-278. https://doi.org/[DOI withheld]"""
    return _div(_s(df, "tg") / TG_MMOL, _s(df, "hdl") / CHOL_MMOL) * whtr(df)

def bri(df):
    """BRI = 364.2 - 365.5 * sqrt(1 - (WC / (2 pi))^2 / (0.5 * height)^2), waist and height in the same unit (cm).
    NaN when the term under the root is negative. OUTCOME_COMPONENT: False.
    Thomas DM, Bredlau C, Bosy-Westphal A, et al. Obesity 2013;21(11):2264-2271. https://doi.org/[DOI withheld]"""
    e = 1 - (_s(df, "waist") / (2 * np.pi)) ** 2 / (0.5 * _s(df, "height")) ** 2
    return 364.2 - 365.5 * np.sqrt(e.where(e >= 0))

def absi(df):
    """ABSI = WC [m] / (BMI^(2/3) * height [m]^(1/2)). Unit m^(11/6) kg^(-2/3); typical value about 0.08.
    OUTCOME_COMPONENT: False.
    Krakauer NY, Krakauer JC. PLoS One 2012;7(7):e39504. https://doi.org/[DOI withheld]"""
    bmi, ht = _s(df, "bmi"), _s(df, "height") / 100
    return _div(_s(df, "waist") / 100, bmi.where(bmi > 0) ** (2 / 3) * np.sqrt(ht.where(ht > 0)))

def wwi(df):
    """WWI = waist [cm] / sqrt(weight [kg]). Unit cm/sqrt(kg). OUTCOME_COMPONENT: False.
    Park Y, Kim NH, Kwon TY, Kim SG. Sci Rep 2018;8:16753. https://doi.org/[DOI withheld]"""
    w = _s(df, "weight")
    return _div(_s(df, "waist"), np.sqrt(w.where(w > 0)))

def whtr(df):
    """WHtR = waist [cm] / height [cm]. Unitless. OUTCOME_COMPONENT: False.
    Ashwell M, Hsieh SD. Int J Food Sci Nutr 2005;56(5):303-307. https://doi.org/[DOI withheld]"""
    return _div(_s(df, "waist"), _s(df, "height"))

def rfm(df):
    """RFM (%) = 64 - 20 * (height / waist) + 12 * female (female = 1, male = 0); height and waist in the same unit.
    OUTCOME_COMPONENT: False.
    Woolcott OO, Bergman RN. Sci Rep 2018;8:10980. https://doi.org/[DOI withheld]"""
    return 64 - 20 * _div(_s(df, "height"), _s(df, "waist")) + 12 * _female(df)


# ================================================================ D4b insulin-based
def homa_ir(df):
    """HOMA-IR = fasting insulin [uU/mL] * fasting glucose [mg/dL] / 405 (= glucose mmol/L * insulin / 22.5).
    No cross-assay insulin calibration (see nhanes_core_variables.md). OUTCOME_COMPONENT: False.
    Matthews DR, Hosker JP, Rudenski AS, et al. Diabetologia 1985;28(7):412-419. https://doi.org/[DOI withheld]"""
    return _s(df, "insulin") * _s(df, "glucose") / 405

def homa_b(df):
    """HOMA-B (%) = 360 * insulin [uU/mL] / (glucose [mg/dL] - 63) (= 20 * insulin / (glucose mmol/L - 3.5)).
    NaN when glucose <= 63 mg/dL. OUTCOME_COMPONENT: False.
    Matthews DR, Hosker JP, Rudenski AS, et al. Diabetologia 1985;28(7):412-419. https://doi.org/[DOI withheld]"""
    return _div(360 * _s(df, "insulin"), _s(df, "glucose") - 63)

def quicki(df):
    """QUICKI = 1 / (log10 insulin [uU/mL] + log10 glucose [mg/dL]). OUTCOME_COMPONENT: False.
    Katz A, Nambi SS, Mather K, et al. J Clin Endocrinol Metab 2000;85(7):2402-2410. https://doi.org/[DOI withheld]"""
    d = _log10(_s(df, "insulin")) + _log10(_s(df, "glucose"))
    return _div(1 + 0 * d, d)

def mets_ir(df):
    """METS-IR = ln(2 * glucose + TG) * BMI / ln(HDL-C), glucose, TG, HDL in mg/dL. NaN when HDL <= 1 mg/dL.
    OUTCOME_COMPONENT: False.
    Bello-Chavolla OY, Almeda-Valdes P, Gomez-Velasco D, et al. Eur J Endocrinol 2018;178(5):533-544.
    https://doi.org/[DOI withheld]"""
    return _div(_log(2 * _s(df, "glucose") + _s(df, "tg")) * _s(df, "bmi"), _log(_s(df, "hdl")))

def egdr(df, ht_col="hypertension_any_140"):
    """eGDR (mg/kg/min), waist-based formula of Williams et al. 2000:
    eGDR = 21.158 - 0.09 * waist [cm] - 3.407 * HT - 0.551 * HbA1c [%], HT = 1 if hypertensive else 0.
    HT default = nhanes_core `hypertension_any_140` (diagnosis, current medication, or mean BP >= 140/90); the
    WHR-based variant (24.31 - 12.22 * WHR - 3.29 * HT - 0.57 * HbA1c) cannot be computed (no hip circumference in
    the downloaded BMX columns used here). OUTCOME_COMPONENT: False.
    Williams KV, Erbey JR, Becker D, et al. Diabetes 2000;49(4):626-632. https://doi.org/[DOI withheld]"""
    return 21.158 - 0.09 * _s(df, "waist") - 3.407 * _s(df, ht_col) - 0.551 * _s(df, "hba1c")


# ================================================================ D3 dietary / lifestyle (stubs)
def life_simple_7(df):
    """Life's Simple 7 (AHA, 7 metrics: smoking, BMI, physical activity, healthy diet score, total cholesterol,
    blood pressure, fasting glucose). OUTCOME_COMPONENT: False.
    Lloyd-Jones DM, Hong Y, Labarthe D, et al. Circulation 2010;121(4):586-613.
    https://doi.org/[DOI withheld]
    STUB: smoking, BMI, TC, BP and glucose are in the harmonised table, but diet and physical activity are not."""
    raise NotImplementedError("Life's Simple 7 needs dietary recall DR1TOT/DR2TOT (healthy diet score), physical "
                              "activity PAQ (and PAQIAF for 1999-2006), and smoking SMQ (present); not downloaded.")

def life_essential_8(df):
    """Life's Essential 8 (AHA, 8 metrics scored 0-100: diet (HEI-2015), physical activity, nicotine exposure, sleep,
    BMI, non-HDL cholesterol, blood glucose or HbA1c, blood pressure). OUTCOME_COMPONENT: False.
    Lloyd-Jones DM, Allen NB, Anderson CAM, et al. Circulation 2022;146(5):e18-e43. https://doi.org/[DOI withheld]
    STUB: BMI, non-HDL, glucose/HbA1c and BP are in the harmonised table; the other four are not."""
    raise NotImplementedError("Life's Essential 8 needs DR1TOT/DR2TOT plus FPED/MPED for HEI-2015, PAQ (activity "
                              "minutes), SLQ (sleep hours, SLD010H/SLD012), SMQ/SMQRTU (nicotine, secondhand smoke), "
                              "and DIQ/BPQ/RXQ_RX for medication terms; not downloaded.")


# ---------------------------------------------------------------- registry
INDICES = {  # name: (family, function, OUTCOME_COMPONENT)
    "TyG": ("D1", tyg, False), "TyG-BMI": ("D1", tyg_bmi, False), "TyG-WC": ("D1", tyg_wc, False), "TyG-WHtR": ("D1", tyg_whtr, False),
    "NLR": ("D2", nlr, False), "PLR": ("D2", plr, False), "MLR": ("D2", mlr, False), "SII": ("D2", sii, False),
    "SIRI": ("D2", siri, False), "AISI": ("D2", aisi, False), "NPAR": ("D2", npar, True), "HALP": ("D2", halp, True),
    "AIP": ("D4a", aip, False), "NHHR": ("D4a", non_hdl_hdl, False), "VAI": ("D4a", vai, False), "LAP": ("D4a", lap, False),
    "CMI": ("D4a", cmi, False), "BRI": ("D4a", bri, False), "ABSI": ("D4a", absi, False), "WWI": ("D4a", wwi, False),
    "WHtR": ("D4a", whtr, False), "RFM": ("D4a", rfm, False),
    "HOMA-IR": ("D4b", homa_ir, False), "HOMA-B": ("D4b", homa_b, False), "QUICKI": ("D4b", quicki, False),
    "METS-IR": ("D4b", mets_ir, False), "eGDR": ("D4b", egdr, False),
}
STUBS = {"LS7": ("D3", life_simple_7, False), "LE8": ("D3", life_essential_8, False)}


def compute_all(df):
    """DataFrame of every implemented index (columns named as in INDICES)."""
    return pd.DataFrame({k: f(df) for k, (_, f, _) in INDICES.items()}, index=df.index)


# ---------------------------------------------------------------- descriptive sanity table (no outcome loaded)
def trunc4(x):
    """Format a number truncated (not rounded) to 4 significant digits."""
    from decimal import Decimal, ROUND_DOWN
    if x is None or not np.isfinite(x): return "NA"
    if x == 0: return "0"
    d = Decimal(repr(float(x))); e = d.adjusted()
    q = d.quantize(Decimal(1).scaleb(e - 3), rounding=ROUND_DOWN)
    return format(q, "f") if -6 < e < 9 else format(q, "E")

def _wq(x, w, q):
    o = np.argsort(x); c = np.cumsum(w[o]) / w.sum()
    return x[o][np.searchsorted(c, q)]

def write_summary(path=os.path.join(os.path.dirname(os.path.abspath(__file__)), "nhanes_indices_summary.md")):
    rows = []
    for c in nc.ALL_CYCLES:
        d = load_blind([c])
        w = nc.pooled_weight(d, "fasting")
        m = ((w > 0) & (d.age >= 20)).to_numpy()
        X = compute_all(d)[m]; W = w[m].to_numpy()
        for k in INDICES:
            v = X[k].to_numpy(float); ok = np.isfinite(v)
            miss = 100 * (1 - ok.mean())
            if ok.sum() < 2:
                rows.append((k, c, int(m.sum()), int(ok.sum()), miss, np.nan, np.nan, np.nan, np.nan)); continue
            x, ww = v[ok], W[ok]; mu = (ww * x).sum() / ww.sum(); sd = np.sqrt((ww * (x - mu) ** 2).sum() / ww.sum())
            rows.append((k, c, int(m.sum()), int(ok.sum()), miss, mu, sd, _wq(x, ww, 0.01), _wq(x, ww, 0.99)))
        print(c, int(m.sum()), flush=True)
    L = ["# nhanes_indices: descriptive sanity table", "",
         "Generated by `python -I tools/nhanes_indices.py --summary` (study I139, prereg_i139_v7_en.md). Descriptive only:",
         "the table is built with `load_blind()`, which never reads serum creatinine, urine albumin or creatinine, BUN or",
         "uric acid, and drops every derived kidney column (eGFR, UACR, CKD flags). No model is fitted.", "",
         "Domain: adults aged 20 or more in the fasting subsample of each release (fasting weight above 0). Each release is",
         "summarised on its own with its own fasting weight (WTSAF2YR, or WTSAFPRP for 2017-March 2020). 2017-2018 is a",
         "subset of 2017-March 2020 and is listed for completeness, not for pooling.", "",
         "Columns: N = adults in the domain; n = with the index computable; Miss % = unweighted share of N without the",
         "index; Mean and SD weighted (SD = weighted population SD); P1 and P99 weighted percentiles (first value whose",
         "cumulative weight share reaches the level). All numbers truncated to 4 significant digits.", "",
         "Height and weight come from BMX (BMXHT, BMXWT). Insulin is not calibrated across assays (see",
         "`nhanes_core_variables.md`). Units and formulas: docstrings of `tools/nhanes_indices.py`.", ""]
    fam = {k: v[0] for k, v in INDICES.items()}
    for k in INDICES:
        L += [f"## {k} ({fam[k]}{', outcome-determined component' if INDICES[k][2] else ''})", "",
              "| Release | N | n | Miss % | Mean | SD | P1 | P99 |", "|---|---|---|---|---|---|---|---|"]
        for r in rows:
            if r[0] != k: continue
            L.append(f"| {r[1]} | {r[2]} | {r[3]} | " + " | ".join(trunc4(v) for v in r[4:]) + " |")
        L.append("")
    L += ["## Not computed", "", "LS7 and LE8 (D3) are stubs: dietary recall, physical activity and sleep files were not downloaded.", ""]
    open(path, "w", encoding="utf-8").write("\n".join(L))


if __name__ == "__main__":
    if "--summary" in sys.argv: write_summary(); print("summary written")
