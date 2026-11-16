"""diet_indices_i139: outcome-blind reconstruction of the dietary/lifestyle indices (D3 family) of the 7 I139 papers.

One function per paper, keyed by rec_id in DIET. Each takes the harmonised person table (nhanes_indices.load_blind,
which never reads a kidney column) and a `reading` (1 = reading_1, 2 = reading_2; see the per-function docstring and
extraction/diet_recon_notes_i139.md) and returns a pandas Series aligned to the table index (NaN = not computable).
No function reads creatinine, urine albumin, eGFR or CKD; no outcome model is fitted.

Data: public CDC/NCHS XPT files in data-external/nhanes-extra/<cycle>/. Day-1 totals DR1TOT_* (2003+), DRXTOT(_B)
(1999-2002, one day only); day-2 totals DR2TOT_* (2003+). A recall is used only when its status is 1 (reliable and met
minimum criteria: DRDDRSTS 1999-2000, DRDDRSTZ 2001-2002, DR1DRSTZ/DR2DRSTZ 2003+).
Third-party processed files (USDA FPED) are not used (prereg section 4).

DII constants (global mean, global SD, overall inflammatory effect score; Shivappa N et al., Public Health Nutr
2014;17(8):1689-1696, doi [DOI withheld]) are copied from the open-source R package `dietaryindex`
(J. Zhan), file R/DII.R, https://raw.githubusercontent.com/jamesjiadazhan/dietaryindex/main/R/DII.R, fetched
2026-11-15, last commit touching the file f478683d1e7130b731c53ad5b078f39481051161 (2024-09-29), file sha256
c6392632713152cb27b48781e64f1105a1396f24f7a7ea1231b262072543dc2a. The package's NHANES mapping
(R/DII_NHANES_FPED.R) is followed for the column choice: n-3 = P183+P184+P205+P225+P226, n-6 = P182+P204,
vitamin A = VARA (RAE), vitamin E = ATOC, folic acid = FA, caffeine mg / 1000 = g. Percentile step: 2*Phi(z) - 1.

Physical activity: GPAQ items (2007-2018), MET-min/week = MET x days x minutes, MET 8 vigorous work, 4 moderate work,
4 walk/bicycle transport, 8 vigorous recreation, 4 moderate recreation (NHANES PAQ analytic notes). 1999-2006: leisure
from the individual activity file PAQIAF (past 30 days) plus walk/bicycle transport from PAQ (see _paq_9906); work
minutes were not collected, so the 1999-2006 total is leisure + transport.

Usage: python -I diet_indices_i139.py   (writes extraction/diet_recon_i139.tsv)
"""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nhanes_core as nc            # sets below-normal priority on import
import nhanes_extra_i139 as nx
import nhanes_indices as ni

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(os.path.dirname(HERE), "extraction")

nx.STEMS["DIET2"] = {"2003-2004": "DR2TOT_C", **{c: "DR2TOT" + s for c, s in nx.SUF.items() if s > "_C"}}

# ---------------------------------------------------------------- dietary totals
# generic name -> column suffix after DR1T/DR2T/DRXT (first one present is used)
NUTR = {"kcal": ["KCAL"], "prot": ["PROT"], "carb": ["CARB"], "fat": ["TFAT"], "sfa": ["SFAT"], "mufa": ["MFAT"],
        "pufa": ["PFAT"], "chol": ["CHOL"], "fiber": ["FIBE"], "vita_rae": ["VARA"], "vb1": ["VB1"], "vb2": ["VB2"],
        "niac": ["NIAC"], "vb6": ["VB6"], "fola": ["FOLA"], "fa": ["FA"], "vb12": ["VB12"], "vc": ["VC"], "vd": ["VD"],
        "ve": ["ATOC", "VE"],      # 1999-2000 has only DRXTVE (mg alpha-tocopherol equivalents)
        "calc": ["CALC"], "magn": ["MAGN"], "iron": ["IRON"], "zinc": ["ZINC"], "copp": ["COPP"], "sele": ["SELE"],
        "caff": ["CAFF"], "alco": ["ALCO"], "acar": ["ACAR"], "bcar": ["BCAR"], "cryp": ["CRYP"], "lyco": ["LYCO"],
        "lz": ["LZ"], "p182": ["P182"], "p183": ["P183"], "p184": ["P184"], "p204": ["P204"], "p205": ["P205"],
        "p225": ["P225"], "p226": ["P226"]}
STATUS = ["DR1DRSTZ", "DR2DRSTZ", "DRDDRSTZ", "DRDDRSTS"]


def _day(cycle, comp, seqn):
    x = nx.read_stem(cycle, comp)
    out = pd.DataFrame(np.nan, index=seqn, columns=list(NUTR))
    if x is None:
        return out
    pre = next((p for p in ("DR1T", "DR2T", "DRXT") if any(c.startswith(p) for c in x.columns)), None)
    st = next((c for c in STATUS if c in x.columns), None)
    ok = x[st] == 1 if st else pd.Series(True, index=x.index)
    for k, sfx in NUTR.items():
        col = next((pre + s for s in sfx if pre + s in x.columns), None)
        if col is not None:
            out[k] = x[col].where(ok).reindex(seqn).to_numpy(float)
    return out


def diet(df, days="mean"):
    """Nutrient totals for every row of df. days='day1' (first recall) or 'mean' (mean of the two recalls when the
    second is valid, otherwise day 1). Adds derived sums n3, n6, carotene_ab (alpha+beta), carotenoids (5 named)."""
    parts = []
    for c in df.cycle.unique():
        m = (df.cycle == c).to_numpy()
        seqn = df.loc[m, "SEQN"].to_numpy()
        d1 = _day(c, "DIET1", seqn)
        if days == "mean":
            d2 = _day(c, "DIET2", seqn)
            has2 = d2["kcal"].notna() & d1["kcal"].notna()
            d1 = d1.mask(has2, (d1 + d2) / 2, axis=0)
        d1.index = df.index[m]
        parts.append(d1)
    d = pd.concat(parts).reindex(df.index)
    d["n3"] = d[["p183", "p184", "p205", "p225", "p226"]].sum(axis=1, min_count=5)
    d["n6"] = d[["p182", "p204"]].sum(axis=1, min_count=2)
    d["carotene_ab"] = d[["acar", "bcar"]].sum(axis=1, min_count=2)
    d["carotenoids"] = d[["acar", "bcar", "cryp", "lyco", "lz"]].sum(axis=1, min_count=5)
    return d


# ---------------------------------------------------------------- lifestyle components
def _paq_9906(c, seqn):
    """1999-2006 physical activity (MET-min/week) from PAQ + PAQIAF (past 30 days, times x minutes, x 7/30).
    Leisure: PAQIAF rows; MET by level (vigorous PADLEVEL 2 = 8, moderate 1 = 4: the GPAQ weights the papers cite)
    or by the activity's own PADMETS. 0 when PAD200 and PAD320 are both answered 2/3 (no leisure activity).
    Transport: PAD020 walk/bicycle; times per week from PAQ050Q/PAQ050U (1 day, 2 week, 3 month), PAD080 min, MET 4.
    Work minutes are not collected in 1999-2006, so the total is leisure + transport."""
    x = nx.read_stem(c, "PAQ").reindex(seqn)
    stem = nx.STEMS["PAQIAF"].get(c)
    a = pd.read_sas(os.path.join(nc.ROOT, c, stem + ".xpt"), format="xport")
    a.columns = [k.upper() for k in a.columns]; a["SEQN"] = a["SEQN"].astype(np.int64)
    t = a.PADTIMES.where(a.PADTIMES.between(1, 3000)); du = a.PADDURAT.where(a.PADDURAT.between(1, 1440))
    lvl = np.where(a.PADLEVEL == 2, 8.0, np.where(a.PADLEVEL == 1, 4.0, np.nan))
    a["lvl"] = lvl * t * du * 7 / 30; a["mets"] = a.PADMETS * t * du * 7 / 30
    g = a.groupby("SEQN")[["lvl", "mets"]].sum(min_count=1)
    g = g.where(a.assign(bad=a[["lvl", "mets"]].isna().any(axis=1)).groupby("SEQN").bad.sum() == 0)
    none = x.PAD200.isin([2, 3]) & x.PAD320.isin([2, 3])
    any_ = x.PAD200.isin([1, 2, 3]) & x.PAD320.isin([1, 2, 3])
    lei = {}
    for k in ("lvl", "mets"):
        v = g[k].reindex(seqn)
        lei[k] = pd.Series(np.where(none, 0.0, np.where(any_, v, np.nan)), index=seqn)
    mult = x.PAQ050U.map({1: 7.0, 2: 1.0, 3: 7.0 / 30})
    q = x.PAQ050Q.where(x.PAQ050Q.between(1, 1000)); mn = x.PAD080.where(x.PAD080.between(1, 1440))
    tr = pd.Series(np.where(x.PAD020.isin([2, 3]), 0.0, np.where(x.PAD020 == 1, 4.0 * q * mult * mn, np.nan)), index=seqn)
    return pd.DataFrame({"pa_total": lei["lvl"] + tr, "pa_total_mets": lei["mets"] + tr, "pa_leisure": lei["lvl"]})


def _paq(df):
    out = pd.DataFrame(index=df.index, columns=["pa_total", "pa_total_mets", "pa_leisure"], dtype=float)
    dom = {"vw": ("PAQ605", "PAQ610", "PAD615", 8.0), "mw": ("PAQ620", "PAQ625", "PAD630", 4.0),
           "tr": ("PAQ635", "PAQ640", "PAD645", 4.0), "vr": ("PAQ650", "PAQ655", "PAD660", 8.0),
           "mr": ("PAQ665", "PAQ670", "PAD675", 4.0)}
    for c in df.cycle.unique():
        x = nx.read_stem(c, "PAQ")
        if x is None:
            continue
        m = (df.cycle == c).to_numpy()
        if "PAQ650" not in x.columns:                    # 1999-2006: PAQ + individual activity file
            if c in nx.STEMS["PAQIAF"]:
                r = _paq_9906(c, df.loc[m, "SEQN"].to_numpy())
                for k in out.columns:
                    out.loc[m, k] = r[k].to_numpy()
            continue
        x = x.reindex(df.loc[m, "SEQN"])
        met = {}
        for k, (yn, dd, mm, w) in dom.items():
            days = x[dd].where(x[dd].between(1, 7)); mins = x[mm].where(x[mm].between(1, 1440))
            met[k] = np.where(x[yn] == 2, 0.0, np.where(x[yn] == 1, w * days * mins, np.nan))
        met = pd.DataFrame(met, index=df.index[m])
        out.loc[m, "pa_total"] = met.sum(axis=1, min_count=5).to_numpy()
        out.loc[m, "pa_total_mets"] = out.loc[m, "pa_total"]
        out.loc[m, "pa_leisure"] = met[["vr", "mr"]].sum(axis=1, min_count=2).to_numpy()
    return out


def lifestyle(df):
    paq = _paq(df)
    return pd.DataFrame({"pa_total": paq.pa_total, "pa_total_mets": paq.pa_total_mets, "pa_leisure": paq.pa_leisure,
                         "cotinine": nx._pull(df, "COT", ["LBXCOT"]), "bmi": ni._s(df, "bmi")}, index=df.index)


# ---------------------------------------------------------------- scoring helpers
def tertile_score(x, ref, sex=None, reverse=False):
    """0/1/2 by tertile; cut points = 1/3 and 2/3 sample quantiles (linear interpolation, R type 7) of x among rows
    with ref True (per sex when sex is given). x <= q1 -> 0, q1 < x <= q2 -> 1, x > q2 -> 2; reversed for
    pro-oxidants. Ties at a cut point (e.g. many zero alcohol intakes) fall in the lower tertile."""
    s = pd.Series(np.nan, index=x.index)
    groups = [slice(None)] if sex is None else [sex == g for g in (1, 2)]
    for g in groups:
        gm = pd.Series(True, index=x.index) if isinstance(g, slice) else g
        r = x[gm & ref & x.notna()]
        if len(r) < 3:
            continue
        q1, q2 = np.quantile(r, [1 / 3, 2 / 3])
        xx = x[gm]
        sc = np.where(xx <= q1, 0.0, np.where(xx <= q2, 1.0, 2.0))
        s[gm] = np.where(xx.notna(), sc, np.nan)
    return 2 - s if reverse else s


def alcohol_cat(alco, female, scores=(2, 1, 0)):
    """Non-drinker (0 g/d), non-heavy (>0 and < 15 g/d women, < 30 g/d men), heavy (>= 15 women, >= 30 men)."""
    heavy = np.where(female == 1, 15.0, 30.0)
    v = np.where(alco <= 0, scores[0], np.where(alco < heavy, scores[1], scores[2])).astype(float)
    return pd.Series(np.where(alco.notna() & female.notna(), v, np.nan), index=alco.index)


def _sum_min(scores, k):
    """Row sum of available component scores when at least k components are present, else NaN."""
    n = scores.notna().sum(axis=1)
    return scores.sum(axis=1).where(n >= k)


OBS16 = [("fiber", +1), ("carotene", +1), ("vb2", +1), ("niac", +1), ("vb6", +1), ("fola", +1), ("vb12", +1),
         ("vc", +1), ("ve", +1), ("calc", +1), ("magn", +1), ("zinc", +1), ("copp", +1), ("sele", +1),
         ("fat", -1), ("iron", -1)]


# ---------------------------------------------------------------- DII
# Shivappa 2014 Table (via dietaryindex R/DII.R, see module docstring): name -> (effect score, global mean, global SD)
DII_STD = {
    "alco": (-0.278, 13.98, 3.72),        # g
    "vb12": (0.106, 5.15, 2.7),           # ug
    "vb6": (-0.365, 1.47, 0.74),          # mg
    "bcar": (-0.584, 3718, 1720),         # ug
    "caff": (-0.11, 8.05, 6.67),          # g
    "carb": (0.097, 272.2, 40),           # g
    "chol": (0.11, 279.4, 51.2),          # mg
    "kcal": (0.18, 2056, 338),            # kcal
    "fat": (0.298, 71.4, 19.4),           # g
    "fiber": (-0.663, 18.8, 4.9),         # g
    "fa": (-0.19, 273, 70.7),             # ug folic acid
    "iron": (0.032, 13.35, 3.71),         # mg
    "magn": (-0.484, 310.1, 139.4),       # mg
    "mufa": (-0.009, 27, 6.1),            # g
    "niac": (-0.246, 25.9, 11.77),        # mg
    "n3": (-0.436, 1.06, 1.06),           # g
    "n6": (-0.159, 10.8, 7.5),            # g
    "prot": (0.021, 79.4, 13.9),          # g
    "pufa": (-0.337, 13.88, 3.76),        # g
    "vb2": (-0.068, 1.7, 0.79),           # mg
    "sfa": (0.373, 28.6, 8),              # g
    "sele": (-0.191, 67, 25.1),           # ug
    "vb1": (-0.098, 1.7, 0.66),           # mg
    "vita_rae": (-0.401, 983.9, 518.6),   # RE in Shivappa; NHANES RAE ug used (as in dietaryindex)
    "vc": (-0.424, 118.2, 43.46),         # mg
    "vd": (-0.446, 6.26, 2.21),           # ug
    "ve": (-0.419, 8.73, 1.49),           # mg
    "zinc": (-0.313, 9.84, 2.19),         # mg
}
DII28 = list(DII_STD)
DII27_9E2E = [k for k in DII28 if k != "vb6"]   # the 27 parameters 9e2e6087fa names (vitamin B6 is not named)


def dii(d, comps, require_all=False):
    """Shivappa DII: sum over comps of effect * (2*Phi((intake - mean)/sd) - 1). Intake columns of d (diet()),
    caffeine converted mg -> g. A parameter absent from a whole cycle (vitamin D before 2007-2008) is skipped
    (dietaryindex behaviour) unless require_all."""
    tot = pd.Series(0.0, index=d.index); valid = d["kcal"].notna()
    for k in comps:
        eff, mu, sd = DII_STD[k]
        v = d[k] / 1000 if k == "caff" else d[k]
        term = eff * (2 * norm.cdf((v - mu) / sd) - 1)
        term = pd.Series(term, index=d.index)
        if require_all:
            valid &= v.notna()
        tot += term.fillna(0.0)
    return tot.where(valid)


# ---------------------------------------------------------------- per paper
def obs_71eedfeafb(df, reading=1):
    """OBS, 20 components (16 dietary mean of 2 recalls + PA, alcohol, cotinine, BMI), sex-specific tertiles 0-2,
    pro-oxidants (fat, iron, cotinine, alcohol, BMI) reversed, valid >= 19 of 20, sum of available scores.
    reading 1: tertiles in the paper's outcome-free population (age >= 20, BMI >= 25, not pregnant);
               'carotenoids' = alpha + beta carotene + beta-cryptoxanthin + lycopene + lutein/zeaxanthin;
               alcohol (g/d, recall) scored by sex-specific tertiles like the other 19 components.
    reading 2: tertiles in all participants of the pooled cycles; carotenoids = beta-carotene only;
               alcohol scored 2/1/0 non/non-heavy/heavy drinker."""
    d = diet(df, "mean"); L = lifestyle(df); sex = df["sex"]
    if reading == 1:
        ref = (df.age >= 20) & (df.bmi >= 25) & (df.pregnant != 1); car = d["carotenoids"]
    else:
        ref = pd.Series(True, index=df.index); car = d["bcar"]
    comp = {k: d[k] for k, _ in OBS16 if k != "carotene"}; comp["carotene"] = car
    sc = {k: tertile_score(comp[k], ref, sex, reverse=(s < 0)) for k, s in OBS16}
    sc["pa"] = tertile_score(L.pa_total, ref, sex)
    sc["cot"] = tertile_score(L.cotinine, ref, sex, reverse=True)
    sc["bmi"] = tertile_score(L.bmi, ref, sex, reverse=True)
    sc["alc"] = tertile_score(d["alco"], ref, sex, reverse=True) if reading == 1 else alcohol_cat(d["alco"], df.female)
    return _sum_min(pd.DataFrame(sc), 19)


def dii_9e2e6087fa(df, reading=1):
    """DII. reading 1: the 27 parameters the paper names, day-1 recall ('during the day'); vitamin D skipped in
    2005-2006 (not in the file). reading 2: 28 parameters (adds vitamin B6) and mean of the two recalls."""
    if reading == 1:
        return dii(diet(df, "day1"), DII27_9E2E)
    return dii(diet(df, "mean"), DII28)


def dii_988df86a79(df, reading=1):
    """DII, 28 named parameters, first-day recall (stated). reading 2: mean of the two recalls."""
    return dii(diet(df, "day1" if reading == 1 else "mean"), DII28)


def obs_fe82b9da6a(df, reading=2):
    """NOT RECONSTRUCTABLE under prereg section 4: the 16 dietary components are not named (only 'OBS comprises
    16 dietary and 4 lifestyle components', citation [20]). reading 2 (sensitivity only): the 16 components of the
    common OBS list (as aca6f864fb names them), day-1 recall, carotene = alpha + beta, tertiles not sex-specific in
    diabetes_any adults >= 20, PA MET-min/week tertiles, cotinine and BMI tertiles reversed, alcohol 2/1/0 by the
    paper's categories, all 20 required."""
    if reading == 1:
        return pd.Series(np.nan, index=df.index)
    d = diet(df, "day1"); L = lifestyle(df)
    ref = (df.age >= 20) & (df.diabetes_any == 1)
    comp = {k: d[k] for k, _ in OBS16 if k != "carotene"}; comp["carotene"] = d["carotene_ab"]
    sc = {k: tertile_score(comp[k], ref, None, reverse=(s < 0)) for k, s in OBS16}
    sc["pa"] = tertile_score(L.pa_total, ref)
    sc["cot"] = tertile_score(L.cotinine, ref, reverse=True)
    sc["bmi"] = tertile_score(L.bmi, ref, reverse=True)
    sc["alc"] = alcohol_cat(d["alco"], df.female)
    return _sum_min(pd.DataFrame(sc), 20)


def obs_aca6f864fb(df, reading=1):
    """OBS, 16 dietary + PA + BMI + alcohol + cotinine. Sex-specific tertiles (0-2, reversed for fat, iron,
    cotinine) for the 16 nutrients, PA (total MET-min/week) and cotinine; alcohol 2/1/0 (non/moderate/heavy);
    BMI 2/1/0 (<25 / 25-<30 / >=30). Score kept when at least 16 components are present (exclusion 'items < 16').
    reading 1: mean of the recalls ('24-h dietary recall interviews'); tertiles among adults >= 18;
               carotene = alpha + beta carotene; folate = total folate (DR1TFOLA).
    reading 2: day-1 recall; tertiles among all participants; carotene = beta-carotene only."""
    d = diet(df, "mean" if reading == 1 else "day1"); L = lifestyle(df); sex = df["sex"]
    ref = (df.age >= 18) if reading == 1 else pd.Series(True, index=df.index)
    comp = {k: d[k] for k, _ in OBS16 if k != "carotene"}
    comp["carotene"] = d["carotene_ab"] if reading == 1 else d["bcar"]
    sc = {k: tertile_score(comp[k], ref, sex, reverse=(s < 0)) for k, s in OBS16}
    sc["pa"] = tertile_score(L.pa_total_mets, ref, sex)   # 'MET of each physical activity' (PADMETS in 2005-2006)
    sc["cot"] = tertile_score(L.cotinine, ref, sex, reverse=True)
    b = L.bmi
    sc["bmi"] = pd.Series(np.where(b >= 30, 0.0, np.where(b >= 25, 1.0, 2.0)), index=df.index).where(b.notna())
    sc["alc"] = alcohol_cat(d["alco"], df.female)
    return _sum_min(pd.DataFrame(sc), 16)


def dobs_9835a67923(df, reading=1):
    """DOBS, 16 dietary components (14 antioxidants, fat and iron pro-oxidants), mean of the two recalls (stated),
    tertile scores 0-2 summed (scoring table is in Supplementary Table 2, not in the main text); all 16 required.
    reading 1: tertiles not sex-specific, computed in the outcome-free population (diabetes_any, not pregnant,
               energy 500-5000 kcal/day); carotene = alpha + beta carotene.
    reading 2: sex-specific tertiles in the same population; carotene = beta-carotene only."""
    d = diet(df, "mean")
    ref = (df.diabetes_any == 1) & (df.pregnant != 1) & d.kcal.between(500, 5000)
    sex = None if reading == 1 else df["sex"]
    comp = {k: d[k] for k, _ in OBS16 if k != "carotene"}
    comp["carotene"] = d["carotene_ab"] if reading == 1 else d["bcar"]
    sc = {k: tertile_score(comp[k], ref, sex, reverse=(s < 0)) for k, s in OBS16}
    return _sum_min(pd.DataFrame(sc), 16)


# b966dc: Table 5 cut points (male, female): antioxidant -> (a, b): <a = 1, [a, b) = 2, >= b = 3;
# pro-oxidant -> (a, b): >= b = 1, [a, b) = 2, < a = 3
B966_CUTS = {
    "vc": (+1, (28.90, 84.50), (26.90, 77.80)), "acar": (+1, (25.00, 152.00), (21.00, 115.00)),
    "bcar": (+1, (434.00, 1614.00), (377.00, 1411.00)), "cryp": (+1, (12.00, 64.00), (12.00, 58.00)),
    "lyco": (+1, (113.00, 3712.00), (144.00, 2891.00)), "lz": (+1, (501.00, 1132.00), (381.00, 931.00)),
    "ve": (+1, (5.01, 8.95), (4.01, 7.44)), "zinc": (+1, (7.73, 12.81), (6.45, 10.33)),
    "sele": (+1, (83.90, 130.70), (67.60, 104.10)),
    "fat": (-1, (56.74, 93.85), (44.53, 73.87)), "pufa": (-1, (11.81, 20.50), (9.63, 16.91)),
    "sfa": (-1, (17.22, 29.48), (13.62, 23.51)), "iron": (-1, (10.72, 16.61), (8.86, 13.83)),
    "pa": (+1, (40.00, 480.00), (40.00, 240.00)), "bmi": (-1, (28.03, 33.26), (29.10, 35.20)),
    "cot": (-1, (0.02, 0.20), (0.02, 0.08)),
}


def obs_b966d3cfdc(df, reading=1):
    """OBS, 13 dietary + 4 lifestyle, scores 1-3 with the paper's own sex-specific cut points (Table 5);
    alcohol 3/2/1 non/non-heavy/heavy; all 17 required.
    reading 1: day-1 recall ('a 24-hour dietary recall interview'); PA = leisure-time MET-min/week (vigorous 8,
               moderate 4).
    reading 2: mean of the two recalls; PA = total MET-min/week over all GPAQ domains."""
    d = diet(df, "day1" if reading == 1 else "mean"); L = lifestyle(df); fem = df.female
    vals = {**{k: d[k] for k in B966_CUTS if k in d.columns}, "bmi": L.bmi, "cot": L.cotinine,
            "pa": L.pa_leisure if reading == 1 else L.pa_total}
    sc = {}
    for k, (dirn, cm, cf) in B966_CUTS.items():
        x = vals[k]
        a = np.where(fem == 1, cf[0], cm[0]); b = np.where(fem == 1, cf[1], cm[1])
        lvl = np.where(x < a, 1.0, np.where(x < b, 2.0, 3.0))
        if dirn < 0:
            lvl = 4.0 - lvl
        sc[k] = pd.Series(lvl, index=df.index).where(x.notna() & fem.notna())
    sc["alc"] = alcohol_cat(d["alco"], fem, scores=(3, 2, 1))
    return _sum_min(pd.DataFrame(sc), 17)


DIET = {"71eedfeafb": obs_71eedfeafb, "9e2e6087fa": dii_9e2e6087fa, "988df86a79": dii_988df86a79,
        "fe82b9da6a": obs_fe82b9da6a, "aca6f864fb": obs_aca6f864fb, "9835a67923": dobs_9835a67923,
        "b966d3cfdc": obs_b966d3cfdc}

CYC = {"71eedfeafb": nc.ALL_CYCLES[:10], "9e2e6087fa": nc.ALL_CYCLES[3:10], "988df86a79": nc.ALL_CYCLES[6:10],
       "fe82b9da6a": nc.ALL_CYCLES[6:8], "aca6f864fb": nc.ALL_CYCLES[3:10], "9835a67923": nc.ALL_CYCLES[1:6],
       "b966d3cfdc": nc.ALL_CYCLES[4:10]}

# components each paper needs (generic names; lifestyle ones checked separately)
NEEDS = {"71eedfeafb": ["fiber", "carotenoids", "vb2", "niac", "vb6", "vb12", "vc", "ve", "fola", "calc", "magn",
                        "zinc", "copp", "sele", "fat", "iron", "alco", "pa_total", "cotinine", "bmi"],
         "9e2e6087fa": DII28, "988df86a79": DII28,
         "fe82b9da6a": [k if k != "carotene" else "carotene_ab" for k, _ in OBS16] + ["alco", "pa_total", "cotinine", "bmi"],
         "aca6f864fb": [k if k != "carotene" else "carotene_ab" for k, _ in OBS16] + ["alco", "pa_total", "cotinine", "bmi"],
         "9835a67923": [k if k != "carotene" else "carotene_ab" for k, _ in OBS16],
         "b966d3cfdc": [k for k in B966_CUTS if k not in ("pa", "bmi", "cot")] + ["alco", "pa_leisure", "cotinine", "bmi"]}

STATUS_NOTE = {
    "71eedfeafb": ("computed", "all 20 components named (Methods, Assessment of OBS); PA 1999-2006 from PAQIAF (MET 8/4 by level) + transport; 1999-2000 lacks carotenoids, scored on 19 of 20",
                   "R1: tertiles in age>=20, BMI>=25, non-pregnant; carotenoids = 5 carotenoids; alcohol by sex-specific tertiles | "
                   "R2: tertiles in all; carotenoids = beta-carotene; alcohol 2/1/0 categories"),
    "9e2e6087fa": ("not_reconstructable", "paper says 28 components but names 27; vitamin B6 not named (prereg section 4 applied strictly)",
                   "R1: not reconstructable | sensitivity row: 27 named parameters, day-1 recall; its R2: 28 (adds B6), mean of 2 recalls"),
    "988df86a79": ("computed", "28 parameters named, first-day recall stated, constants from Shivappa 2014",
                   "R1: day-1 recall | R2: mean of 2 recalls"),
    "fe82b9da6a": ("not_reconstructable", "the 16 dietary components are not named anywhere in the text, tables or XML (only cited to ref 20, Ke et al. 2023)",
                   "R1: not reconstructable | R2 (sensitivity): common 16-nutrient list, day 1, non-sex-specific tertiles in diabetes_any age>=20"),
    "aca6f864fb": ("computed", "all 20 components named with categories (section 2.3); PA 2005-2006 from PAQIAF with each activity's PADMETS",
                   "R1: mean of recalls; tertiles in age>=18; carotene alpha+beta | R2: day 1; tertiles in all; beta-carotene"),
    "9835a67923": ("computed", "16 components named in Table 1; scoring rule only in Supplementary Table 2 (not available locally); tertile 0-2 assumed",
                   "R1: non-sex-specific tertiles in diabetes_any non-pregnant 500-5000 kcal; carotene alpha+beta | R2: sex-specific; beta-carotene"),
    "b966d3cfdc": ("computed", "17 components named, sex-specific cut points printed in Table 5",
                   "R1: day-1 recall; PA leisure MET-min/week | R2: mean of recalls; PA all domains"),
}


SENS = {"9e2e6087fa": (1, 2), "fe82b9da6a": (2, None)}
SENS_NOTE = {"9e2e6087fa": "sensitivity: the 27 named parameters, day-1 recall; n_index_r2 = 28 parameters (adds B6), mean of 2 recalls",
             "fe82b9da6a": "sensitivity: common 16-nutrient OBS list (as aca6f864fb), day 1, non-sex-specific tertiles in diabetes_any age>=20"}


def unavailable(df, rec):
    d = diet(df, "day1"); L = lifestyle(df); allc = pd.concat([d, L], axis=1)
    out = []
    for c in CYC[rec]:
        m = df.cycle == c
        miss = [k for k in NEEDS[rec] if allc.loc[m, k].notna().sum() == 0]
        if miss:
            out.append(f"{c}: {','.join(miss)}")
    return "; ".join(out) if out else "none"


def main():
    full = ni.load_blind(nc.ALL_CYCLES[:10])
    rows = []
    for rec, f in DIET.items():
        df = full[full.cycle.isin(CYC[rec])].copy()
        stat, reason, rd = STATUS_NOTE[rec]
        res = {}
        for r in (1, 2):
            s = f(df, reading=r)
            a = s[(df.age >= 18)].dropna()
            res[r] = (len(a), a.mean() if len(a) else np.nan, a.std() if len(a) > 1 else np.nan)
            print(f"{rec} reading_{r}: n_adult={res[r][0]} mean={ni.trunc4(res[r][1])} sd={ni.trunc4(res[r][2])}")
        un = unavailable(df, rec)
        print(f"  unavailable: {un}")
        if stat == "computed":
            rows.append({"rec_id": rec, "status": stat, "reason": reason, "n_index": res[1][0],
                         "mean": ni.trunc4(res[1][1]), "sd": ni.trunc4(res[1][2]),
                         "n_index_r2": res[2][0], "mean_r2": ni.trunc4(res[2][1]), "sd_r2": ni.trunc4(res[2][2]),
                         "components_unavailable": un, "readings_decided": rd})
            continue
        rows.append({"rec_id": rec, "status": stat, "reason": reason, "n_index": 0, "mean": "NA", "sd": "NA",
                     "n_index_r2": "NA", "mean_r2": "NA", "sd_r2": "NA", "components_unavailable": un, "readings_decided": rd})
        sens = SENS[rec]   # (reading used as the sensitivity main, its reading_2 or None)
        r2 = res[sens[1]] if sens[1] else (np.nan, np.nan, np.nan)
        rows.append({"rec_id": rec, "status": "sensitivity_only", "reason": SENS_NOTE[rec], "n_index": res[sens[0]][0],
                     "mean": ni.trunc4(res[sens[0]][1]), "sd": ni.trunc4(res[sens[0]][2]),
                     "n_index_r2": r2[0] if sens[1] else "NA", "mean_r2": ni.trunc4(r2[1]), "sd_r2": ni.trunc4(r2[2]),
                     "components_unavailable": un, "readings_decided": rd})
    pd.DataFrame(rows).to_csv(os.path.join(OUT_DIR, "diet_recon_i139.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
