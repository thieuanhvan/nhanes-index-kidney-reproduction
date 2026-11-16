"""I139: extra NHANES components not in nhanes_core (CRP, cotinine, copper/zinc, PSA, DXA trunk, bilirubin,
dietary recall totals, alcohol, physical activity, PHQ-9, sleep). Public CDC/NCHS XPT files in
data-external/nhanes-extra/<cycle>/ (sha256 in download_log.json). No outcome is read here.

add_extra(df) adds columns matched on cycle and SEQN:
  crp_mgL (LBXCRP mg/dL * 10, or LBXHSCRP mg/L), cotinine (LBXCOT ng/mL), copper (LBXSCU ug/dL), zinc (LBXSZN ug/dL),
  psa_total (LBXP1 ng/mL), psa_free (LBXP2 ng/mL), trunk_fat_g (DXXTRFAT), trunk_lean_g (DXXTRLI, lean excl. BMC),
  bilirubin_umol (LBDSTBSI umol/L, BIOPRO), and presence flags has_diet1, has_alq, has_paq, has_dpq, has_slq.
diet_columns(cycle) lists the day-1 nutrient total columns (DR1T* or DRXT*) of a cycle.
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nhanes_core as nc

SUF = {"1999-2000": "", "2001-2002": "_B", "2003-2004": "_C", "2005-2006": "_D", "2007-2008": "_E", "2009-2010": "_F",
       "2011-2012": "_G", "2013-2014": "_H", "2015-2016": "_I", "2017-2018": "_J", "2021-2023": "_L"}

STEMS = {
    "CRP": {"1999-2000": "LAB11", "2001-2002": "L11_B", "2003-2004": "L11_C", "2005-2006": "CRP_D", "2007-2008": "CRP_E",
            "2009-2010": "CRP_F", "2015-2016": "HSCRP_I", "2017-2018": "HSCRP_J", "2017-2020": "P_HSCRP"},
    "COT": {"1999-2000": "LAB06", "2001-2002": "L06_B", "2003-2004": "L06COT_C", "2005-2006": "COT_D", "2007-2008": "COTNAL_E",
            "2009-2010": "COTNAL_F", "2011-2012": "COTNAL_G", "2013-2014": "COT_H", "2015-2016": "COT_I", "2017-2018": "COT_J",
            "2017-2020": "P_COT"},
    "CUSEZN": {"2011-2012": "CUSEZN_G", "2013-2014": "CUSEZN_H", "2015-2016": "CUSEZN_I"},
    "PSA": {"2001-2002": "L11PSA_B", "2003-2004": "L11PSA_C", "2005-2006": "PSA_D", "2007-2008": "PSA_E", "2009-2010": "PSA_F"},
    "DXX": {"2011-2012": "DXX_G", "2013-2014": "DXX_H", "2015-2016": "DXX_I", "2017-2018": "DXX_J"},
    "DIET1": {"1999-2000": "DRXTOT", "2001-2002": "DRXTOT_B", "2003-2004": "DR1TOT_C", **{c: "DR1TOT" + s for c, s in SUF.items() if s > "_C"}, "2017-2020": "P_DR1TOT"},
    "ALQ": {**{c: "ALQ" + s for c, s in SUF.items()}, "2017-2020": "P_ALQ"},
    "PAQ": {**{c: "PAQ" + s for c, s in SUF.items()}, "2017-2020": "P_PAQ"},
    "DPQ": {**{c: "DPQ" + s for c, s in SUF.items() if s >= "_D"}, "2017-2020": "P_DPQ"},
    "SLQ": {**{c: "SLQ" + s for c, s in SUF.items() if s >= "_D"}, "2017-2020": "P_SLQ"},
    # individual leisure activities, one row per activity: read the file directly (read_stem keeps one row per SEQN)
    "PAQIAF": {"1999-2000": "PAQIAF", "2001-2002": "PAQIAF_B", "2003-2004": "PAQIAF_C", "2005-2006": "PAQIAF_D"},
}

_cache = {}


def read_stem(cycle, comp):
    stem = STEMS.get(comp, {}).get(cycle)
    if stem is None:
        return None
    p = os.path.join(nc.ROOT, cycle, stem + ".xpt")
    if not os.path.exists(p):
        return None
    if p not in _cache:
        x = pd.read_sas(p, format="xport")
        x.columns = [c.upper() for c in x.columns]
        num = x.select_dtypes("number").columns
        x[num] = x[num].mask(x[num].abs() < 1e-70, 0.0)
        x["SEQN"] = x["SEQN"].astype(np.int64)
        if comp == "DXX" and "_MULT_" in x.columns:   # multiple-imputation files: keep imputation 1
            x = x[x["_MULT_"] == 1]
        _cache[p] = x.drop_duplicates("SEQN").set_index("SEQN")
    return _cache[p]


def _pull(df, comp, col_options, scale=1.0):
    out = pd.Series(np.nan, index=df.index)
    for c in df.cycle.unique():
        x = read_stem(c, comp)
        if x is None:
            continue
        col = next((k for k in col_options if k in x.columns), None)
        if col is None:
            continue
        m = (df.cycle == c).to_numpy()
        out[m] = x[col].reindex(df.loc[m, "SEQN"]).to_numpy(float) * (scale(col) if callable(scale) else scale)
    return out


def _has(df, comp):
    out = pd.Series(0.0, index=df.index)
    for c in df.cycle.unique():
        x = read_stem(c, comp)
        if x is None:
            continue
        m = (df.cycle == c).to_numpy()
        out[m] = df.loc[m, "SEQN"].isin(x.index).to_numpy(float)
    return out


def add_extra(df):
    df = df.copy()
    df["crp_mgL"] = _pull(df, "CRP", ["LBXHSCRP", "LBXCRP"], scale=lambda col: 1.0 if col == "LBXHSCRP" else 10.0)
    df["cotinine"] = _pull(df, "COT", ["LBXCOT"])
    df["copper"] = _pull(df, "CUSEZN", ["LBXSCU"])
    df["zinc"] = _pull(df, "CUSEZN", ["LBXSZN"])
    df["psa_total"] = _pull(df, "PSA", ["LBXP1"])
    df["psa_free"] = _pull(df, "PSA", ["LBXP2"])
    df["trunk_fat_g"] = _pull(df, "DXX", ["DXXTRFAT"])
    df["trunk_lean_g"] = _pull(df, "DXX", ["DXXTRLI"])   # trunk lean excluding bone mineral content
    bil = pd.Series(np.nan, index=df.index)
    for c in df.cycle.unique():
        x = nc.read(c, "BIOPRO")
        if x is None or "LBDSTBSI" not in x.columns:
            continue
        m = (df.cycle == c).to_numpy()
        bil[m] = x.drop_duplicates("SEQN").set_index("SEQN")["LBDSTBSI"].reindex(df.loc[m, "SEQN"]).to_numpy(float)
    df["bilirubin_umol"] = bil
    for comp, name in (("DIET1", "has_diet1"), ("ALQ", "has_alq"), ("PAQ", "has_paq"), ("DPQ", "has_dpq"), ("SLQ", "has_slq")):
        df[name] = _has(df, comp)
    return df


def diet_columns(cycle):
    x = read_stem(cycle, "DIET1")
    return [] if x is None else [c for c in x.columns if c.startswith(("DR1T", "DRXT"))]
