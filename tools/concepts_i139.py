"""I139 prereg 5.4: map a covariate name, as written by a paper, to concept codes (equivalence table, code form).
Locked before any covariate label is assigned; its sha256 is logged in verdict 9.11.
A name can map to several codes (for example "hypertension and diabetes history").
"""
import re

RULES = [
    ("age", r"\bage\b|\bages\b"),
    ("sex", r"\bsex\b|gender"),
    ("race", r"\brace\b|ethnic|racial"),
    ("education", r"educat|schooling"),
    ("pir", r"poverty|\bpir\b|income"),
    ("smoking", r"smok|cotinine|tobacco|cigarette"),
    ("alcohol", r"alcohol|drink"),
    ("physical_activity", r"physical activit|vigorous|moderate activit|exercise|\bmet\b|\bmets?[- ]?min|activity status|work activity|leisure"),
    ("bmi_waist", r"\bbmi\b|body mass|waist|\bwc\b|whtr|obes"),
    ("lipid_drugs", r"statin|lipid[- ]lowering|lipoprotein[- ]lowering|cholesterol[- ]lowering|antihyperlipid"),
    ("diabetes_dx", r"diabet|\bdm\b|\bt2dm\b"),
    ("glycaemia_measured", r"glucose|hba1c|glycated|glycohemoglobin|glycosylated|\bfpg\b|insulin\b|homa"),
    ("hypertension_dx", r"hypertension|\bhtn\b"),
    ("bp_measured", r"blood pressure|\bsbp\b|\bdbp\b|systolic|diastolic"),
    ("diabetes_drugs", r"metformin|hypoglyc|antihyperglyc|antidiabet|sglt|glp-?1"),
    ("bp_drugs", r"antihypertensive|\bace\b|\barb\b|diuretic|\brasi\b|renin angiotensin"),
    ("lipids", r"\bhdl|\bldl|lipoprotein|cholesterol|triglycer|\btg\b|\btc\b|hyperlipid|dyslipid|lipid profile|blood lipid"),
    ("inflammation", r"\bcrp\b|c-reactive|white blood|\bwbc\b|leukocyte|neutrophil|lymphocyte|platelet"),
    ("uric_acid", r"uric|urate|\bua\b"),
    ("kidney_measure", r"egfr|glomerular|creatinine|cystatin|uacr|\bacr\b|albuminuria|urinary albumin|urea|\bbun\b|kidney function|renal function|\bs?cr\b|\bucr\b"),
    ("cvd", r"cardiovascular|\bcvds?\b|ascvd|\bchf\b|heart failure|coronary|heart|stroke|myocard|angina"),
    ("cycle", r"cycle|survey year|period"),
]


def concepts(name):
    s = str(name).lower()
    out = {c for c, rx in RULES if re.search(rx, s)}
    if "diabetes_dx" in out and re.search(r"glucose|hba1c|fasting", s) and not re.search(r"history|diagnos|self", s):
        out.add("diabetes_measured")
    return out


def concept_set(lst):
    if isinstance(lst, dict):  # field form {value, quote, loc}
        lst = lst.get("value") or []
    s = set()
    for x in lst or []:
        s |= concepts(x)
    return s
