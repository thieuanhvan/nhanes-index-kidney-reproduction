"""Freeze the sampling frame for the full I139 study (no sampling, no screening, no data access).
Input: search/sys_i139/P3_epmc.json (Europe PMC, query P3, retrieved 14/11/2026).
Frame rule: record mentions NHANES / National Health and Nutrition Examination Survey in title or abstract; dedupe by DOI then title.
Stratum: index family from title keywords (first match wins, order fixed below); outcome is NOT assigned here (decided at eligibility).
Random order: within each stratum, sort by sha256(seed + ':' + key) -> screening order is fixed before anyone reads a record.
Output: frame/frame_i139.tsv + frame/frame_i139.tsv.sha256. python -I tools/build_frame_i139.py"""
import json, re, os, hashlib, csv
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SEED = "20261114-i139-v4"
FAM = [
    ("F1_tyg", r"triglyceride[- ]glucose|\btyg"),
    ("F2_inflam", r"neutrophil|lymphocyte|platelet|monocyte|systemic immune|inflammation[- ]response|\bsii\b|\bsiri\b|\bnlr\b|\bplr\b|\bmlr\b|\bnpar\b|\bhalp\b|c-reactive|crp|albumin ratio|myeloid|aggregate index of systemic"),
    ("F3_diet_life", r"diet|dietary|nutrition|eating|food|healthy eating|\bhei\b|\bdii\b|\bcdai\b|oxidative balance|life'?s (essential|simple)|\ble8\b|physical activity|sleep|lifestyle|depress"),
    ("F4_lipid_anthro", r"lipid|cholesterol|\bhdl\b|\bldl\b|atherogenic|remnant|\baip\b|mets-ir|visceral|adiposity|body roundness|weight-adjusted waist|waist|\bwwi\b|\bbri\b|\babsi\b|lipid accumulation|cardiometabolic index|insulin resistance|homa"),
]
def fam(title):
    t = title.lower()
    for name, pat in FAM:
        if re.search(pat, t): return name
    return "F5_other"
recs = json.load(open(os.path.join(ROOT, "search", "sys_i139", "P3_epmc.json"), encoding="utf-8"))
seen, rows = set(), []
for r in recs:
    title = (r.get("title") or "").strip(); txt = (title + " " + (r.get("abstract") or "")).lower()
    if not re.search(r"nhanes|national health and nutrition", txt): continue
    key = (r.get("doi") or "").lower() or re.sub(r"[^a-z0-9]", "", title.lower())[:120]
    if key in seen: continue
    seen.add(key)
    rows.append({"key": key, "doi": r.get("doi") or "", "pmid": r.get("id") or "", "year": r.get("year") or "", "family": fam(title),
                 "order_hash": hashlib.sha256(f"{SEED}:{key}".encode()).hexdigest(), "title": re.sub(r"\s+", " ", title)})
rows.sort(key=lambda x: (x["family"], x["order_hash"]))
for f in sorted({x["family"] for x in rows}):
    for i, x in enumerate([y for y in rows if y["family"] == f], 1): x["rank_in_stratum"] = i
os.makedirs(os.path.join(ROOT, "frame"), exist_ok=True)
out = os.path.join(ROOT, "frame", "frame_i139.tsv")
cols = ["family", "rank_in_stratum", "order_hash", "doi", "pmid", "year", "title", "key"]
with open(out, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, cols, delimiter="\t"); w.writeheader(); [w.writerow({c: x[c] for c in cols}) for x in rows]
h = hashlib.sha256(open(out, "rb").read()).hexdigest()
open(out + ".sha256", "w").write(f"{h} *frame/frame_i139.tsv\n")
from collections import Counter
print(len(rows), dict(sorted(Counter(x["family"] for x in rows).items())), h)
