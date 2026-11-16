"""KD5 competing-work scan before submission (prereg v7 section 9): Europe PMC and Semantic Scholar on three lines,
(a) reproduction and survey-design checks of NHANES, (b) covariate sets of published papers against a DAG or for
over-adjustment with re-estimated effects, (c) computational reproduction on open data. Fixed queries, publication
years 2024 onward, top 50 per query per source. Writes search/kd5_presubmission_<date>.tsv.
Usage: python -I tools/kd5_presubmission_scan_i139.py
"""
import csv, datetime, json, os, sys, time, urllib.parse, urllib.request

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
UA = {"User-Agent": "meta-research-competing-work-scan/1.0"}
Q = {
    "a": ['NHANES AND (reproducibility OR reproduction OR reanalysis OR "re-analysis") AND (published OR articles OR studies)',
          'NHANES AND ("survey weights" OR "sampling weights" OR "complex survey design") AND (misuse OR errors OR reporting OR quality)'],
    "b": ['("directed acyclic graph" OR DAG OR "over-adjustment" OR overadjustment) AND ("published studies" OR "systematic review" OR meta-research) AND (reanalysis OR "re-estimate" OR "effect estimates")',
          'NHANES AND (confounding OR "covariate selection" OR "adjustment set") AND (meta-research OR "systematic review" OR audit)'],
    "c": ['("computational reproducibility" OR "analytic reproducibility") AND ("open data" OR "public data") AND (epidemiology OR health)',
          '(NHANES OR "public health survey") AND ("paper mill" OR formulaic OR "redundant publication" OR "data dredging")'],
}
S2Q = {"a": ["reproducibility of published NHANES analyses", "survey weight errors in published NHANES studies"],
       "b": ["over-adjustment directed acyclic graph published observational studies reanalysis",
             "confounder adjustment sets NHANES published papers audit"],
       "c": ["computational reproducibility epidemiology public data reanalysis", "formulaic NHANES papers"]}


def get(url):
    for k in range(4):
        try:
            return json.load(urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=40))
        except Exception as e:
            if "429" in str(e):
                time.sleep(6 * (k + 1)); continue
            return {"error": str(e)}
    return {"error": "429"}


rows = []
for line, qs in Q.items():
    for q in qs:
        url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&pageSize=50&resultType=lite&query="
               + urllib.parse.quote(q + " AND PUB_YEAR:[2024 TO 2026]"))
        r = get(url)
        for w in (r.get("resultList") or {}).get("result", []):
            rows.append(["europepmc", line, q, w.get("pubYear"), w.get("doi", ""), w.get("title", "")])
        if "error" in r: rows.append(["europepmc", line, q, "", "", "ERROR " + r["error"]])
        time.sleep(1)
for line, qs in S2Q.items():
    for q in qs:
        r = get("https://api.semanticscholar.org/graph/v1/paper/search?limit=50&year=2024-&fields=title,year,externalIds&query="
                + urllib.parse.quote(q))
        for w in r.get("data", []) or []:
            rows.append(["s2", line, q, w.get("year"), (w.get("externalIds") or {}).get("DOI", ""), w.get("title")])
        if "error" in r: rows.append(["s2", line, q, "", "", "ERROR " + r["error"]])
        time.sleep(4)
out = os.path.join(ROOT, "search", "kd5_presubmission_%s.tsv" % datetime.date.today().strftime("%Y%m%d"))
with open(out, "w", encoding="utf-8", newline="") as fh:
    csv.writer(fh, delimiter="\t").writerows([["source", "line", "query", "year", "doi", "title"]] + rows)
seen = set()
for r in rows:
    k = (r[4] or r[5]).lower()
    if k in seen: continue
    seen.add(k); print(r[1], r[0], r[3], r[4], "|", r[5])
print("rows", len(rows), "unique", len(seen), "->", out)
