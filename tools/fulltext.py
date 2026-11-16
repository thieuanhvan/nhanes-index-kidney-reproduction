"""Europe PMC full text (OA) by DOI: resolve PMCID, fetch fullTextXML, save text, print sentences matching regex.
python -I -u fulltext.py <outdir> <regex> doi..."""
import sys, os, re, json, time, urllib.parse, urllib.request
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
UA = "Mozilla/5.0 (compatible; literature-search-script/1.0)"
def get(u, js=True):
    for w in (0, 10, 30):
        time.sleep(w)
        try:
            with urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": UA}), timeout=90) as r:
                b = r.read(); return json.loads(b) if js else b.decode("utf-8", "replace")
        except Exception: pass
    return {} if js else ""
out, RX = sys.argv[1], re.compile(sys.argv[2], re.I); os.makedirs(out, exist_ok=True)
for doi in sys.argv[3:]:
    d = get("https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({"query": f'DOI:"{doi}"', "format": "json", "resultType": "core"})); time.sleep(1.5)
    rs = d.get("resultList", {}).get("result", []); r = rs[0] if rs else {}
    pmc = r.get("pmcid"); fn = os.path.join(out, re.sub(r"[^\w.-]", "_", doi) + ".txt")
    if pmc and r.get("isOpenAccess") == "Y":
        x = get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML", js=False); time.sleep(1.5)
        txt = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", x)); kind = "fulltext"
    else:
        txt = re.sub(r"<[^>]+>", " ", r.get("abstractText") or ""); kind = "abstract_only"
    open(fn, "w", encoding="utf-8").write(txt)
    print(f"=== {doi} | {kind} | {len(txt)} chars | {(r.get('title') or '')[:100]}")
    sents = [s for s in re.split(r"(?<=[.!?])\s+", txt) if RX.search(s)]
    for s in sents[:14]: print("   -", s[:330])
