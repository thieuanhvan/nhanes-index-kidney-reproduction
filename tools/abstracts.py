"""Fetch abstracts from Europe PMC by DOI. python -I abstracts.py out.jsonl doi..."""
import sys, json, time, urllib.parse, urllib.request
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
UA = "Mozilla/5.0 (compatible; literature-search-script/1.0)"
def get(u):
    for w in (0, 10, 30):
        time.sleep(w)
        try:
            with urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": UA}), timeout=60) as r: return json.loads(r.read())
        except Exception as e: err = e
    return {}
with open(sys.argv[1], "a", encoding="utf-8") as f:
    for doi in sys.argv[2:]:
        d = get("https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({"query": f'DOI:"{doi}"', "format": "json", "resultType": "core"}))
        rs = d.get("resultList", {}).get("result", [])
        r = rs[0] if rs else {}
        rec = {"doi": doi, "title": r.get("title"), "year": r.get("pubYear"), "journal": (r.get("journalInfo") or {}).get("journal", {}).get("title"),
               "authors": r.get("authorString"), "abstract": r.get("abstractText"), "isOA": r.get("isOpenAccess"), "pmcid": r.get("pmcid")}
        f.write(json.dumps(rec, ensure_ascii=False) + "\n")
        print("=== ", doi, "|", rec["year"], "|", rec["journal"], "|", (rec["authors"] or "")[:80]); print(rec["title"]); print((rec["abstract"] or "NO ABSTRACT")[:1800]); print()
        time.sleep(1.5)
