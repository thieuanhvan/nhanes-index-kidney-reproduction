"""Fetch abstracts by DOI from Europe PMC (core). python -I -u abs_by_doi.py out.jsonl doi..."""
import json, sys, time, urllib.parse, urllib.request, urllib.error, ctypes, re
try: ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x4000)
except Exception: pass
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
UA = "Mozilla/5.0 (compatible; literature-search-script/1.0)"
out = open(sys.argv[1], "w", encoding="utf-8")
for doi in sys.argv[2:]:
    d = {}
    for w in (0, 15, 30):
        time.sleep(w)
        try:
            u = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({"query": f'DOI:"{doi}"', "format": "json", "resultType": "core"})
            with urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": UA}), timeout=60) as r: d = json.loads(r.read()); break
        except urllib.error.HTTPError as e:
            if e.code != 429: break
        except Exception: pass
    rs = d.get("resultList", {}).get("result", [])
    r = rs[0] if rs else {}
    rec = {"doi": doi, "title": r.get("title"), "journal": (r.get("journalInfo") or {}).get("journal", {}).get("title"), "year": r.get("pubYear"),
           "authors": r.get("authorString"), "abstract": re.sub(r"<[^>]+>", " ", r.get("abstractText") or "")}
    out.write(json.dumps(rec, ensure_ascii=False) + "\n"); out.flush()
    print("##", doi, "|", rec["title"], "|", rec["journal"], rec["year"], "|", (rec["authors"] or "")[:80]); print(rec["abstract"][:1800]); print()
    time.sleep(2.5)
