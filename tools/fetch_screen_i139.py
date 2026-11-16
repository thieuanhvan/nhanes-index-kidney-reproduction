"""Fetch title + abstract for every frame record from Europe PMC REST (one process, >= 1.5 s between calls).
Usage: python -I tools/fetch_screen_i139.py <root>
Writes screening_full/raw/<rec_id>.json and appends screening_full/fetch_log.tsv. Stops on HTTP 403/429."""
import csv, datetime, hashlib, json, os, sys, time, urllib.error, urllib.parse, urllib.request
root = sys.argv[1]
UA = "Mozilla/5.0 (compatible; literature-screening-script/1.0)"
API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?"
out_dir = os.path.join(root, "screening_full", "raw")
log_path = os.path.join(root, "screening_full", "fetch_log.tsv")
rows = list(csv.DictReader(open(os.path.join(root, "frame", "frame_i139.tsv"), encoding="utf-8"), delimiter="\t"))
new_log = not os.path.exists(log_path)
log = open(log_path, "a", encoding="utf-8", newline="")
if new_log:
    log.write("rec_id\tquery\tutc_time\thttp_status\thit_count\n")
last = [0.0]

def call(q):
    wait = 1.5 - (time.time() - last[0])
    if wait > 0:
        time.sleep(wait)
    u = API + urllib.parse.urlencode({"query": q, "format": "json", "resultType": "core", "pageSize": 5})
    last[0] = time.time()
    try:
        with urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": UA}), timeout=60) as r:
            return r.status, json.loads(r.read())
    except urllib.error.HTTPError as e:
        if e.code in (403, 429):
            log.write(f"-\t{q}\t{datetime.datetime.utcnow().isoformat()}Z\t{e.code}\t-\n"); log.close()
            sys.exit(f"STOP: HTTP {e.code}")
        return e.code, {}

def queries(r):
    pid, doi = r["pmid"].strip(), r["doi"].strip()
    qs = []
    if pid.isdigit():
        qs.append(f"EXT_ID:{pid} AND SRC:MED")
    elif pid.upper().startswith("PMC"):
        qs.append(f"PMCID:{pid}")
    elif pid:
        qs.append(f"EXT_ID:{pid}")
    if doi:
        qs.append(f'DOI:"{doi}"')
    return qs

for r in rows:
    rid = hashlib.sha256(r["key"].encode("utf-8")).hexdigest()[:10]
    fp = os.path.join(out_dir, rid + ".json")
    if os.path.exists(fp):
        continue
    saved = None
    for q in queries(r):
        st, d = call(q)
        n = d.get("hitCount", 0) if d else 0
        log.write(f"{rid}\t{q}\t{datetime.datetime.utcnow().isoformat()}Z\t{st}\t{n}\n"); log.flush()
        res = (d.get("resultList") or {}).get("result", []) if d else []
        if res:
            saved = {"query": q, "response": d}
            if any(x.get("abstractText") for x in res):
                break
    with open(fp, "w", encoding="utf-8") as f:
        json.dump(saved or {"query": None, "response": None}, f, ensure_ascii=False)
    print(rid, "ok" if saved else "MISSING", flush=True)
log.close()
