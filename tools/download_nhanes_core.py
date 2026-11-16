"""Download the I139 core NHANES XPT files (public CDC/NCHS) into data-external/nhanes-extra/<cycle>/.
Appends to the existing download_log.json (same record format as the 07/11 downloader).
Rules: one process, >= 1.5 s between requests, generic User-Agent (no email), stop on HTTP 403/429 (never bypass).
Files already present (in the target folder) are skipped. A file that CDC does not serve is logged with status "absent (...)".
Usage: python -I download_nhanes_core.py"""
import datetime, hashlib, json, os, sys, time, urllib.error, urllib.request

ROOT = "<workdir>/data-external/nhanes-extra"
UA = "Mozilla/5.0 (compatible; data-download-script/1.0)"
BASE = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{y}/DataFiles/{f}"

STD = ["DEMO", "ALB_CR", "BIOPRO", "GLU", "TRIGLY", "HDL", "TCHOL", "GHB", "BMX", "BPX", "DIQ", "BPQ", "SMQ", "MCQ", "CBC"]
# (folder, year folder on CDC, list of file stems)
PLAN = [
    ("1999-2000", "1999", ["DEMO", "LAB16", "LAB18", "LAB10AM", "LAB13AM", "LAB13", "LAB10", "BMX", "BPX", "DIQ", "BPQ", "SMQ", "MCQ", "LAB25"]),
    ("2001-2002", "2001", [s + "_B" for s in ["DEMO", "L16", "L40", "L10AM", "L13AM", "L13", "L10", "BMX", "BPX", "DIQ", "BPQ", "SMQ", "MCQ", "L25"]]),
    ("2003-2004", "2003", [s + "_C" for s in ["DEMO", "L16", "L40", "L10AM", "L13AM", "L13", "L10", "BMX", "BPX", "DIQ", "BPQ", "SMQ", "MCQ", "L25"]]),
    ("2005-2006", "2005", [s + "_D" for s in STD]),
    ("2007-2008", "2007", [s + "_E" for s in STD]),
    ("2009-2010", "2009", [s + "_F" for s in STD]),
    ("2011-2012", "2011", [s + "_G" for s in STD]),
    ("2013-2014", "2013", [s + "_H" for s in STD + ["INS"]]),
    ("2015-2016", "2015", [s + "_I" for s in STD + ["INS"]]),
    ("2017-2018", "2017", [s + "_J" for s in STD + ["INS", "BPXO"]]),
    ("2017-2020", "2017", ["P_" + s for s in [x for x in STD if x != "BPX"] + ["INS", "BPXO"]]),
    ("2021-2023", "2021", [s + "_L" for s in [x for x in STD if x != "BPX"] + ["INS", "BPXO"]]),
]
EXTRA = sys.argv[1:]  # optional "folder,year,STEM" triples for follow-up requests
for e in EXTRA:
    fo, yr, st = e.split(","); PLAN.append((fo, yr, [st]))

logp = os.path.join(ROOT, "download_log.json")
log = json.load(open(logp, encoding="utf-8"))
now = lambda: datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def save():
    tmp = logp + ".tmp"
    json.dump(log, open(tmp, "w", encoding="utf-8"), indent=1); os.replace(tmp, logp)

n_ok = 0
for folder, y, stems in PLAN:
    os.makedirs(os.path.join(ROOT, folder), exist_ok=True)
    for stem in stems:
        f = stem + ".xpt"; dest = os.path.join(ROOT, folder, f); url = BASE.format(y=y, f=f)
        if os.path.exists(dest):
            continue
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                data = r.read()
            if data[:200].lower().lstrip().startswith((b"<!doctype", b"<html")) or not data.startswith(b"HEADER RECORD"):
                st, h, n = "absent (CDC returned HTML, not XPT)", None, 0
            else:
                tmp = dest + ".part"; open(tmp, "wb").write(data); os.replace(tmp, dest)
                st, h, n = "ok", hashlib.sha256(data).hexdigest(), len(data); n_ok += 1
        except urllib.error.HTTPError as e:
            if e.code in (403, 429):
                log.append({"cycle": folder, "file": f, "url": url, "status": f"HTTP {e.code} - stopped", "sha256": None, "bytes": 0, "time": now()})
                save(); print("STOP", e.code, url, flush=True); sys.exit(2)
            st, h, n = f"absent (HTTP {e.code})", None, 0
        except Exception as e:
            st, h, n = f"error: {str(e)[:80]}", None, 0
        log.append({"cycle": folder, "file": f, "url": url, "status": st, "sha256": h, "bytes": n, "time": now()})
        save(); print(folder, f, st, n, flush=True)
        time.sleep(1.6)
print("downloaded", n_ok)
