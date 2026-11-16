"""Build the blinded title/abstract screening packet from frame + Europe PMC raw JSON.
Usage: python -I tools/build_packet_i139.py <root>
Outputs screening_full/packet_screen.tsv (rec_id, title, abstract, year; sorted by rec_id),
screening_full/packet_screen_key.tsv (rec_id -> frame key, doi, pmid, matched source id, title check)."""
import csv, hashlib, html, json, os, re, sys, unicodedata
root = sys.argv[1]
sf = os.path.join(root, "screening_full")
rows = list(csv.DictReader(open(os.path.join(root, "frame", "frame_i139.tsv"), encoding="utf-8"), delimiter="\t"))


def clean(s):
    s = html.unescape(re.sub(r"<[^>]+>", " ", s or ""))
    return re.sub(r"\s+", " ", s).strip()


def norm_title(s):
    s = unicodedata.normalize("NFKD", clean(s)).lower()
    return re.sub(r"[^a-z0-9]", "", s)


def pick(res, r):
    pid, doi = r["pmid"].strip(), r["doi"].strip().lower()

    def match(x):
        return bool((pid and pid in (x.get("pmid"), x.get("pmcid"), x.get("id")))
                    or (doi and (x.get("doi") or "").lower() == doi))
    for cond in (lambda x: match(x) and x.get("abstractText"), match):
        for x in res:
            if cond(x):
                return x
    return None


packet, key, missing = [], [], []
for r in rows:
    rid = hashlib.sha256(r["key"].encode("utf-8")).hexdigest()[:10]
    d = json.load(open(os.path.join(sf, "raw", rid + ".json"), encoding="utf-8"))
    res = ((d.get("response") or {}).get("resultList") or {}).get("result", [])
    x = pick(res, r) or {}
    title = clean(r["title"])
    ab = clean(x.get("abstractText"))
    tcheck = "match" if x and norm_title(x.get("title")) == norm_title(r["title"]) else ("no_record" if not x else "differs")
    if not ab:
        missing.append(rid)
    packet.append((rid, title, ab, r["year"].strip()))
    key.append((rid, r["key"], r["doi"], r["pmid"], f'{x.get("source", "")}:{x.get("id", "")}' if x else "",
                d.get("query") or "", tcheck, "1" if ab else "0"))

packet.sort()
key.sort()
with open(os.path.join(sf, "packet_screen.tsv"), "w", encoding="utf-8", newline="") as f:
    f.write("rec_id\ttitle\tabstract\tyear\n")
    for p in packet:
        f.write("\t".join(c.replace("\t", " ") for c in p) + "\n")
with open(os.path.join(sf, "packet_screen_key.tsv"), "w", encoding="utf-8", newline="") as f:
    f.write("rec_id\tframe_key\tdoi\tpmid\tmatched_epmc_id\tquery_used\ttitle_check\thas_abstract\n")
    for k in key:
        f.write("\t".join(k) + "\n")
print("records", len(packet), "missing_abstract", len(missing), " ".join(missing))
print("title_check", {v: sum(1 for k in key if k[6] == v) for v in ("match", "differs", "no_record")})
