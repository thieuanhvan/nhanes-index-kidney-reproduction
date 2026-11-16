"""I139: fetch open-access full texts for the sampling queue (prereg section 2.3, step 2).
Source order: Europe PMC fullTextXML (OA subset) -> OpenAlex OA locations (PDF/landing, plain GET only).
No login, no paywall or Cloudflare bypass; a 403/challenge page is logged as 'blocked'.
Usage: python -I -u tools/fetch_fulltext_i139.py screening_full/sampling_queue.tsv fulltext/
"""
import sys, os, csv, json, time, hashlib, urllib.parse, urllib.request
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
UA = "Mozilla/5.0 (compatible; literature-search-script/1.0)"
GAP = 1.6


def get(u, accept=None):
    h = {"User-Agent": UA}
    if accept:
        h["Accept"] = accept
    try:
        with urllib.request.urlopen(urllib.request.Request(u, headers=h), timeout=90) as r:
            return r.status, r.headers.get("Content-Type", ""), r.read()
    except urllib.error.HTTPError as e:
        return e.code, "", b""
    except Exception as e:
        return -1, type(e).__name__, b""


def epmc(doi, pmid):
    q = f'DOI:"{doi}"' if doi else f"EXT_ID:{pmid} AND SRC:MED"
    st, _, b = get("https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(
        {"query": q, "format": "json", "resultType": "core"}))
    time.sleep(GAP)
    try:
        rs = json.loads(b).get("resultList", {}).get("result", [])
    except Exception:
        rs = []
    return rs[0] if rs else {}


def openalex_locs(doi):
    st, _, b = get("https://api.openalex.org/works/doi:" + urllib.parse.quote(doi))
    time.sleep(GAP)
    try:
        w = json.loads(b)
    except Exception:
        return []
    locs = [w.get("best_oa_location") or {}] + (w.get("oa_locations") or [])
    out = []
    for l in locs:
        for k in ("pdf_url", "landing_page_url"):
            u = l.get(k)
            if u and u not in out:
                out.append(u)
    return out


def main(queue, outdir):
    os.makedirs(outdir, exist_ok=True)
    rows = list(csv.DictReader(open(queue, encoding="utf-8"), delimiter="\t"))
    logp = os.path.join(outdir, "fetch_log.tsv")
    done = set()
    if os.path.exists(logp):
        done = {r["rec_id"] for r in csv.DictReader(open(logp, encoding="utf-8"), delimiter="\t") if r["status"] == "ok"}
    new = not os.path.exists(logp)
    log = open(logp, "a", encoding="utf-8", newline="")
    w = csv.writer(log, delimiter="\t", lineterminator="\n")
    if new:
        w.writerow(["rec_id", "family", "queue_pos", "role", "doi", "status", "source", "file", "sha256", "bytes", "note"])
    for r in rows:
        if r["rec_id"] in done:
            continue
        doi, rid = r["doi"], r["rec_id"]
        status, source, fn, note = "no_oa", "", "", ""
        e = epmc(doi, r["pmid"])
        pmc = e.get("pmcid")
        if pmc and e.get("isOpenAccess") == "Y":
            st, ct, b = get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML")
            time.sleep(GAP)
            if st == 200 and b"<body" in b:
                fn, source, status = f"{rid}.xml", f"epmc:{pmc}", "ok"
            else:
                note = f"epmc_xml_{st}"
        elif pmc:
            note = f"pmc_not_oa:{pmc}"
        if status != "ok" and doi:
            for u in openalex_locs(doi):
                st, ct, b = get(u, accept="application/pdf,text/html;q=0.8")
                time.sleep(GAP)
                if st == 200 and b[:5] == b"%PDF-":
                    fn, source, status = f"{rid}.pdf", u, "ok"
                    break
                if st == 200 and b"<html" in b[:2000].lower() and len(b) > 60000 and b"cf-challenge" not in b and b"Just a moment" not in b:
                    fn, source, status = f"{rid}.html", u, "ok_html"
                    break
                if st in (401, 403) or b"Just a moment" in b[:5000]:
                    note += f" blocked:{urllib.parse.urlparse(u).netloc}"
            else:
                pass
        h = n = ""
        if fn:
            open(os.path.join(outdir, fn), "wb").write(b)
            h, n = hashlib.sha256(b).hexdigest(), len(b)
        w.writerow([rid, r["family"], r["queue_pos"], r["role"], doi, status, source, fn, h, n, note.strip()])
        log.flush()
        print(rid, r["family"], r["role"], status, source[:60], note.strip()[:80])


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
