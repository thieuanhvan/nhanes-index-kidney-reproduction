"""I139: convert fetched full texts (JATS XML / PDF / HTML) to plain text with tables kept as rows.
Usage: python -I tools/to_text_i139.py fulltext/ fulltext_txt/
Reads fulltext/fetch_log.tsv; writes <rec_id>.txt and text_log.tsv (sha256 of each text).
"""
import sys, os, csv, hashlib, re


def jats(path):
    from lxml import etree
    t = etree.parse(path, etree.XMLParser(recover=True, resolve_entities=False, no_network=True))
    out = []

    def txt(e):
        return re.sub(r"\s+", " ", "".join(e.itertext())).strip()

    root = t.getroot()
    for ab in root.iter("abstract"):
        out.append("## ABSTRACT\n" + txt(ab))
    body = root.find(".//body")
    if body is not None:
        for el in body.iter():
            if el.tag == "title":
                out.append("\n## " + txt(el))
            elif el.tag == "p" and el.getparent() is not None and el.getparent().tag not in ("table-wrap-foot", "fn", "caption"):
                out.append(txt(el))
    for tw in root.iter("table-wrap"):
        lab = tw.find("label")
        cap = tw.find("caption")
        out.append("\n## TABLE " + (txt(lab) if lab is not None else "") + " " + (txt(cap) if cap is not None else ""))
        for tr in tw.iter("tr"):
            out.append(" | ".join(txt(c) for c in tr if c.tag in ("td", "th")))
        foot = tw.find("table-wrap-foot")
        if foot is not None:
            out.append("TABLE FOOTNOTE: " + txt(foot))
    return "\n".join(out)


def pdf(path):
    import pymupdf
    d = pymupdf.open(path)
    return "\n".join(f"\n## PAGE {i + 1}\n" + p.get_text() for i, p in enumerate(d))


def html(path):
    from bs4 import BeautifulSoup
    s = BeautifulSoup(open(path, "rb").read(), "lxml")
    for x in s(["script", "style", "nav", "header", "footer"]):
        x.decompose()
    for tr in s.find_all("tr"):
        tr.replace_with(" | ".join(c.get_text(" ", strip=True) for c in tr.find_all(["td", "th"])) + "\n")
    return re.sub(r"\n\s*\n+", "\n\n", s.get_text("\n"))


def main(src, dst):
    os.makedirs(dst, exist_ok=True)
    rows = list(csv.DictReader(open(os.path.join(src, "fetch_log.tsv"), encoding="utf-8"), delimiter="\t"))
    lg = csv.writer(open(os.path.join(dst, "text_log.tsv"), "w", encoding="utf-8", newline=""), delimiter="\t", lineterminator="\n")
    lg.writerow(["rec_id", "src_file", "kind", "chars", "sha256", "note"])
    for r in rows:
        if not r["file"]:
            continue
        p = os.path.join(src, r["file"])
        ext = r["file"].rsplit(".", 1)[1]
        try:
            t = {"xml": jats, "pdf": pdf, "html": html}[ext](p)
            note = ""
        except Exception as e:
            t, note = "", f"error:{type(e).__name__}"
        b = t.encode("utf-8")
        open(os.path.join(dst, r["rec_id"] + ".txt"), "wb").write(b)
        lg.writerow([r["rec_id"], r["file"], ext, len(t), hashlib.sha256(b).hexdigest(), note])


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
