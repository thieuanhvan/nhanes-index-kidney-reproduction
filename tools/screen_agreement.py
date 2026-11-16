"""Agreement between two title/abstract screening passes (I139, prereg v7 sections 2.2 and 3.2).
Usage: python -I tools/screen_agreement.py screener_A.tsv screener_B.tsv packet_screen.tsv [--out disagreements.tsv]
Reports: file validity, Cohen's kappa (3-class decision; include vs not include; exclude vs not exclude;
outcome-component flag; exposure-descendant label), the disagreement list, and a verbatim-quote check
against the packet title + abstract. Exit code 1 if either file is malformed."""
import argparse, csv, re, sys, unicodedata
from collections import Counter

DECISIONS = ("include", "exclude", "unclear")
REASONS = {
    "include": {"INC_ALL_MET"},
    "exclude": {"EX_PUB_TYPE", "EX_RETRACTED", "EX_NOT_NHANES", "EX_KIDNEY_EXPOSURE", "EX_NOT_INDEX",
                "EX_RESTRICTED_DATA", "EX_CKD_MORTALITY_PROGRESSION", "EX_OUTCOME_NOT_KIDNEY",
                "EX_NOT_CROSS_SECTIONAL", "EX_NO_ADJUSTED_ESTIMATE"},
    "unclear": {"UN_DATA_SOURCE", "UN_EXPOSURE", "UN_OUTCOME", "UN_DESIGN", "UN_ESTIMATE", "UN_NO_ABSTRACT"},
}
COLS = ["rec_id", "decision", "reason_code", "flag_outcome_component", "flag_exposure_descendant", "quote"]
DASHES = dict.fromkeys(map(ord, "‐‑‒–—―−"), "-")
QUOTES = {ord("‘"): "'", ord("’"): "'", ord("“"): '"', ord("”"): '"'}


def norm(s):
    s = unicodedata.normalize("NFKC", s or "").translate(DASHES).translate(QUOTES)
    return re.sub(r"\s+", " ", s).strip()


def read_tsv(path):
    with open(path, encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))


def kappa(a, b, labels):
    n = len(a)
    if n == 0:
        return float("nan")
    po = sum(x == y for x, y in zip(a, b)) / n
    ca, cb = Counter(a), Counter(b)
    pe = sum(ca[l] * cb[l] for l in labels) / (n * n)
    return float("nan") if pe == 1 else (po - pe) / (1 - pe)  # undefined when both raters use one label only


def load_screen(path, ids, tag):
    rows, problems = {}, []
    data = read_tsv(path)
    if data and set(COLS) - set(data[0].keys()):
        problems.append(f"{tag}: missing columns {sorted(set(COLS) - set(data[0].keys()))}")
    for i, r in enumerate(data, 2):
        rid = (r.get("rec_id") or "").strip()
        d = (r.get("decision") or "").strip().lower()
        rc = (r.get("reason_code") or "").strip()
        if rid in rows:
            problems.append(f"{tag} line {i}: duplicate rec_id {rid}")
        if rid not in ids:
            problems.append(f"{tag} line {i}: rec_id {rid!r} not in packet")
        if d not in DECISIONS:
            problems.append(f"{tag} line {i}: bad decision {d!r}")
        elif rc not in REASONS[d]:
            problems.append(f"{tag} line {i}: reason_code {rc!r} not allowed with decision {d}")
        for fl in ("flag_outcome_component", "flag_exposure_descendant"):
            if (r.get(fl) or "").strip() not in ("0", "1"):
                problems.append(f"{tag} line {i}: {fl} must be 0 or 1")
        rows[rid] = {"decision": d, "reason_code": rc, "foc": (r.get("flag_outcome_component") or "").strip(),
                     "fed": (r.get("flag_exposure_descendant") or "").strip(), "quote": r.get("quote") or ""}
    for rid in sorted(ids - rows.keys()):
        problems.append(f"{tag}: rec_id {rid} missing")
    return rows, problems


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("packet")
    ap.add_argument("--out", help="write the disagreement list as TSV")
    g = ap.parse_args()
    sys.stdout.reconfigure(encoding="utf-8")
    packet = {r["rec_id"]: r for r in read_tsv(g.packet)}
    ids = set(packet)
    A, pa = load_screen(g.a, ids, "A")
    B, pb = load_screen(g.b, ids, "B")
    print(f"packet records: {len(ids)}; A rows: {len(A)}; B rows: {len(B)}")
    for p in pa + pb:
        print("PROBLEM", p)
    common = sorted(ids & A.keys() & B.keys())
    da = [A[i]["decision"] for i in common]
    db = [B[i]["decision"] for i in common]
    print(f"\nrecords compared: {len(common)}")
    print("decision counts A:", dict(Counter(da)), " B:", dict(Counter(db)))
    print("cross-tab (A rows x B columns):")
    print("          " + "".join(f"{c:>10}" for c in DECISIONS))
    for x in DECISIONS:
        print(f"{x:>10}" + "".join(f"{sum(1 for p, q in zip(da, db) if p == x and q == y):>10}" for y in DECISIONS))
    agree = sum(x == y for x, y in zip(da, db))
    print(f"raw agreement 3-class: {agree}/{len(common)}")
    print(f"Cohen kappa 3-class (include/exclude/unclear): {kappa(da, db, DECISIONS):.4f}")
    ba = ["include" if x == "include" else "not" for x in da]
    bb = ["include" if x == "include" else "not" for x in db]
    print(f"Cohen kappa include vs not include: {kappa(ba, bb, ('include', 'not')):.4f}")
    ea = ["exclude" if x == "exclude" else "keep" for x in da]
    eb = ["exclude" if x == "exclude" else "keep" for x in db]
    print(f"Cohen kappa exclude vs not exclude (include + unclear go to full text): {kappa(ea, eb, ('exclude', 'keep')):.4f}")
    for fl, name in (("foc", "flag_outcome_component"), ("fed", "flag_exposure_descendant")):
        fa = [A[i][fl] for i in common]
        fb = [B[i][fl] for i in common]
        print(f"Cohen kappa {name}: {kappa(fa, fb, ('0', '1')):.4f} (A=1: {fa.count('1')}, B=1: {fb.count('1')})")

    dis = []
    for i in common:
        a, b = A[i], B[i]
        kinds = [k for k, f in (("decision", "decision"), ("reason", "reason_code"), ("flag_oc", "foc"), ("flag_ed", "fed"))
                 if a[f] != b[f]]
        if kinds:
            dis.append((i, ",".join(kinds), a["decision"], a["reason_code"], a["foc"], a["fed"],
                        b["decision"], b["reason_code"], b["foc"], b["fed"]))
    print(f"\ndisagreements: {len(dis)} records (decision differs: {sum(1 for d in dis if 'decision' in d[1])})")
    hdr = ("rec_id", "differs_in", "A_decision", "A_reason", "A_flag_oc", "A_flag_ed",
           "B_decision", "B_reason", "B_flag_oc", "B_flag_ed")
    print("\t".join(hdr))
    for d in dis:
        print("\t".join(d))
    if g.out:
        with open(g.out, "w", encoding="utf-8", newline="") as f:
            f.write("\t".join(hdr) + "\n")
            for d in dis:
                f.write("\t".join(d) + "\n")

    print("\nverbatim-quote check (whitespace, dashes and curly quotes normalised; case-sensitive):")
    for tag, S in (("A", A), ("B", B)):
        fails = []
        for i, r in S.items():
            if i not in packet:
                continue
            q = norm(r["quote"]).strip('"')
            src = norm(packet[i]["title"]) + " || " + norm(packet[i]["abstract"])
            if not q or q not in src:
                fails.append(i)
        print(f"{tag}: {len(S) - len(fails)}/{len(S)} quotes found verbatim; failures: {' '.join(sorted(fails)) or 'none'}")
    return 1 if (pa or pb) else 0


if __name__ == "__main__":
    sys.exit(main())
