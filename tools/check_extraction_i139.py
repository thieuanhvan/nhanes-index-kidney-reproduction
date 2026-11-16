"""I139 prereg 3.2 step 3 (text-level part): automatic checks on one extraction JSON against its full text.
Checks: verbatim quote present; estimate and covariates from the same table; OR inside CI; CI symmetric on log scale.
(N vs CDC file and XPT variable catalogue checks run later in the reconstruction stage.)
Usage: python -I tools/check_extraction_i139.py <fulltext.txt> <extraction.json> [more pairs...]
Prints one line per failed check and a summary line per paper; exit code 0 always.
"""
import sys, json, re, math
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

DASHES = "‐‑‒–—―−﹣－"


def norm(s):
    s = s or ""
    for d in DASHES:
        s = s.replace(d, "-")
    s = s.replace("­", "")
    for a, b in (("ﬁ", "fi"), ("ﬂ", "fl"), ("ﬀ", "ff"), ("ﬃ", "ffi"), ("ﬄ", "ffl"), ("&lt;", "<"), ("&gt;", ">"), ("&amp;", "&")):
        s = s.replace(a, b)
    s = s.replace(" ", " ").replace(" ", " ").replace(" ", " ")
    return re.sub(r"\s+", " ", s).strip().lower()


def fields(obj, path=""):
    if isinstance(obj, dict):
        if "quote" in obj and "value" in obj:
            yield path, obj
        for k, v in obj.items():
            if isinstance(v, (dict, list)) and not ("quote" in obj and k == "value"):
                yield from fields(v, f"{path}.{k}" if path else k)
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from fields(v, f"{path}[{i}]")


def table_of(loc):
    m = re.search(r"table\s*([0-9]+|[ivx]+|s[0-9]+)", (loc or "").lower())
    return m.group(1) if m else None


def num(f):
    try:
        return float(f["value"])
    except Exception:
        return None


def check(txt_path, js_path):
    T = norm(open(txt_path, encoding="utf-8", errors="replace").read())
    T_nospace = T.replace(" ", "")
    J = json.load(open(js_path, encoding="utf-8"))
    fails, n = [], 0
    for path, f in fields(J):
        if f.get("value") in (None, "", []):
            continue
        n += 1
        q = norm(f.get("quote"))
        if not q:
            fails.append(f"NOQUOTE {path}")
        elif q not in T and q.replace(" ", "") not in T_nospace:
            fails.append(f"VERBATIM {path}: {f.get('quote')[:90]!r}")
        v = f.get("value")
        if isinstance(v, (int, float)) and not isinstance(v, bool) and q:
            digits = re.sub(r"[^0-9.]", "", q.replace(",", ""))
            sv = ("%f" % abs(v)).rstrip("0").rstrip(".").lstrip("0") or "0"
            if sv not in digits:
                fails.append(f"VALUE_NOT_IN_QUOTE {path}: {v} not in quote")
    me = J.get("main_estimate", {})
    t_est = table_of(me.get("source_table", {}).get("loc") or me.get("estimate", {}).get("loc"))
    t_cov = table_of(J.get("covariates_footnote", {}).get("loc"))
    if J.get("covariates_footnote", {}).get("value") and t_est and t_cov and t_est != t_cov:
        fails.append(f"LOCATION estimate table {t_est} vs covariate footnote table {t_cov}")
    est, lo, hi = (num(me.get(k, {})) for k in ("estimate", "ci_low", "ci_high"))
    et = str(me.get("effect_type", {}).get("value") or "").upper()
    if None not in (est, lo, hi):
        if not (lo <= est <= hi):
            fails.append(f"NUMERIC estimate {est} outside CI ({lo}, {hi})")
        if et in ("OR", "PR", "HR", "RR") and min(est, lo, hi) > 0:
            mid = (math.log(lo) + math.log(hi)) / 2
            half = (math.log(hi) - math.log(lo)) / 2
            if half > 0 and abs(math.log(est) - mid) > max(0.02, 0.06 * half):
                fails.append(f"SYMMETRY log CI not centred on estimate: est {est}, CI ({lo}, {hi})")
    print(f"== {js_path}: {n} filled fields, {len(fails)} failed")
    for x in fails:
        print("   " + x)
    return n, len(fails)


if __name__ == "__main__":
    a = sys.argv[1:]
    tot_n = tot_f = 0
    for i in range(0, len(a), 2):
        n, f = check(a[i], a[i + 1])
        tot_n += n
        tot_f += f
    if tot_n:
        print(f"TOTAL {tot_f}/{tot_n} fields failed ({int(1000000 * tot_f / tot_n) / 10000}%)")
