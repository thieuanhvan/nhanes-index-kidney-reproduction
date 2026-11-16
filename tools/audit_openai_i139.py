"""I139 blind three-field audit by a non-Claude model (OpenAI API), replacing the dropped author audit (prereg step 10)
and meeting the "two model families" safeguard for the fields that decide E3a.
Usage:
  python -I tools/audit_openai_i139.py models            list models available to the key (no paper sent)
  python -I tools/audit_openai_i139.py run <model>       audit every item of extraction/openai_audit/packet_openai_audit_i139.json
The API key is read from the environment variable OPENAI_API_KEY, or from the Windows user registry (setx); it is never
printed or written. Only open-access full texts of published papers and the extracted fields are sent.
Outputs: extraction/openai_audit/raw/<rec_id>.json (full response), extraction/openai_audit/audit_openai.tsv
"""
import os, sys, json, time, csv, urllib.request, urllib.error
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "extraction", "openai_audit")
RAW = os.path.join(D, "raw")
os.makedirs(RAW, exist_ok=True)

INSTR = """You check data extracted from a published paper against the paper's full text. You are an auditor; you do not
know why the data were extracted. Definition of the paper's "main estimate": the most fully adjusted estimate stated in
the ABSTRACT, preferring the continuous index over categories. If the abstract has several outcomes: CKD first, then
albuminuria, then low eGFR; within the same type, the first one stated. The value is then found in the body (usually a
regression table). If NHANES is one of several data sources, only the NHANES analysis counts.

Check three fields independently from the paper:
1. model: is this the label of the model (as named in the paper) that produced the main estimate?
2. covariates: is this the covariate list of that model as stated in the paper (table footnote or Methods)? Judge on
   content, not wording (synonyms and abbreviations are fine); a missing or extra variable is an error. Ignore the
   exposure itself.
3. estimate: are the effect estimate and its 95% CI lower and upper bounds those of the main estimate?
Some fields may contain errors; most do not.

Answer ONLY with a JSON object:
{"model": {"verdict": "correct|wrong|cannot_tell", "correct_value": "...", "evidence": "short quote"},
 "covariates": {"verdict": "...", "missing": [...], "extra": [...], "evidence": "short quote"},
 "estimate": {"verdict": "...", "correct_value": "...", "evidence": "short quote"}}"""


def api_key():
    k = os.environ.get("OPENAI_API_KEY")
    if not k and sys.platform == "win32":
        import winreg
        try:
            with winreg.OpenKey(winreg.HKEY_CURRENT_USER, "Environment") as h:
                k = winreg.QueryValueEx(h, "OPENAI_API_KEY")[0]
        except OSError:
            k = None
    if not k:
        sys.exit("OPENAI_API_KEY not found (environment or user registry)")
    return k


def call(path, payload=None, timeout=300):
    req = urllib.request.Request("https://api.openai.com/v1/" + path, data=None if payload is None else json.dumps(payload).encode(),
                                 headers={"Authorization": "Bearer " + api_key(), "Content-Type": "application/json"},
                                 method="GET" if payload is None else "POST")
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def main():
    if sys.argv[1] == "models":
        ms = sorted(m["id"] for m in call("models")["data"])
        print("\n".join(m for m in ms if m.startswith(("gpt", "o1", "o3", "o4"))))
        return
    model = sys.argv[2]
    items = json.load(open(os.path.join(D, "packet_openai_audit_i139.json"), encoding="utf-8"))
    rows = []
    for it in items:
        rid = it["rec_id"]
        rp = os.path.join(RAW, rid + ".json")
        if not os.path.exists(rp):
            text = open(os.path.join(ROOT, "fulltext_txt", rid + ".txt"), encoding="utf-8", errors="replace").read()
            fields = {k: it[k] for k in ("index", "outcome", "model", "covariates", "estimate", "exposure_scale")}
            msg = [{"role": "system", "content": INSTR},
                   {"role": "user", "content": "EXTRACTED FIELDS:\n" + json.dumps(fields, ensure_ascii=False) + "\n\nFULL TEXT:\n" + text}]
            for attempt in range(4):
                try:
                    r = call("chat/completions", {"model": model, "messages": msg, "response_format": {"type": "json_object"}, "seed": 20261114})
                    break
                except urllib.error.HTTPError as e:
                    body = e.read().decode(errors="replace")[:300]
                    print(rid, "HTTP", e.code, body, flush=True)
                    if e.code in (429, 500, 502, 503):
                        time.sleep(20 * (attempt + 1)); continue
                    raise
            json.dump(r, open(rp, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
        r = json.load(open(rp, encoding="utf-8"))
        a = json.loads(r["choices"][0]["message"]["content"])
        for f in ("model", "covariates", "estimate"):
            v = a.get(f, {})
            cv = v.get("correct_value", "") if f != "covariates" else f"missing: {v.get('missing', [])}; extra: {v.get('extra', [])}"
            rows.append({"item": it["item"], "rec_id": rid, "field": f, "verdict": v.get("verdict", "na"), "correct_value": cv,
                         "evidence": str(v.get("evidence", ""))[:300], "model": r.get("model", model)})
        print(it["item"], rid, [a.get(f, {}).get("verdict") for f in ("model", "covariates", "estimate")], flush=True)
    with open(os.path.join(D, "audit_openai.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
