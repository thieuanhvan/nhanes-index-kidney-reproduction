"""English DAG figures and dagitty code for I139, built from table 5.2 of the registered protocol prereg_i139_v7_en.md.
The registered DAG files (dag/, built from prereg_i139_v6_vi.md and archived in dag_i139.zip on OSF) are not changed;
this script writes an English rendering to dag/en/ and checks that every node role equals the registered manifest.
python -I build_dag_i139_en.py   (CPU only, below-normal priority)"""
import hashlib, json, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_dag_i139 as base  # registered builder: edge encoding, layout and backdoor check are reused unchanged
import networkx as nx
import matplotlib.pyplot as plt

HERE = base.HERE
SRC = os.path.join(HERE, "osf", "upload", "prereg_i139_v7_en.md")
if not os.path.exists(SRC):  # repository layout: the registered protocol sits at the root
    SRC = os.path.join(HERE, "prereg_i139_v7_en.md")
OUT = os.path.join(HERE, "dag", "en")
os.makedirs(OUT, exist_ok=True)
FAM_EN = {"D1": "D1 TyG and derivatives", "D2": "D2 inflammation and haematology", "D3": "D3 diet and lifestyle",
          "D4a": "D4a pure lipid", "D4b": "D4b glucose, insulin or anthropometry", "D5": "D5 other"}
SHORT_EN = ["Age, sex", "Race/ethnicity", "Education, PIR", "Smoking", "Alcohol", "Physical activity", "BMI, waist",
            "Diagnosed diabetes", "HbA1c, glucose,\ninsulin", "Diagnosed\nhypertension", "Measured BP",
            "Lipid-lowering\ndrugs", "Diabetes drugs", "Antihypertensive\ndrugs", "HDL, LDL, TC", "CRP, WBC",
            "Uric acid", "eGFR, creatinine,\nUACR...", "Cardiovascular\ndisease", "NHANES cycle"]
X, Y, IDS = base.X, base.Y, base.IDS


def parse_table(md):
    lines = md.splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith("### 5.2"))
    rows = []
    for l in lines[i:]:
        if l.startswith("### 5.3"):
            break
        if l.startswith("| ") and not l.startswith("| Node") and not l.startswith("|---"):
            rows.append([c.strip() for c in l.strip().strip("|").split("|")])
    assert len(rows) == len(IDS), (len(rows), len(IDS))
    return rows


def role(cell):
    """Same rule as the registered builder, English wording: 'K-C or K-M' -> K-M; 'K-C if a component, T if not' -> T."""
    c = cell.strip()
    bplus = "B+" in c.split("(")[0]
    m = re.match(r"(K-[CMOX]|G|T|P)", c)
    r = m.group(1)
    if r == "K-C" and "or K-M" in c:
        r = "K-M"
    if r == "K-C" and re.search(r"if a component, T if not", c):
        r = "T"
    return r, bplus, c[m.end():].strip(" ,")


def draw(ax, G, pos, groups, title):
    gof = {n: g for g, ns in groups.items() for n in ns}
    gof[X] = "X"
    gof[Y] = "Y"
    lab = dict(zip(IDS, SHORT_EN))
    lab[X] = "Index\n(exposure)"
    lab[Y] = "Kidney disease\n(outcome)"
    for u, v in G.edges():
        (x0, y0), (x1, y1) = pos[u], pos[v]
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color="#555555", lw=0.7, shrinkA=17, shrinkB=17, alpha=0.8))
    for n, (x, y) in pos.items():
        fc, tc = base.STYLE[gof[n]]
        ax.add_patch(base.FancyBboxPatch((x - 0.78, y - 0.32), 1.56, 0.64, boxstyle="round,pad=0.04",
                                         fc=fc, ec="#333333", lw=0.6, zorder=3))
        ax.text(x, y, lab[n], ha="center", va="center", fontsize=6.2, color=tc, zorder=4)
    ax.set_xlim(-2.3, 12.3)
    ax.set_ylim(-4.6, 5.0)
    ax.axis("off")
    ax.set_title(title, fontsize=10)


def main():
    sha_src = hashlib.sha256(open(SRC, "rb").read()).hexdigest()
    rows = parse_table(open(SRC, encoding="utf-8").read())
    reg = json.load(open(os.path.join(HERE, "dag", "dag_manifest_i139.json"), encoding="utf-8"))
    manifest = {"source": "osf/upload/prereg_i139_v7_en.md", "source_sha256": sha_src,
                "translation_of": "dag/dag_manifest_i139.json (registered, built from prereg_i139_v6_vi.md)",
                "families": {}}
    mismatches = []
    for j, fam in enumerate(base.FAM):
        roles = {IDS[k]: role(rows[k][1 + j]) for k in range(len(IDS))}
        fig, axes = plt.subplots(1, 2, figsize=(16, 6.2))
        manifest["families"][fam] = {}
        for ax, variant in zip(axes, ["A", "B"]):
            G = nx.DiGraph()
            G.add_nodes_from([X, Y] + IDS)
            G.add_edge(X, Y)
            Z = []
            for nid, (r, bp, _) in roles.items():
                es, adj = base.edges_for(nid, r, bp, variant)
                G.add_edges_from(es)
                if adj:
                    Z.append(nid)
            assert nx.is_directed_acyclic_graph(G)
            ok = base.backdoor_ok(G, Z)
            minimal = all(not base.backdoor_ok(G, [z for z in Z if z != w]) for w in Z)
            pos, groups = base.layout(roles, variant)
            got = {k: v[0] + (" B+" if v[1] else "") for k, v in roles.items()}
            want = reg["families"][fam][variant]["roles"]
            mismatches += [(fam, variant, k, want[k], got[k]) for k in IDS if want[k] != got[k]]
            if Z != reg["families"][fam][variant]["adjustment_set"]:
                mismatches.append((fam, variant, "adjustment_set", reg["families"][fam][variant]["adjustment_set"], Z))
            L = [f"// I139 Appendix A, family {FAM_EN[fam]}, DAG {variant}. Built by tools/build_dag_i139_en.py from table 5.2 of",
                 f"// prereg_i139_v7_en.md (sha256 {sha_src}). English rendering of the registered DAG (dag/, OSF egjy3).",
                 "// Roles: G (and T with B+ in DAG B) = confounder, in the adjustment set; T/P outside the set = risk factor of the outcome only;",
                 "// K-C = component (-> index); K-M = mediator or descendant of the exposure; K-O = descendant of the outcome; K-X = collider.",
                 "// Family-level DAG: only edges that decide the adjustment set are drawn; edges between covariates are not drawn.",
                 "dag {",
                 f'{X} [exposure,pos="{pos[X][0]:.2f},{-pos[X][1]:.2f}"]',
                 f'{Y} [outcome,pos="{pos[Y][0]:.2f},{-pos[Y][1]:.2f}"]']
            for nid in IDS:
                r, bp, note = roles[nid]
                attr = "adjusted," if nid in Z else ""
                L.append(f'{nid} [{attr}pos="{pos[nid][0]:.2f},{-pos[nid][1]:.2f}"]  // {r}{" B+" if bp else ""}{(" | " + note) if note else ""}')
            L += [f"{u} -> {v}" for u, v in sorted(G.edges())] + ["}"]
            p = os.path.join(OUT, f"dag_i139_{fam}_{variant}_en.dagitty")
            open(p, "w", encoding="utf-8").write("\n".join(L) + "\n")
            manifest["families"][fam][variant] = {"adjustment_set": Z, "backdoor_valid": ok, "minimal": minimal, "roles": got,
                                                  "dagitty_sha256": hashlib.sha256(open(p, "rb").read()).hexdigest()}
            draw(ax, G, pos, groups, f"DAG {variant}: adjustment set of {len(Z)} nodes; backdoor criterion {'met' if ok else 'NOT met'}")
        fig.suptitle(f"I139, Appendix A, family {FAM_EN[fam]} (from table 5.2, registered protocol)", fontsize=12)
        handles = [plt.Rectangle((0, 0), 1, 1, fc=base.STYLE[k][0], ec="#333333") for k in ["adj", "tonly", "comp", "med", "out"]]
        fig.legend(handles, ["In the adjustment set", "Contested or precision variable, outside the set", "Index component (K-C)",
                             "Mediator or descendant of the exposure (K-M)", "Descendant of the outcome or collider (K-O, K-X)"],
                   loc="lower center", ncol=5, fontsize=8, frameon=False)
        fig.tight_layout(rect=(0, 0.05, 1, 0.95))
        for ext in ("png", "svg"):
            fig.savefig(os.path.join(OUT, f"dag_i139_{fam}_en.{ext}"), dpi=200 if ext == "png" else None)
        plt.close(fig)
    manifest["role_check_vs_registered"] = "identical" if not mismatches else [list(m) for m in mismatches]
    json.dump(manifest, open(os.path.join(OUT, "dag_manifest_i139_en.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("role check vs registered:", "identical" if not mismatches else mismatches)
    for fam, d in manifest["families"].items():
        for v, r in d.items():
            print(fam, v, len(r["adjustment_set"]), "valid", r["backdoor_valid"], "minimal", r["minimal"])


if __name__ == "__main__":
    main()
