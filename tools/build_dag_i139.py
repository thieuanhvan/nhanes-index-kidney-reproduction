"""Build Appendix A DAGs (dagitty code + figures) for I139 from table 5.2 of prereg_i139_v6_vi.md.
The table is the single source: cells are parsed by code, nothing is typed by hand except the edge
encoding of each role and the node identifiers.
python -I build_dag_i139.py   (CPU only, below-normal priority)"""
import ctypes, hashlib, json, os, re, sys
try: ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x4000)
except Exception: pass
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(HERE, "prereg_i139_v6_vi.md")
OUT = os.path.join(HERE, "dag")
os.makedirs(OUT, exist_ok=True)
FAM = ["D1", "D2", "D3", "D4a", "D4b", "D5"]
FAM_VI = {"D1": "D1 TyG và dẫn xuất", "D2": "D2 viêm và huyết học", "D3": "D3 ăn uống và lối sống",
          "D4a": "D4a lipid thuần", "D4b": "D4b có glucose, insulin hoặc nhân trắc", "D5": "D5 khác"}
# node identifiers (dagitty) in table order; label = first cell of the table row
IDS = ["Age_Sex", "Race_ethnicity", "Education_PIR", "Smoking", "Alcohol", "Physical_activity", "BMI_waist",
       "Diagnosed_diabetes", "Glycaemia_HbA1c_insulin", "Diagnosed_hypertension", "Measured_BP",
       "Lipid_lowering_drugs", "Diabetes_drugs", "Antihypertensive_drugs", "HDL_LDL_TC", "CRP_WBC",
       "Uric_acid", "Kidney_markers", "CVD", "NHANES_cycle"]
SHORT_VI = ["Tuổi, giới", "Chủng tộc", "Học vấn, PIR", "Hút thuốc", "Rượu", "Thể lực", "BMI, vòng eo",
            "ĐTĐ đã chẩn đoán", "HbA1c, glucose,\ninsulin", "THA đã chẩn đoán", "Huyết áp đo",
            "Thuốc hạ lipid", "Thuốc ĐTĐ", "Thuốc hạ áp", "HDL, LDL, TC", "CRP, bạch cầu",
            "Acid uric", "eGFR, creatinin,\nUACR…", "Bệnh tim mạch", "Chu kỳ NHANES"]
X, Y = "Index", "Kidney_disease"

def parse_table(md):
    lines = md.splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith("### 5.2"))
    rows = []
    for l in lines[i:]:
        if l.startswith("### 5.3"): break
        if l.startswith("| ") and not l.startswith("| Nút") and not l.startswith("|---"):
            cells = [c.strip() for c in l.strip().strip("|").split("|")]
            rows.append(cells)
    assert len(rows) == len(IDS), (len(rows), len(IDS))
    return rows

def role(cell):
    """main role of a cell + B+ flag + note. 'K-C hoặc K-M' -> K-M (drawn as mediator); conditional
    alternatives in parentheses or after 'nếu' are kept as note only."""
    c = cell.strip()
    bplus = "B+" in c.split("(")[0]
    m = re.match(r"(K-[CMOX]|G|T|P)", c)
    r = m.group(1)
    if r == "K-C" and "hoặc K-M" in c: r = "K-M"
    if r == "K-C" and re.search(r"nếu là thành phần, T nếu không", c): r = "T"
    note = c[m.end():].strip(" ,")
    return r, bplus, note

def edges_for(nid, r, bplus, variant):
    """edge encoding of each role; variant A or B"""
    if r == "G" or (r == "T" and bplus and variant == "B"): return [(nid, X), (nid, Y)], True
    if r in ("T", "P"): return [(nid, Y)], False          # T outside the DAG: risk factor of outcome only
    if r == "K-C": return [(nid, X)], False               # component of the index
    if r == "K-M": return [(X, nid), (nid, Y)], False     # mediator / descendant of exposure
    if r == "K-O": return [(Y, nid)], False               # descendant of outcome
    if r == "K-X": return [(X, nid), (Y, nid)], False     # collider
    raise ValueError(r)

def backdoor_ok(G, Z):
    if any(z in nx.descendants(G, X) for z in Z): return False
    H = G.copy(); H.remove_edges_from(list(G.out_edges(X)))
    return nx.is_d_separator(H, {X}, {Y}, set(Z))

def layout(roles, variant):
    pos = {X: (0.0, 0.0), Y: (10.0, 0.0)}
    groups = {"adj": [], "tonly": [], "comp": [], "med": [], "out": [], "col": []}
    for nid, (r, bp, _) in roles.items():
        if r == "G" or (r == "T" and bp and variant == "B"): groups["adj"].append(nid)
        elif r in ("T", "P"): groups["tonly"].append(nid)
        elif r == "K-C": groups["comp"].append(nid)
        elif r == "K-M": groups["med"].append(nid)
        elif r == "K-O": groups["out"].append(nid)
        else: groups["col"].append(nid)
    def row(ids, y, x0, x1):
        n = len(ids)
        for k, nid in enumerate(ids): pos[nid] = (x0 + (x1 - x0) * (k + 0.5) / max(n, 1), y)
    adj = groups["adj"]
    row(adj[: (len(adj) + 1) // 2], 4.2, 0.0, 10.0); row(adj[(len(adj) + 1) // 2:], 3.0, 0.0, 10.0)
    row(groups["med"], -1.3, 2.0, 8.0)
    row(groups["comp"], -2.6, -1.5, 1.5)
    row(groups["out"] + groups["col"], -3.8, 4.0, 11.5)
    tonly = groups["tonly"]
    row(tonly[: (len(tonly) + 1) // 2], 1.6, 1.5, 11.5); row(tonly[(len(tonly) + 1) // 2:], -2.6, 2.5, 11.5) if not groups["comp"] else row(tonly[(len(tonly) + 1) // 2:], -2.6, 2.6, 11.5)
    return pos, groups

STYLE = {"adj": ("#1f4e79", "white"), "tonly": ("#d9d9d9", "black"), "comp": ("#f4b183", "black"),
         "med": ("#c5e0b4", "black"), "out": ("#ffe699", "black"), "col": ("#ffe699", "black"),
         "X": ("#7030a0", "white"), "Y": ("#c00000", "white")}

def draw(ax, G, pos, groups, title):
    gof = {n: g for g, ns in groups.items() for n in ns}; gof[X] = "X"; gof[Y] = "Y"
    lab = dict(zip(IDS, SHORT_VI)); lab[X] = "Chỉ số\n(phơi nhiễm)"; lab[Y] = "Bệnh thận\n(kết cục)"
    for u, v in G.edges():
        (x0, y0), (x1, y1) = pos[u], pos[v]
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color="#555555", lw=0.7, shrinkA=17, shrinkB=17, alpha=0.8))
    for n, (x, y) in pos.items():
        fc, tc = STYLE[gof[n]]
        ax.add_patch(FancyBboxPatch((x - 0.78, y - 0.32), 1.56, 0.64, boxstyle="round,pad=0.04", fc=fc, ec="#333333", lw=0.6, zorder=3))
        ax.text(x, y, lab[n], ha="center", va="center", fontsize=6.2, color=tc, zorder=4)
    ax.set_xlim(-2.3, 12.3); ax.set_ylim(-4.6, 5.0); ax.axis("off"); ax.set_title(title, fontsize=10)

def main():
    md = open(SRC, encoding="utf-8").read()
    sha_src = hashlib.sha256(open(SRC, "rb").read()).hexdigest()
    rows = parse_table(md)
    manifest = {"source": "prereg_i139_v6_vi.md", "source_sha256": sha_src, "families": {}}
    for j, fam in enumerate(FAM):
        roles = {IDS[k]: role(rows[k][1 + j]) for k in range(len(IDS))}
        fig, axes = plt.subplots(1, 2, figsize=(16, 6.2))
        manifest["families"][fam] = {}
        for ax, variant in zip(axes, ["A", "B"]):
            G = nx.DiGraph(); G.add_nodes_from([X, Y] + IDS); G.add_edge(X, Y)
            Z = []
            for nid, (r, bp, note) in roles.items():
                es, adj = edges_for(nid, r, bp, variant); G.add_edges_from(es)
                if adj: Z.append(nid)
            assert nx.is_directed_acyclic_graph(G)
            ok = backdoor_ok(G, Z)
            minimal = all(not backdoor_ok(G, [z for z in Z if z != w]) for w in Z)
            pos, groups = layout(roles, variant)
            # dagitty code
            L = [f"// I139 Phụ lục A, họ {FAM_VI[fam]}, DAG {variant}. Sinh bằng tools/build_dag_i139.py từ bảng 5.2 của",
                 f"// prereg_i139_v6_vi.md (sha256 {sha_src}). Bảng là nguồn duy nhất; nếu khác nhau thì bảng đúng.",
                 "// Mã hoá vai trò: G (và T có B+ ở DAG B) = gây nhiễu, vào bộ hiệu chỉnh; T/P ngoài bộ = chỉ yếu tố nguy cơ của kết cục;",
                 "// K-C = thành phần (-> chỉ số); K-M = trung gian/hậu duệ phơi nhiễm; K-O = hậu duệ kết cục; K-X = va chạm.",
                 "// DAG ở mức họ chỉ số, chỉ vẽ cạnh quyết định bộ hiệu chỉnh; cạnh giữa các hiệp biến không vẽ.",
                 "dag {"]
            L.append(f'{X} [exposure,pos="{pos[X][0]:.2f},{-pos[X][1]:.2f}"]')
            L.append(f'{Y} [outcome,pos="{pos[Y][0]:.2f},{-pos[Y][1]:.2f}"]')
            for nid in IDS:
                r, bp, note = roles[nid]
                attr = "adjusted," if nid in Z else ""
                L.append(f'{nid} [{attr}pos="{pos[nid][0]:.2f},{-pos[nid][1]:.2f}"]  // {r}{" B+" if bp else ""}{(" | " + note) if note else ""}')
            for u, v in sorted(G.edges()): L.append(f"{u} -> {v}")
            L.append("}")
            p = os.path.join(OUT, f"dag_i139_{fam}_{variant}.dagitty")
            open(p, "w", encoding="utf-8").write("\n".join(L) + "\n")
            manifest["families"][fam][variant] = {"adjustment_set": Z, "backdoor_valid": ok, "minimal": minimal,
                                                  "roles": {k: v[0] + (" B+" if v[1] else "") for k, v in roles.items()},
                                                  "dagitty_sha256": hashlib.sha256(open(p, "rb").read()).hexdigest()}
            draw(ax, G, pos, groups, f"DAG {variant}: bộ hiệu chỉnh {len(Z)} nút; tiêu chí cửa sau {'đạt' if ok else 'KHÔNG đạt'}")
        fig.suptitle(f"I139, Phụ lục A, họ {FAM_VI[fam]} (từ bảng 5.2, bản 6)", fontsize=12)
        handles = [plt.Rectangle((0, 0), 1, 1, fc=STYLE[k][0], ec="#333333") for k in ["adj", "tonly", "comp", "med", "out"]]
        fig.legend(handles, ["Vào bộ hiệu chỉnh", "Tranh cãi hoặc biến chính xác, ngoài bộ", "Thành phần chỉ số (K-C)",
                             "Trung gian/hậu duệ phơi nhiễm (K-M)", "Hậu duệ kết cục hoặc va chạm (K-O, K-X)"],
                   loc="lower center", ncol=5, fontsize=8, frameon=False)
        fig.tight_layout(rect=(0, 0.05, 1, 0.95))
        for ext in ("png", "svg"):
            fig.savefig(os.path.join(OUT, f"dag_i139_{fam}.{ext}"), dpi=200 if ext == "png" else None)
        plt.close(fig)
    json.dump(manifest, open(os.path.join(OUT, "dag_manifest_i139.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    for fam, d in manifest["families"].items():
        for v, r in d.items(): print(fam, v, len(r["adjustment_set"]), "valid", r["backdoor_valid"], "minimal", r["minimal"], r["adjustment_set"])

if __name__ == "__main__": main()
