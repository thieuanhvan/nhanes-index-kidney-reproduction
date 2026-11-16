"""I139 blind human extraction form (second human coder, 15/11/2026).
Same 32 papers as the OpenAI audit (22 with a "missing DAG-A confounder" label + 10 random), shuffled with a fixed seed so
the coder cannot tell the two groups apart. The form shows only where to look (DOI, table, outcome, contrast); none of our
extracted values. The coder fills the most-adjusted model, its estimate and a 12-concept covariate checklist.
v2 (16/11/2026, round-3 review point M5): adds 8 of the 12 reporting-checklist items (2, 3, 4, 5, 6, 8, 9, 11; items 1, 7
and 10 were yes in 47/47 and item 12 not applicable in 40/47), with "not reported" kept apart from "not done" for items
3, 6 and 11, and an explicit not-applicable rule for item 11.
Key (rec_id per row) is written separately and hashed; the form itself carries only a row number.
Usage: python -I tools/build_human_audit_i139.py <out_dir> <stamp>
"""
import os, sys, csv, json, hashlib, random
import openpyxl
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.worksheet.datavalidation import DataValidation

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT, STAMP = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)

ids = [p["rec_id"] for p in json.load(open(os.path.join(HERE, "extraction/openai_audit/packet_openai_audit_i139.json"), encoding="utf-8"))]
spec = {r["rec_id"]: r for r in csv.DictReader(open(os.path.join(HERE, "extraction/spec_table_i139.tsv"), encoding="utf-8"), delimiter="\t")}
title = {r["doi"].lower(): r["title"] for r in csv.DictReader(open(os.path.join(HERE, "frame/frame_i139.tsv"), encoding="utf-8"), delimiter="\t")}
random.Random("20261115-i139-human-audit").shuffle(ids)

CONCEPTS = [("age", "Tuổi"), ("sex", "Giới"), ("race", "Chủng tộc"), ("education", "Học vấn"), ("pir", "Thu nhập (PIR)"),
            ("smoking", "Hút thuốc"), ("alcohol", "Rượu bia"), ("physical_activity", "Hoạt động thể lực"),
            ("bmi_waist", "BMI hoặc vòng eo"), ("diabetes_dx", "Đái tháo đường"), ("hypertension_dx", "Tăng huyết áp"),
            ("lipid_drugs", "Thuốc hạ mỡ máu")]

CHECK = [("c2", "BK2 Nêu TÊN biến trọng số", "yn"),
         ("c3", "BK3 Dùng trọng số của mẫu con nhỏ nhất", "ynr"),
         ("c4", "BK4 Chia trọng số khi gộp chu kỳ", "ynr"),
         ("c5", "BK5 Dùng strata và PSU", "ynr"),
         ("c6", "BK6 Ước lượng theo miền cho nhóm con", "ynra"),
         ("c8", "BK8 Nêu cách xử lý thiếu dữ liệu", "yn"),
         ("c9", "BK9 Hiệu chuẩn xét nghiệm giữa các chu kỳ", "ynra"),
         ("c11", "BK11 Độ tin cậy ước lượng ở nhóm nhỏ", "ynra")]

wb = openpyxl.Workbook()
h = wb.active
h.title = "Huong dan"
guide = [
    "PHIẾU SOÁT ĐỘC LẬP: MÔ HÌNH HIỆU CHỈNH TRONG 32 BÀI NHANES",
    "",
    "Mục đích: người thứ hai đọc độc lập để đối chiếu với bản trích của người thứ nhất. Phiếu không hiện số đã trích, nên cô cứ đọc bài và điền theo bài.",
    "Thời gian: khoảng 8–12 phút mỗi bài, tổng 4–6 giờ; có thể làm nhiều lần, mỗi lần vài bài.",
    "",
    "Mỗi dòng ở trang 'Phieu' là một bài. Cột B có đường dẫn DOI, bấm để mở toàn văn (tất cả đều mở miễn phí).",
    "Cột E cho biết bảng chứa kết quả chính; cột F, G cho biết kết cục và phép so sánh cần đọc.",
    "",
    "Cô điền:",
    "1. Cột H: tên mô hình hiệu chỉnh NHIỀU BIẾN NHẤT của bảng đó (ví dụ 'Model 3'), đúng như trong bài.",
    "2. Cột I, J, K: OR (hoặc HR/β nếu bài dùng) và khoảng tin cậy 95 % của phép so sánh ở cột G, trong mô hình đó.",
    "3. Cột L đến W: mô hình đó CÓ hiệu chỉnh biến này không. Danh sách biến thường ở chú thích dưới bảng hoặc đoạn 'Statistical analysis'.",
    "   - Chọn 'Có' hoặc 'Không' trong danh sách thả xuống.",
    "   - Chủng tộc: chọn '≥ 3 nhóm' nếu bài dùng biến chủng tộc NHANES nhiều nhóm; '2 nhóm' nếu chỉ da trắng/khác; 'Không' nếu không có.",
    "   - Thu nhập: tính cả 'poverty income ratio', 'PIR', 'family income'.",
    "   - BMI hoặc vòng eo: một trong hai là 'Có'.",
    "   - Đái tháo đường, tăng huyết áp: 'Có' nếu mô hình đưa chẩn đoán hoặc tình trạng bệnh vào (không tính đo đường huyết hay huyết áp thô).",
    "4. Cột X đến AE: cách bài báo cáo thiết kế khảo sát (bảng kiểm BK, đọc phần Methods). Chọn trong danh sách thả xuống:",
    "   - 'Có': bài nói rõ đã làm. 'Không nêu': bài không nói gì. 'Nói là không làm': bài nói rõ đã không làm (ví dụ 'unweighted').",
    "   - 'Không áp dụng': chỉ dùng ở BK6, BK9, BK11 theo quy tắc dưới.",
    "   - BK2: 'Có' chỉ khi bài nêu tên biến (WTMEC2YR, WTSAF2YR, WTDRD1…); chỉ ghi 'weighted' là 'Không nêu'.",
    "   - BK3: trọng số của mẫu con nhỏ nhất trong các biến bài dùng (ví dụ trọng số mẫu nhịn đói khi có đường huyết đói, TG, LDL; trọng số ăn uống khi dùng chỉ số khẩu phần). Ghi 'Có' cả khi bài chỉ nói 'fasting subsample weights' mà không nêu tên biến.",
    "   - BK4: ví dụ chia trọng số 2 năm cho số chu kỳ, dùng 3,2 năm cho 2017–3/2020.",
    "   - BK6: 'Không áp dụng' nếu bài phân tích toàn bộ người lớn, chỉ loại người thiếu dữ liệu; nếu bài chỉ lấy một nhóm (ví dụ chỉ người đái tháo đường) mà không nói 'subpopulation/domain' thì 'Không nêu'.",
    "   - BK9: 'Không áp dụng' nếu mọi chu kỳ từ 2007 trở đi; nếu dùng creatinine 1999–2000 hoặc 2005–2006 mà không nói hiệu chuẩn thì 'Không nêu'.",
    "   - BK11: 'Không áp dụng' nếu bài không báo ước lượng theo nhóm con (không có bảng phân nhóm); nếu có bảng phân nhóm thì 'Có' khi bài bàn sai số tương đối/ước lượng không ổn định, ngược lại 'Không nêu'.",
    "5. Cột AF: ghi chú nếu không chắc, hoặc nếu bài không có mô hình/bảng như mô tả.",
    "",
    "Không cần đọc toàn bài. Không cần tra thêm nguồn nào khác.",
    "Xong thì gửi lại tệp này (giữ nguyên tên, thêm _co-nhan-soat ở cuối).",
]
for i, t in enumerate(guide, 1):
    c = h.cell(row=i, column=1, value=t)
    c.alignment = Alignment(wrap_text=True, vertical="top")
    if i == 1:
        c.font = Font(bold=True, size=14)
h.column_dimensions["A"].width = 120

ws = wb.create_sheet("Phieu")
head = ["STT", "DOI (bấm để mở)", "Tên bài", "Chỉ số", "Bảng", "Kết cục", "So sánh cần đọc",
        "Mô hình hiệu chỉnh nhiều nhất", "OR", "CI thấp", "CI cao"] + [n for _, n in CONCEPTS] + [n for _, n, _ in CHECK] + ["Ghi chú"]
ws.append(head)
thin = Side(style="thin", color="999999")
for j in range(1, len(head) + 1):
    c = ws.cell(row=1, column=j)
    c.font = Font(bold=True)
    c.alignment = Alignment(wrap_text=True, vertical="center", horizontal="center")
    c.fill = PatternFill("solid", fgColor="DDEBF7" if j <= 7 else ("E2EFDA" if 11 + len(CONCEPTS) < j < len(head) else "FFF2CC"))
    c.border = Border(top=thin, bottom=thin, left=thin, right=thin)
key = []
for i, rid in enumerate(ids, 1):
    s = spec[rid]
    doi = s["doi"].strip()
    row = [i, f"https://doi.org/{doi}", title.get(doi.lower(), ""), s["index_name"], s["source_table"], s["outcome_main"],
           s["exposure_scale"]] + [""] * (4 + len(CONCEPTS) + len(CHECK) + 1)
    ws.append(row)
    ws.cell(row=i + 1, column=2).hyperlink = f"https://doi.org/{doi}"
    ws.cell(row=i + 1, column=2).font = Font(color="0563C1", underline="single")
    for j in range(1, len(head) + 1):
        c = ws.cell(row=i + 1, column=j)
        c.alignment = Alignment(wrap_text=True, vertical="top")
        c.border = Border(top=thin, bottom=thin, left=thin, right=thin)
    key.append({"row": i, "rec_id": rid})
n = len(ids) + 1
yn = DataValidation(type="list", formula1='"Có,Không"', allow_blank=True)
race = DataValidation(type="list", formula1='"≥ 3 nhóm,2 nhóm,Không"', allow_blank=True)
ws.add_data_validation(yn); ws.add_data_validation(race)
OPTS = {"yn": '"Có,Không nêu"', "ynr": '"Có,Không nêu,Nói là không làm"', "ynra": '"Có,Không nêu,Nói là không làm,Không áp dụng"'}
dvs = {k: DataValidation(type="list", formula1=v, allow_blank=True) for k, v in OPTS.items()}
for dv in dvs.values():
    ws.add_data_validation(dv)
first = 12
for k, (code, _) in enumerate(CONCEPTS):
    col = openpyxl.utils.get_column_letter(first + k)
    (race if code == "race" else yn).add(f"{col}2:{col}{n}")
for k, (code, _, kind) in enumerate(CHECK):
    col = openpyxl.utils.get_column_letter(first + len(CONCEPTS) + k)
    dvs[kind].add(f"{col}2:{col}{n}")
widths = [5, 22, 40, 16, 9, 18, 22, 14, 7, 7, 7] + [9] * len(CONCEPTS) + [12] * len(CHECK) + [30]
for j, wdt in enumerate(widths, 1):
    ws.column_dimensions[openpyxl.utils.get_column_letter(j)].width = wdt
ws.freeze_panes = "H2"
ws.row_dimensions[1].height = 75

base = f"thieuanhvan_diabetesxai_paper6_phieu-soat_vi_v2_{STAMP}"
xp = os.path.join(OUT, base + ".xlsx")
wb.save(xp)
kp = os.path.join(HERE, "extraction", "human_audit_key_i139.json")
json.dump({"form": os.path.basename(xp), "concepts": [c for c, _ in CONCEPTS], "checklist_items": [c for c, _, _ in CHECK], "rows": key}, open(kp, "w", encoding="utf-8"), indent=1)
hs = hashlib.sha256(open(kp, "rb").read()).hexdigest()
open(os.path.join(HERE, "extraction", "human_audit_key_i139.sha256"), "w").write(f"{hs}  human_audit_key_i139.json\n")
print(xp, len(ids), "rows; key", hs[:12])
