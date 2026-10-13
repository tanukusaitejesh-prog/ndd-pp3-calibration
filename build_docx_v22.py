"""Build the Genetics in Medicine submission file (MANUSCRIPT_v22.docx) from MANUSCRIPT_v22.md.
Double-spaced Times New Roman 12 pt, page numbers, Word tables, figures embedded after the legends.
Run: python build_docx_v21.py"""
import os, re
from docx import Document
from docx.enum.text import WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "MANUSCRIPT_v22.md"), encoding="utf8").read()
src = src.split("\n", 1)[1]                                   # title handled separately
title = open(os.path.join(HERE, "MANUSCRIPT_v22.md"), encoding="utf8").readline().lstrip("# ").strip()
src = re.sub(r"^\*Draft \d+.*?\*\n", "", src, flags=re.M)        # drop the internal draft note
FIGS = {"Figure 1.": "fig1_design.png", "Figure 2.": "figures/fig2_dependence_v21.png", "Figure 3.": "figures/fig3_cost_v21.png"}

doc = Document()
st = doc.styles["Normal"]
st.font.name = "Times New Roman"; st.font.size = Pt(12)
st.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
pf = st.paragraph_format; pf.line_spacing_rule = WD_LINE_SPACING.DOUBLE; pf.space_after = Pt(0)
for s in doc.sections:
    s.left_margin = s.right_margin = s.top_margin = s.bottom_margin = Inches(1)
    fp = s.footer.paragraphs[0]; fp.alignment = 1
    r = fp.add_run()
    for tag, txt in (("begin", None), (None, "PAGE"), ("end", None)):
        if tag:
            e = OxmlElement("w:fldChar"); e.set(qn("w:fldCharType"), tag); r._r.append(e)
        else:
            e = OxmlElement("w:instrText"); e.set(qn("xml:space"), "preserve"); e.text = txt; r._r.append(e)

TOK = re.compile(r"(\*\*[^*]+\*\*|\*[^*]+\*|`[^`]+`)")


def runs(par, text):
    for part in TOK.split(text):
        if not part:
            continue
        if not TOK.fullmatch(part):
            par.add_run(part)
        elif part.startswith("**"):
            par.add_run(part[2:-2]).bold = True
        elif part.startswith("*"):
            par.add_run(part[1:-1]).italic = True
        elif part.startswith("`"):
            par.add_run(part[1:-1])
        else:
            par.add_run(part)


h = doc.add_paragraph(); h.alignment = 1
r = h.add_run(title); r.bold = True; r.font.size = Pt(14)

lines = src.split("\n")
i = 0
while i < len(lines):
    l = lines[i]
    if not l.strip() or l.strip() == "---":
        i += 1; continue
    if l.startswith("|"):                                     # table
        rows = []
        while i < len(lines) and lines[i].startswith("|"):
            cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
            if not all(re.fullmatch(r":?-+:?", c) for c in cells):
                rows.append(cells)
            i += 1
        ncol = max(len(r_) for r_ in rows)
        tb = doc.add_table(rows=len(rows), cols=ncol); tb.style = "Table Grid"
        for ri, row in enumerate(rows):
            for ci in range(ncol):
                cell = tb.cell(ri, ci); cell.text = ""
                p = cell.paragraphs[0]; p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
                runs(p, row[ci] if ci < len(row) else "")
                for rr in p.runs:
                    rr.font.size = Pt(9)
                    if ri == 0:
                        rr.bold = True
        doc.add_paragraph()
        continue
    m = re.match(r"^(#{2,3}) (.*)", l)
    if m:
        if m.group(2) in ("Abstract", "Introduction", "Tables", "Figure legends", "Data Availability", "References"):
            doc.add_page_break()
        hp = doc.add_paragraph(); rr = hp.add_run(m.group(2)); rr.bold = True
        if m.group(1) == "##":
            rr.font.size = Pt(13)
        i += 1; continue
    if re.match(r"^(\d+\.|-) ", l):
        p = doc.add_paragraph(style="List Number" if re.match(r"^\d+\. ", l) and "## References" not in src[:0] else "List Bullet")
        text = re.sub(r"^(\d+\.|-) ", "", l)
        if re.match(r"^\d+\. ", l) and i > 0 and any(lines[j].startswith("## References") for j in range(max(0, i - 40), i)):
            p.style = doc.styles["Normal"]; text = l               # keep reference numbers as written
        runs(p, text)
        i += 1; continue
    p = doc.add_paragraph(); runs(p, l.replace("\\*", "*"))
    for key, fn in FIGS.items():
        if l.startswith(f"**{key}"):
            path = os.path.join(HERE, fn)
            if os.path.exists(path):
                doc.add_picture(path, width=Inches(6.5))
    i += 1

out = os.path.join(HERE, "MANUSCRIPT_v22.docx")
doc.save(out)
print("written", out)
