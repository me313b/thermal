"""report_export.py - turn the app's report sections into downloadable
Word (.docx), PowerPoint (.pptx) and PDF files.

Everything here is pure Python (python-docx, python-pptx, fpdf2) so the
export runs wherever the app runs. Figures are rasterised from the
live Plotly objects with kaleido when available; if kaleido is missing
the documents are built text-only and a note is added instead.

The input contract matches render_report_tab:
    secs : list of (title, markdown, [figure_keys])
    figs : dict  figure_key -> plotly Figure
    meta : dict  version, date, spec, ok, T_gov, T_limit, kwh, mass,
                 Cmax, chil_el
"""
import io
import re

INK = "0F172A"
MUTED = "475569"
ACCENT = "6366F1"
GOOD = "16A34A"
BAD = "DC2626"


# ------------------------------------------------------------------ #
#  Shared markdown parsing                                            #
# ------------------------------------------------------------------ #
def _clean_inline(t: str) -> str:
    """Strip the light markdown/LaTeX the report bodies use down to
    plain text (bold markers are handled separately for docx runs)."""
    t = re.sub(r"<br\s*/?>", "\n", t)
    t = re.sub(r"</?su[bp]>", "", t)
    t = re.sub(r"\[([^\]]+)\]\([^)]*\)", r"\1", t)      # links -> text
    t = re.sub(r"\$([^$]*)\$", lambda m: m.group(1), t)  # $x$ -> x
    t = t.replace("\\;", " ").replace("\\,", " ")
    t = re.sub(r"\\[a-zA-Z]+", "", t)                    # \frac etc.
    t = t.replace("{", "").replace("}", "")
    return t


def _bold_segments(t: str):
    """Split 'a **b** c' -> [('a ',False),('b',True),(' c',False)]."""
    out, bold = [], False
    for part in t.split("**"):
        if part:
            out.append((part, bold))
        bold = not bold
    return out


def parse_md(md: str):
    """Yield blocks: ('h', level, text) | ('bullets', [items]) |
    ('table', header_row, rows) | ('p', text)."""
    lines = md.split("\n")
    i, blocks = 0, []
    while i < len(lines):
        ln = lines[i].rstrip()
        if not ln.strip():
            i += 1
            continue
        if ln.lstrip().startswith("|"):
            tbl = []
            while i < len(lines) and lines[i].lstrip().startswith("|"):
                row = [c.strip() for c in lines[i].strip().strip("|")
                       .split("|")]
                if not all(re.fullmatch(r":?-{2,}:?", c) for c in row):
                    tbl.append(row)
                i += 1
            if tbl:
                blocks.append(("table", tbl[0], tbl[1:]))
            continue
        m = re.match(r"^(#{2,6})\s+(.*)", ln.strip())
        if m:
            blocks.append(("h", len(m.group(1)), m.group(2)))
            i += 1
            continue
        if re.match(r"^\s*[-*]\s+", ln):
            items = []
            while i < len(lines) and re.match(r"^\s*[-*]\s+", lines[i]):
                items.append(re.sub(r"^\s*[-*]\s+", "", lines[i]).strip())
                i += 1
            blocks.append(("bullets", items))
            continue
        para = [ln.strip()]
        i += 1
        while i < len(lines) and lines[i].strip() and \
                not lines[i].lstrip().startswith(("|", "-", "*", "#")):
            para.append(lines[i].strip())
            i += 1
        blocks.append(("p", " ".join(para)))
    return blocks


# ------------------------------------------------------------------ #
#  Figures -> PNG                                                     #
# ------------------------------------------------------------------ #
def figs_to_png(figs: dict, keys=None) -> dict:
    """Rasterise the requested figures. Returns {} if kaleido is not
    available so callers can degrade gracefully."""
    out = {}
    want = keys if keys is not None else list(figs)
    for k in want:
        f = figs.get(k)
        if f is None:
            continue
        try:
            h = int(getattr(f.layout, "height", None) or 420)
            out[k] = f.to_image(format="png", width=1000, height=h,
                                scale=2)
        except Exception:
            return {}
    return out


# ------------------------------------------------------------------ #
#  Word                                                               #
# ------------------------------------------------------------------ #
def build_docx(secs, meta, pngs=None) -> bytes:
    from docx import Document
    from docx.shared import Pt, Inches, RGBColor
    from docx.enum.text import WD_ALIGN_PARAGRAPH

    pngs = pngs or {}
    doc = Document()
    st_ = doc.styles["Normal"]
    st_.font.name = "Calibri"
    st_.font.size = Pt(10.5)

    def _p(text="", size=None, bold=False, color=None, align=None,
           space_after=6):
        p = doc.add_paragraph()
        p.paragraph_format.space_after = Pt(space_after)
        if align == "c":
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        for seg, b in _bold_segments(text):
            r = p.add_run(_clean_inline(seg))
            r.bold = b or bold
            if size:
                r.font.size = Pt(size)
            if color:
                r.font.color.rgb = RGBColor.from_string(color)
        return p

    # ---- cover ----
    _p("", space_after=60)
    _p("Immersion Pack Lab", 30, bold=True, color=INK, align="c",
       space_after=2)
    _p("Design report", 20, color=ACCENT, align="c", space_after=18)
    _p(meta["spec"], 12, color=MUTED, align="c", space_after=4)
    _p(f"{meta['version']}  ·  {meta['date']}", 10, color=MUTED,
       align="c", space_after=18)
    verdict = ("WITHIN LIMIT" if meta["ok"] else "OVER LIMIT")
    _p(f"{verdict} - governing temperature {meta['T_gov']:.1f} °C "
       f"against a {meta['T_limit']:.0f} °C limit", 12, bold=True,
       color=(GOOD if meta["ok"] else BAD), align="c", space_after=24)
    t = doc.add_table(rows=2, cols=4)
    t.style = "Table Grid"
    heads = ["Energy", "Pack mass", "Max continuous", "Chiller"]
    vals = [f"{meta['kwh']:.1f} kWh", f"{meta['mass']:.0f} kg",
            f"{meta['Cmax']:.2f} C", f"{meta['chil_el']/1000:.2f} kW el"]
    for j, (h_, v_) in enumerate(zip(heads, vals)):
        c = t.rows[0].cells[j].paragraphs[0].add_run(h_)
        c.bold = True
        c.font.size = Pt(9)
        t.rows[1].cells[j].paragraphs[0].add_run(v_).font.size = Pt(11)
    doc.add_page_break()

    # ---- sections ----
    for si, (title, md, fkeys) in enumerate(secs, 1):
        doc.add_heading(f"{si} · {title}", level=1)
        for blk in parse_md(md):
            if blk[0] == "h":
                doc.add_heading(_clean_inline(blk[2]),
                                level=min(blk[1] - 2, 3) or 2)
            elif blk[0] == "bullets":
                for it in blk[1]:
                    p = doc.add_paragraph(style="List Bullet")
                    for seg, b in _bold_segments(it):
                        r = p.add_run(_clean_inline(seg))
                        r.bold = b
            elif blk[0] == "table":
                head, rows = blk[1], blk[2]
                tb = doc.add_table(rows=1 + len(rows), cols=len(head))
                tb.style = "Table Grid"
                for j, h_ in enumerate(head):
                    r = tb.rows[0].cells[j].paragraphs[0].add_run(
                        _clean_inline(h_.replace("**", "")))
                    r.bold = True
                    r.font.size = Pt(8.5)
                for ri, row in enumerate(rows, 1):
                    for j, cell in enumerate(row[:len(head)]):
                        for seg, b in _bold_segments(cell):
                            r = tb.rows[ri].cells[j].paragraphs[0] \
                                .add_run(_clean_inline(seg))
                            r.bold = b
                            r.font.size = Pt(8.5)
                _p("", space_after=4)
            else:
                _p(blk[1])
        for fk in fkeys:
            if fk in pngs:
                doc.add_picture(io.BytesIO(pngs[fk]), width=Inches(6.2))
                _p(f"Figure: {fk}", 8, color=MUTED, space_after=10)
    buf = io.BytesIO()
    doc.save(buf)
    return buf.getvalue()


# ------------------------------------------------------------------ #
#  PowerPoint                                                         #
# ------------------------------------------------------------------ #
def _bullets_from_md(md: str, cap=6, width=110):
    """Pull the most slide-worthy lines out of a section body."""
    out = []
    for blk in parse_md(md):
        if blk[0] == "bullets":
            out += [_clean_inline(b.replace("**", "")) for b in blk[1]]
        elif blk[0] == "p":
            txt = _clean_inline(blk[1].replace("**", ""))
            first = re.split(r"(?<=[.!?])\s", txt, 1)[0]
            if 25 < len(first):
                out.append(first)
    seen, res = set(), []
    for b in out:
        b = b.strip()
        if len(b) > width:
            b = b[:width - 1].rsplit(" ", 1)[0] + "…"
        if b and b not in seen:
            seen.add(b)
            res.append(b)
        if len(res) >= cap:
            break
    return res


def build_pptx(secs, meta, pngs=None) -> bytes:
    from pptx import Presentation
    from pptx.util import Inches, Pt
    from pptx.dml.color import RGBColor
    from pptx.enum.text import PP_ALIGN

    pngs = pngs or {}
    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)
    blank = prs.slide_layouts[6]

    def _box(s, x, y, w, h):
        tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w),
                                  Inches(h))
        tb.text_frame.word_wrap = True
        return tb.text_frame

    def _run(p, text, size, color=INK, bold=False):
        r = p.add_run()
        r.text = text
        r.font.size = Pt(size)
        r.font.bold = bold
        r.font.color.rgb = RGBColor.from_string(color)
        return r

    # ---- title slide ----
    s = prs.slides.add_slide(blank)
    tf = _box(s, 0.9, 2.1, 11.5, 2.6)
    _run(tf.paragraphs[0], "Immersion Pack Lab", 46, INK, True)
    p2 = tf.add_paragraph()
    _run(p2, "Design report", 28, ACCENT, False)
    p3 = tf.add_paragraph()
    p3.space_before = Pt(16)
    _run(p3, _clean_inline(meta["spec"]), 16, MUTED)
    p4 = tf.add_paragraph()
    _run(p4, f"{meta['version']}  ·  {meta['date']}", 12, MUTED)
    tf2 = _box(s, 0.9, 5.6, 11.5, 0.7)
    verdict = "WITHIN LIMIT" if meta["ok"] else "OVER LIMIT"
    _run(tf2.paragraphs[0],
         f"{verdict}  ·  {meta['T_gov']:.1f} °C vs "
         f"{meta['T_limit']:.0f} °C limit", 15,
         GOOD if meta["ok"] else BAD, True)

    # ---- at-a-glance ----
    s = prs.slides.add_slide(blank)
    tf = _box(s, 0.9, 0.55, 11.5, 0.8)
    _run(tf.paragraphs[0], "At a glance", 30, INK, True)
    tiles = [("Energy", f"{meta['kwh']:.1f} kWh"),
             ("Pack mass", f"{meta['mass']:.0f} kg"),
             ("Max continuous", f"{meta['Cmax']:.2f} C"),
             ("Governing T", f"{meta['T_gov']:.1f} °C"),
             ("Limit", f"{meta['T_limit']:.0f} °C"),
             ("Chiller", f"{meta['chil_el']/1000:.2f} kW el")]
    for i, (k, v) in enumerate(tiles):
        x = 0.9 + (i % 3) * 4.05
        y = 1.9 + (i // 3) * 2.3
        tf = _box(s, x, y, 3.7, 1.8)
        _run(tf.paragraphs[0], k.upper(), 12, MUTED, True)
        pv = tf.add_paragraph()
        pv.space_before = Pt(6)
        _run(pv, v, 30, INK, True)

    # ---- one slide per section ----
    for si, (title, md, fkeys) in enumerate(secs, 1):
        s = prs.slides.add_slide(blank)
        img = next((k for k in fkeys if k in pngs), None)
        tw = 6.2 if img else 11.5
        tf = _box(s, 0.9, 0.5, tw, 0.5)
        _run(tf.paragraphs[0], f"SECTION {si}", 11, ACCENT, True)
        tf = _box(s, 0.9, 0.95, tw, 1.0)
        _run(tf.paragraphs[0], _clean_inline(title), 27, INK, True)
        bl = _bullets_from_md(md)
        tf = _box(s, 0.9, 2.05, tw, 4.9)
        for j, b in enumerate(bl):
            p = tf.paragraphs[0] if j == 0 else tf.add_paragraph()
            p.space_after = Pt(10)
            _run(p, "·  ", 14, ACCENT, True)
            _run(p, b, 14, INK)
        if img:
            s.shapes.add_picture(io.BytesIO(pngs[img]), Inches(7.4),
                                 Inches(1.3), width=Inches(5.4))
        ftf = _box(s, 0.9, 7.02, 11.5, 0.35)
        _run(ftf.paragraphs[0],
             f"Immersion Pack Lab {meta['version']}  ·  "
             f"{_clean_inline(meta['spec'])}", 9, MUTED)

    # ---- closing ----
    s = prs.slides.add_slide(blank)
    tf = _box(s, 0.9, 2.8, 11.5, 1.6)
    _run(tf.paragraphs[0], "Everything here is live in the app", 28,
         INK, True)
    p2 = tf.add_paragraph()
    p2.space_before = Pt(10)
    _run(p2, "Design, Duty, Results, Cockpit, Zones, Improve, Ideas, "
             "Safety, Compare, Learn, Validate, System and Report "
             "tabs - every number regenerates when the design moves.",
         14, MUTED)
    buf = io.BytesIO()
    prs.save(buf)
    return buf.getvalue()


# ------------------------------------------------------------------ #
#  PDF                                                                #
# ------------------------------------------------------------------ #
_PDF_MAP = {"Δ": "d", "×": "x", "→": "->", "·": "-", "±": "+/-",
            "η": "eta", "μ": "mu", "Ω": "ohm", "ν": "nu", "β": "beta",
            "α": "alpha", "ρ": "rho", "λ": "lambda", "θ": "theta",
            "σ": "sigma", "√": "sqrt", "≈": "~", "≤": "<=", "≥": ">=",
            "…": "...", "’": "'", "‘": "'", "“": '"', "”": '"',
            "–": "-", "—": "-", "″": '"', "′": "'", "⁻": "-", "₃₂": "32"}


def _pdf_txt(t: str) -> str:
    t = _clean_inline(t).replace("**", "")
    for k, v in _PDF_MAP.items():
        t = t.replace(k, v)
    return t.encode("latin-1", "replace").decode("latin-1")


def build_pdf(secs, meta, pngs=None) -> bytes:
    from fpdf import FPDF

    pngs = pngs or {}
    pdf = FPDF(format="A4")
    pdf.set_auto_page_break(True, margin=16)
    pdf.set_margins(16, 14, 16)
    W = 178

    def _rgb(hexs):
        return tuple(int(hexs[i:i + 2], 16) for i in (0, 2, 4))

    pdf.add_page()
    pdf.ln(46)
    pdf.set_font("Helvetica", "B", 30)
    pdf.set_text_color(*_rgb(INK))
    pdf.cell(W, 12, "Immersion Pack Lab", align="C", ln=1)
    pdf.set_font("Helvetica", "", 19)
    pdf.set_text_color(*_rgb(ACCENT))
    pdf.cell(W, 11, "Design report", align="C", ln=1)
    pdf.ln(7)
    pdf.set_font("Helvetica", "", 11)
    pdf.set_text_color(*_rgb(MUTED))
    pdf.multi_cell(W, 6, _pdf_txt(meta["spec"]), align="C")
    pdf.cell(W, 6, _pdf_txt(f"{meta['version']}  -  {meta['date']}"),
             align="C", ln=1)
    pdf.ln(6)
    pdf.set_font("Helvetica", "B", 12)
    pdf.set_text_color(*_rgb(GOOD if meta["ok"] else BAD))
    verdict = "WITHIN LIMIT" if meta["ok"] else "OVER LIMIT"
    pdf.cell(W, 7, _pdf_txt(
        f"{verdict} - {meta['T_gov']:.1f} degC vs "
        f"{meta['T_limit']:.0f} degC limit"), align="C", ln=1)
    pdf.ln(8)
    pdf.set_text_color(*_rgb(INK))
    kv = [("Energy", f"{meta['kwh']:.1f} kWh"),
          ("Pack mass", f"{meta['mass']:.0f} kg"),
          ("Max continuous", f"{meta['Cmax']:.2f} C"),
          ("Chiller", f"{meta['chil_el']/1000:.2f} kW el")]
    cw = W / 4
    pdf.set_font("Helvetica", "B", 9)
    for k, _ in kv:
        pdf.cell(cw, 6, k, border=1, align="C")
    pdf.ln()
    pdf.set_font("Helvetica", "", 11)
    for _, v in kv:
        pdf.cell(cw, 8, v, border=1, align="C")
    pdf.ln()

    def _table(head, rows):
        n = len(head)
        lens = [max(len(_pdf_txt(head[j])),
                    *(len(_pdf_txt(r[j])) if j < len(r) else 0
                      for r in rows)) if rows else len(_pdf_txt(head[j]))
                for j in range(n)]
        tot = sum(lens) or 1
        ws = [max(W * l / tot, 16) for l in lens]
        sc = W / sum(ws)
        ws = [w * sc for w in ws]
        pdf.set_font("Helvetica", "B", 7)
        pdf.set_fill_color(238, 240, 244)
        for j in range(n):
            pdf.cell(ws[j], 6, _pdf_txt(head[j])[:60], border=1,
                     fill=True)
        pdf.ln()
        pdf.set_font("Helvetica", "", 7)
        for r in rows:
            if pdf.get_y() > 268:
                pdf.add_page()
            for j in range(n):
                cell = _pdf_txt(r[j]) if j < len(r) else ""
                pdf.cell(ws[j], 6, cell[:int(ws[j] / 1.28)], border=1)
            pdf.ln()
        pdf.ln(2)

    for si, (title, md, fkeys) in enumerate(secs, 1):
        pdf.add_page()
        pdf.set_font("Helvetica", "B", 16)
        pdf.set_text_color(*_rgb(INK))
        pdf.multi_cell(W, 8, _pdf_txt(f"{si} - {title}"))
        pdf.ln(2)
        for blk in parse_md(md):
            pdf.set_text_color(*_rgb(INK))
            if blk[0] == "h":
                pdf.ln(1)
                pdf.set_font("Helvetica", "B", 12)
                pdf.multi_cell(W, 6.5, _pdf_txt(blk[2]))
            elif blk[0] == "bullets":
                pdf.set_font("Helvetica", "", 9.5)
                for it in blk[1]:
                    pdf.set_x(20)
                    pdf.multi_cell(W - 6, 5.4, "- " + _pdf_txt(it))
            elif blk[0] == "table":
                _table(blk[1], blk[2])
            else:
                pdf.set_font("Helvetica", "", 9.5)
                pdf.multi_cell(W, 5.4, _pdf_txt(blk[1]))
                pdf.ln(1)
        for fk in fkeys:
            if fk in pngs:
                import tempfile
                import os
                fd, pth = tempfile.mkstemp(suffix=".png")
                with os.fdopen(fd, "wb") as f:
                    f.write(pngs[fk])
                if pdf.get_y() > 175:
                    pdf.add_page()
                pdf.image(pth, w=W)
                os.unlink(pth)
                pdf.set_font("Helvetica", "I", 8)
                pdf.set_text_color(*_rgb(MUTED))
                pdf.cell(W, 5, _pdf_txt(f"Figure: {fk}"), ln=1)
    out = pdf.output()
    return bytes(out)
