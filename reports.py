"""
reports.py
----------
Export utilities for the Genome Sequence Comparison Tool.

Provides:
  - CSV export of the mutation table.
  - Multi-page PDF report (ReportLab) with metadata, metrics, alignment
    preview (and an embedded Plotly alignment-map chart when kaleido is
    available), and mutation summary / detail tables.
"""

from __future__ import annotations

import io
import csv
from datetime import datetime
from typing import Optional

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
    PageBreak,
    KeepTogether,
    Image as RLImage,
)
from reportlab.lib.enums import TA_LEFT, TA_CENTER

# Kaleido check — needed for Plotly static image export.
# We import it here so the module loads even when kaleido is absent;
# the PDF generator will skip chart embedding and add a fallback note.
try:
    import plotly.io as pio
    _KALEIDO_AVAILABLE = True
except Exception:
    _KALEIDO_AVAILABLE = False


# ---------------------------------------------------------------------------
# Colour palette (ReportLab HexColor wrappers)
# ---------------------------------------------------------------------------

COL_HEADER   = colors.HexColor("#1a252f")
COL_ODD      = colors.HexColor("#f2f3f4")
COL_EVEN     = colors.white
COL_ACCENT   = colors.HexColor("#2980b9")
COL_RED      = colors.HexColor("#e74c3c")
COL_GREEN    = colors.HexColor("#27ae60")
COL_ORANGE   = colors.HexColor("#e67e22")
COL_PURPLE   = colors.HexColor("#8e44ad")

MUT_COLORS = {
    "Substitution": COL_RED,
    "Insertion":    COL_ACCENT,
    "Deletion":     COL_ORANGE,
    "Duplication":  COL_PURPLE,
}


# ---------------------------------------------------------------------------
# CSV export
# ---------------------------------------------------------------------------

def export_mutations_csv(df: pd.DataFrame) -> bytes:
    """
    Serialise the mutation DataFrame to CSV bytes (UTF-8).

    Parameters
    ----------
    df : pd.DataFrame
        Output of :func:`mutation_detection.detect_mutations`.

    Returns
    -------
    bytes
        Raw CSV content ready to pass to ``st.download_button``.
    """
    buf = io.StringIO()
    df.to_csv(buf, index=False, quoting=csv.QUOTE_NONNUMERIC)
    return buf.getvalue().encode("utf-8")


# ---------------------------------------------------------------------------
# PDF report
# ---------------------------------------------------------------------------
# Kaleido chart-to-PNG helper
# ---------------------------------------------------------------------------

def _try_embed_chart(fig, width_cm: float = 15.0) -> Optional[RLImage]:
    """
    Attempt to render a Plotly figure to PNG bytes using kaleido and return
    a ReportLab Image flowable sized to *width_cm* centimetres.

    Returns ``None`` (silently) when kaleido is not installed, so callers
    can emit a fallback text note instead of crashing.

    Parameters
    ----------
    fig
        A ``plotly.graph_objects.Figure`` instance.
    width_cm : float
        Desired width of the embedded image in the PDF (centimetres).

    Returns
    -------
    RLImage | None
    """
    if not _KALEIDO_AVAILABLE or fig is None:
        return None
    try:
        png_bytes = pio.to_image(fig, format="png", width=900, height=400, scale=1.5)
        img_io = io.BytesIO(png_bytes)
        # Maintain aspect ratio: 900x400 px => width_cm x (width_cm * 400/900)
        w = width_cm * cm
        h = w * (400 / 900)
        return RLImage(img_io, width=w, height=h)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# PDF report generator
# ---------------------------------------------------------------------------

def generate_pdf_report(
    seq1_id: str,
    seq2_id: str,
    seq_type: str,
    alignment_mode: str,
    metrics: dict[str, object],
    alignment_preview: str,
    mutations_df: pd.DataFrame,
    mutation_summary: dict[str, int],
    score_params: Optional[dict[str, float]] = None,
    alignment_fig=None,        # Optional plotly Figure for the alignment map
) -> bytes:
    """
    Build a multi-page A4 PDF report using ReportLab Platypus.

    Parameters
    ----------
    seq1_id : str
        Identifier of the first sequence.
    seq2_id : str
        Identifier of the second sequence.
    seq_type : str
        'DNA', 'RNA', or 'Protein'.
    alignment_mode : str
        'global' or 'local'.
    metrics : dict
        Output of :func:`similarity.calculate_similarity_metrics`.
    alignment_preview : str
        Formatted alignment string (first 500 chars will be used).
    mutations_df : pd.DataFrame
        Output of :func:`mutation_detection.detect_mutations`.
    mutation_summary : dict[str, int]
        Output of :func:`mutation_detection.mutation_summary`.
    score_params : dict, optional
        Scoring parameters (match, mismatch, open gap, extend gap).

    Returns
    -------
    bytes
        PDF file content.
    """
    buf = io.BytesIO()

    # ---- Document layout -----------------------------------------------
    doc = BaseDocTemplate(
        buf,
        pagesize=A4,
        rightMargin=2 * cm,
        leftMargin=2 * cm,
        topMargin=2.5 * cm,
        bottomMargin=2 * cm,
    )

    frame = Frame(
        doc.leftMargin,
        doc.bottomMargin,
        doc.width,
        doc.height,
        id="normal",
    )

    def _header_footer(canvas, doc_obj):
        canvas.saveState()
        # Header bar
        canvas.setFillColor(COL_HEADER)
        canvas.rect(0, A4[1] - 1.5 * cm, A4[0], 1.5 * cm, fill=True, stroke=False)
        canvas.setFillColor(colors.white)
        canvas.setFont("Helvetica-Bold", 10)
        canvas.drawString(2 * cm, A4[1] - 1.0 * cm, "Genome Sequence Comparison Tool — Report")
        # Footer
        canvas.setFillColor(COL_HEADER)
        canvas.setFont("Helvetica", 8)
        canvas.drawRightString(A4[0] - 2 * cm, 0.75 * cm, f"Page {doc_obj.page}")
        canvas.drawString(2 * cm, 0.75 * cm, f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        canvas.restoreState()

    doc.addPageTemplates(
        [PageTemplate(id="main", frames=[frame], onPage=_header_footer)]
    )

    # ---- Styles --------------------------------------------------------
    styles = getSampleStyleSheet()
    title_style   = ParagraphStyle("Title2",   parent=styles["Title"],   fontSize=18, spaceAfter=6)
    h1_style      = ParagraphStyle("H1",       parent=styles["Heading1"], fontSize=13, textColor=COL_ACCENT, spaceAfter=4)
    h2_style      = ParagraphStyle("H2",       parent=styles["Heading2"], fontSize=11, textColor=COL_HEADER, spaceAfter=3)
    body_style    = ParagraphStyle("Body",     parent=styles["Normal"],   fontSize=9,  leading=13)
    mono_style    = ParagraphStyle("Mono",     parent=styles["Code"],     fontSize=7.5, leading=11, fontName="Courier")

    # ---- Story ---------------------------------------------------------
    story = []

    # Cover / Title
    story.append(Spacer(1, 1.5 * cm))
    story.append(Paragraph("Genome Sequence Comparison Report", title_style))
    story.append(Spacer(1, 0.3 * cm))

    meta_data = [
        ["Report Date",      datetime.now().strftime("%Y-%m-%d %H:%M:%S")],
        ["Sequence 1 ID",    seq1_id],
        ["Sequence 2 ID",    seq2_id],
        ["Sequence Type",    seq_type],
        ["Alignment Mode",   alignment_mode.capitalize()],
    ]
    if score_params:
        meta_data += [
            ["Match Score",    str(score_params.get("match_score",    2.0))],
            ["Mismatch Score", str(score_params.get("mismatch_score", -1.0))],
            ["Open Gap Score", str(score_params.get("open_gap_score", -2.0))],
            ["Extend Gap Score", str(score_params.get("extend_gap_score", -0.5))],
        ]

    meta_table = _make_simple_table(meta_data, col_widths=[5 * cm, 10 * cm])
    story.append(meta_table)
    story.append(Spacer(1, 0.5 * cm))

    # Page 1 — Similarity Metrics
    story.append(Paragraph("1. Similarity Metrics", h1_style))
    story.append(Spacer(1, 0.2 * cm))

    metric_labels = {
        "seq1_length":      "Sequence 1 Length",
        "seq2_length":      "Sequence 2 Length",
        "alignment_length": "Alignment Length",
        "alignment_score":  "Alignment Score",
        "match_count":      "Matches",
        "mismatch_count":   "Mismatches",
        "gap_count":        "Total Gap Positions",
        "gap_open_count":   "Gap Open Events",
        "percent_identity": "Percent Identity (%)",
        "gc_content_seq1":  "GC Content — Seq 1 (%)",
        "gc_content_seq2":  "GC Content — Seq 2 (%)",
    }

    metric_rows = []
    for key, label in metric_labels.items():
        val = metrics.get(key)
        if val is None:
            continue
        display = f"{val:.4f}" if isinstance(val, float) else str(val)
        metric_rows.append([label, display])

    metric_table = _make_simple_table(metric_rows, col_widths=[8 * cm, 7 * cm])
    story.append(metric_table)
    story.append(Spacer(1, 0.5 * cm))

    # Alignment Preview
    story.append(Paragraph("2. Alignment Preview (first 1 000 characters)", h1_style))
    story.append(Spacer(1, 0.2 * cm))
    preview_text = alignment_preview[:1000].replace("\n", "<br/>")
    story.append(Paragraph(preview_text, mono_style))
    story.append(Spacer(1, 0.3 * cm))

    # Alignment Map Chart (embedded via kaleido if available)
    story.append(Paragraph("3. Alignment Map Chart", h1_style))
    story.append(Spacer(1, 0.2 * cm))
    chart_img = _try_embed_chart(alignment_fig, width_cm=15.0)
    if chart_img is not None:
        story.append(chart_img)
        story.append(Spacer(1, 0.2 * cm))
    else:
        note = (
            "<i>Interactive alignment map chart not embedded — install "
            "<b>kaleido</b> (<code>pip install kaleido</code>) to enable "
            "static chart export in PDF reports.</i>"
        )
        story.append(Paragraph(note, body_style))
    story.append(Spacer(1, 0.3 * cm))

    story.append(PageBreak())

    # Page 2 — Mutation Summary
    story.append(Paragraph("4. Mutation Summary", h1_style))
    story.append(Spacer(1, 0.2 * cm))

    summary_rows = [["Mutation Type", "Count"]]
    for mtype, count in mutation_summary.items():
        summary_rows.append([mtype, str(count)])

    summary_table = Table(summary_rows, colWidths=[8 * cm, 4 * cm])
    summary_table.setStyle(_summary_table_style())
    story.append(summary_table)
    story.append(Spacer(1, 0.5 * cm))

    # Detailed Mutation Table
    if not mutations_df.empty:
        story.append(Paragraph("5. Detailed Mutation Table", h1_style))
        story.append(Spacer(1, 0.2 * cm))

        detail_rows = [["Pos", "Seq1", "Seq2", "Type", "Description"]]
        for _, row in mutations_df.iterrows():
            detail_rows.append([
                str(row["Position"]),
                str(row["Seq1_Base"]),
                str(row["Seq2_Base"]),
                str(row["Mutation_Type"]),
                str(row["Description"])[:80],
            ])

        det_table = Table(
            detail_rows,
            colWidths=[1.5 * cm, 1.5 * cm, 1.5 * cm, 3 * cm, 9.5 * cm],
            repeatRows=1,
        )
        det_table.setStyle(_detail_table_style(mutations_df))
        story.append(det_table)
    else:
        story.append(
            Paragraph(
                "No mutations detected — the sequences are identical within the aligned region.",
                body_style,
            )
        )

    doc.build(story)
    return buf.getvalue()


# ---------------------------------------------------------------------------
# ReportLab table helpers
# ---------------------------------------------------------------------------

def _make_simple_table(
    rows: list[list[str]],
    col_widths: list[float],
) -> Table:
    """Build a two-column key/value table."""
    table = Table(rows, colWidths=col_widths)
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND",  (0, 0),  (0, -1), COL_ODD),
                ("FONTNAME",    (0, 0),  (0, -1), "Helvetica-Bold"),
                ("FONTSIZE",    (0, 0),  (-1, -1), 9),
                ("GRID",        (0, 0),  (-1, -1), 0.5, colors.lightgrey),
                ("ROWBACKGROUNDS", (0, 0), (-1, -1), [COL_ODD, COL_EVEN]),
                ("VALIGN",      (0, 0),  (-1, -1), "MIDDLE"),
                ("TOPPADDING",  (0, 0),  (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )
    return table


def _summary_table_style() -> TableStyle:
    return TableStyle(
        [
            ("BACKGROUND",  (0, 0), (-1, 0), COL_HEADER),
            ("TEXTCOLOR",   (0, 0), (-1, 0), colors.white),
            ("FONTNAME",    (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE",    (0, 0), (-1, -1), 9),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [COL_ODD, COL_EVEN]),
            ("GRID",        (0, 0), (-1, -1), 0.5, colors.lightgrey),
            ("ALIGN",       (1, 0), (1, -1), "CENTER"),
            ("TOPPADDING",  (0, 0), (-1, -1), 4),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ]
    )


def _detail_table_style(df: pd.DataFrame) -> TableStyle:
    """Build a TableStyle with per-row mutation-type colouring."""
    base = [
        ("BACKGROUND",    (0, 0),  (-1, 0), COL_HEADER),
        ("TEXTCOLOR",     (0, 0),  (-1, 0), colors.white),
        ("FONTNAME",      (0, 0),  (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE",      (0, 0),  (-1, -1), 8),
        ("GRID",          (0, 0),  (-1, -1), 0.4, colors.lightgrey),
        ("ROWBACKGROUNDS",(0, 1),  (-1, -1), [COL_ODD, COL_EVEN]),
        ("TOPPADDING",    (0, 0),  (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0),  (-1, -1), 3),
        ("WORDWRAP",      (4, 1),  (4, -1), True),
    ]

    for i, (_, row) in enumerate(df.iterrows(), start=1):
        col = MUT_COLORS.get(row["Mutation_Type"])
        if col:
            base.append(("TEXTCOLOR", (3, i), (3, i), col))
            base.append(("FONTNAME",  (3, i), (3, i), "Helvetica-Bold"))

    return TableStyle(base)
