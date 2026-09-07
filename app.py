"""
app.py
------
Main Streamlit dashboard for the Genome Sequence Comparison Tool.

Layout:
  Sidebar  — Mode, algorithm, scoring, and sample loaders.
  Tab 1    — Sequence Input & Validation.
  Tab 2    — Pairwise Alignment View.
  Tab 3    — Similarity Metrics & Mutation Table.
  Tab 4    — Visualizations & Interactive Plots.
  Tab 5    — Download & Report Export.
"""

from __future__ import annotations

import os
import pathlib
from typing import Optional, Literal

import pandas as pd
import streamlit as st

from typing import Optional, Literal, Any
from sequence_reader import read_sequence_from_upload, read_sequence_from_text, SequenceRecord
from alignment import run_pairwise_alignment, AlignmentResult
from similarity import calculate_similarity_metrics, metrics_to_display_rows
from mutation_detection import detect_mutations, mutation_summary
from visualization import (
    plot_alignment_map,
    plot_dot_plot,
    plot_base_composition,
    plot_mutation_distribution,
)
from reports import export_mutations_csv, generate_pdf_report


# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="Genome Sequence Comparison Tool",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


# ---------------------------------------------------------------------------
# Dataset directory
# ---------------------------------------------------------------------------

DATASETS_DIR = pathlib.Path(__file__).parent / "datasets"


def _load_sample_fasta(filename: str) -> bytes:
    """Read a sample FASTA file from the datasets directory."""
    fpath = DATASETS_DIR / filename
    if fpath.exists():
        return fpath.read_bytes()
    return b""


SAMPLE_FILES: dict[str, dict[str, str]] = {
    "DNA": {
        "Sample 1": "dna_sample_1.fasta",
        "Sample 2": "dna_sample_2.fasta",
    },
    "RNA": {
        "Sample 1": "rna_sample_1.fasta",
        "Sample 2": "rna_sample_2.fasta",
    },
    "Protein": {
        "Sample 1": "protein_sample_1.fasta",
        "Sample 2": "protein_sample_2.fasta",
    },
}


# ---------------------------------------------------------------------------
# Session-state initialisation
# ---------------------------------------------------------------------------

def _init_state() -> None:
    defaults: dict[str, object] = {
        "seq_type":         "DNA",
        "align_mode":       "global",
        "match_score":      2.0,
        "mismatch_score":   -1.0,
        "open_gap_score":   -2.0,
        "extend_gap_score": -0.5,
        "record1":          None,
        "record2":          None,
        "alignment":        None,
        "metrics":          None,
        "mutations_df":     None,
        "mut_summary":      None,
        "run_done":         False,
        "uploader_key":     0,   # Used to force-clear file uploaders
        "text1":            "",  # Tracks text area 1
        "text2":            "",  # Tracks text area 2
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------

def _render_sidebar() -> None:
    st.sidebar.image(
        "https://upload.wikimedia.org/wikipedia/commons/thumb/4/47/PNG_transparency_demonstration_1.png/280px-PNG_transparency_demonstration_1.png",
        width=40,
        use_container_width=False,
    )
    st.sidebar.title("🧬 Configuration")
    st.sidebar.markdown("---")

    # Sequence type
    st.sidebar.subheader("Sequence Type")
    seq_type = st.sidebar.radio(
        "Select sequence type:",
        options=["DNA", "RNA", "Protein"],
        index=["DNA", "RNA", "Protein"].index(st.session_state.seq_type),
        horizontal=True,
        key="seq_type_radio",
        label_visibility="collapsed",
    )
    st.session_state.seq_type = seq_type

    st.sidebar.markdown("---")

    # Alignment mode
    st.sidebar.subheader("Alignment Algorithm")
    align_mode = st.sidebar.selectbox(
        "Mode:",
        options=["global", "local"],
        index=0 if st.session_state.align_mode == "global" else 1,
        format_func=lambda x: (
            "🌍 Global (Needleman-Wunsch)" if x == "global"
            else "🔎 Local (Smith-Waterman)"
        ),
        key="align_mode_sel",
    )
    st.session_state.align_mode = align_mode

    st.sidebar.markdown("---")

    # Scoring parameters
    st.sidebar.subheader("Scoring Parameters")

    if seq_type == "Protein":
        st.sidebar.info(
            "**BLOSUM62** substitution matrix is used automatically for protein "
            "alignment — match/mismatch sliders are disabled.  "
            "Gap penalties still apply."
        )
    else:
        st.session_state.match_score = st.sidebar.number_input(
            "Match score", value=float(st.session_state.match_score), step=0.5, key="ms"
        )
        st.session_state.mismatch_score = st.sidebar.number_input(
            "Mismatch penalty", value=float(st.session_state.mismatch_score), step=0.5, key="mm"
        )

    st.session_state.open_gap_score = st.sidebar.number_input(
        "Gap open penalty", value=float(st.session_state.open_gap_score), step=0.5, key="go"
    )
    st.session_state.extend_gap_score = st.sidebar.number_input(
        "Gap extend penalty", value=float(st.session_state.extend_gap_score), step=0.25, key="ge"
    )

    st.sidebar.markdown("---")

    # Sample loaders
    st.sidebar.subheader("📂 Load Sample Sequences")
    col_s1, col_s2 = st.sidebar.columns(2)
    with col_s1:
        if st.button("Load Seq 1", key="load_s1", use_container_width=True):
            fname = SAMPLE_FILES[seq_type]["Sample 1"]
            fbytes = _load_sample_fasta(fname)
            if fbytes:
                rec = read_sequence_from_upload(fbytes, seq_type)
                st.session_state.record1 = rec
                st.session_state.run_done = False
                st.session_state.uploader_key += 1  # Clear file uploaders
                st.session_state.text1 = ""         # Clear text area 1
                st.toast(f"Loaded {fname} ✓", icon="📄")
            else:
                st.sidebar.error(f"Sample file '{fname}' not found.")
    with col_s2:
        if st.button("Load Seq 2", key="load_s2", use_container_width=True):
            fname = SAMPLE_FILES[seq_type]["Sample 2"]
            fbytes = _load_sample_fasta(fname)
            if fbytes:
                rec = read_sequence_from_upload(fbytes, seq_type)
                st.session_state.record2 = rec
                st.session_state.run_done = False
                st.session_state.uploader_key += 1  # Clear file uploaders
                st.session_state.text2 = ""         # Clear text area 2
                st.toast(f"Loaded {fname} ✓", icon="📄")
            else:
                st.sidebar.error(f"Sample file '{fname}' not found.")

    if st.sidebar.button("🔄 Reset All", key="reset_all", use_container_width=True):
        for key in [
            "record1", "record2", "alignment", "metrics",
            "mutations_df", "mut_summary", "run_done",
        ]:
            st.session_state[key] = None if key != "run_done" else False
        
        # Clear widget states
        st.session_state.uploader_key += 1
        st.session_state.text1 = ""
        st.session_state.text2 = ""
        st.rerun()

    st.sidebar.markdown("---")
    st.sidebar.caption(
        "Genome Sequence Comparison Tool  \n"
        "Built with Streamlit · Biopython · Plotly  \n"
        "© 2024"
    )


# ---------------------------------------------------------------------------
# Tab 1 — Sequence Input & Validation
# ---------------------------------------------------------------------------

def _render_input_tab() -> None:
    st.header("📥 Sequence Input & Validation")
    st.markdown(
        "Upload `.fasta` / `.fa` / `.txt` files **or** paste sequences directly. "
        "Accepted types for the selected mode are shown below."
    )

    seq_type: str = st.session_state.seq_type
    char_map = {"DNA": "A C G T N", "RNA": "A C G U N", "Protein": "20 standard amino acids"}
    st.info(f"**{seq_type} mode** — allowed characters: `{char_map[seq_type]}`")

    col1, col2 = st.columns(2, gap="large")

    with col1:
        st.subheader("Sequence 1")
        upload1 = st.file_uploader(
            "Upload FASTA / FA / TXT",
            type=["fasta", "fa", "txt"],
            key=f"upload1_{st.session_state.uploader_key}",  # Dynamic key for resetting
            help="Upload a file containing Sequence 1.",
        )
        text1 = st.text_area(
            "Or paste sequence (FASTA or raw):",
            height=130,
            key="text1",  # Tied to session_state
            placeholder=">seq1\nATGCATGCATGCATGC...",
        )
        
        # Only parse if new input is provided
        if upload1 is not None:
            # FIX: Use .getvalue() instead of .read() to prevent exhausting the byte buffer
            record1 = read_sequence_from_upload(upload1.getvalue(), seq_type)
            st.session_state.record1 = record1
            st.session_state.run_done = False
        elif text1.strip():
            record1 = read_sequence_from_text(text1, seq_type)
            st.session_state.record1 = record1
            st.session_state.run_done = False

        _display_validation(st.session_state.record1, "1")

    with col2:
        st.subheader("Sequence 2")
        upload2 = st.file_uploader(
            "Upload FASTA / FA / TXT",
            type=["fasta", "fa", "txt"],
            key=f"upload2_{st.session_state.uploader_key}",  # Dynamic key for resetting
            help="Upload a file containing Sequence 2.",
        )
        text2 = st.text_area(
            "Or paste sequence (FASTA or raw):",
            height=130,
            key="text2",  # Tied to session_state
            placeholder=">seq2\nATGCATGCATGCATGC...",
        )
        
        # Only parse if new input is provided
        if upload2 is not None:
            # FIX: Use .getvalue() instead of .read()
            record2 = read_sequence_from_upload(upload2.getvalue(), seq_type)
            st.session_state.record2 = record2
            st.session_state.run_done = False
        elif text2.strip():
            record2 = read_sequence_from_text(text2, seq_type)
            st.session_state.record2 = record2
            st.session_state.run_done = False

        _display_validation(st.session_state.record2, "2")

    st.markdown("---")

    # Run alignment button
    rec1: Optional[SequenceRecord] = st.session_state.record1
    rec2: Optional[SequenceRecord] = st.session_state.record2
    both_valid = (
        rec1 is not None and rec1.is_valid
        and rec2 is not None and rec2.is_valid
    )

    run_col, _ = st.columns([1, 3])
    with run_col:
        run_btn = st.button(
            "▶ Run Alignment",
            key="run_btn",
            disabled=not both_valid,
            type="primary",
            use_container_width=True,
        )

    if not both_valid:
        st.caption("⬆ Provide and validate both sequences to enable alignment.")

    if run_btn and both_valid:
        _run_alignment(rec1, rec2)  # type: ignore[arg-type]
        st.success("Alignment complete! Navigate to the other tabs to explore results.")


def _display_validation(record: Optional[SequenceRecord], label: str) -> None:
    """Render a validation status badge for a sequence record."""
    if record is None:
        st.caption("_No sequence loaded yet._")
        return
    if record.is_valid:
        st.success(
            f"✅ **{record.seq_id}** — {record.length:,} residues — valid {record.seq_type}"
        )
        with st.expander("Preview (first 200 chars)"):
            st.code(record.sequence[:200], language=None)
    else:
        st.error(f"❌ Validation failed: {record.error_message}")
        if record.invalid_chars:
            st.caption(f"Invalid characters: {', '.join(record.invalid_chars)}")


# ---------------------------------------------------------------------------
# Alignment runner (with @st.cache_data for performance)
# ---------------------------------------------------------------------------

@st.cache_data(show_spinner=False)
def _compute_alignment(
    seq1: str,
    seq2: str,
    seq1_id: str,
    seq2_id: str,
    mode: Literal["global", "local"],
    seq_type: str,
    match_score: float,
    mismatch_score: float,
    open_gap_score: float,
    extend_gap_score: float,
) -> AlignmentResult:
    """
    Pure cached alignment computation.

    Decorated with ``@st.cache_data`` so the expensive PairwiseAligner call
    is only re-executed when any of the input parameters actually change.
    Switching between tabs or adjusting unrelated UI widgets will NOT
    trigger a re-alignment.

    Parameters mirror :func:`alignment.run_pairwise_alignment`.
    """
    return run_pairwise_alignment(
        seq1=seq1,
        seq2=seq2,
        seq1_id=seq1_id,
        seq2_id=seq2_id,
        mode=mode,
        seq_type=seq_type,
        match_score=match_score,
        mismatch_score=mismatch_score,
        open_gap_score=open_gap_score,
        extend_gap_score=extend_gap_score,
    )


def _run_alignment(rec1: SequenceRecord, rec2: SequenceRecord) -> None:
    """Execute the full analysis pipeline and store results in session state."""
    with st.spinner("Running pairwise alignment…"):
        try:
            result: AlignmentResult = _compute_alignment(
                seq1=rec1.sequence,
                seq2=rec2.sequence,
                seq1_id=rec1.seq_id,
                seq2_id=rec2.seq_id,
                mode=st.session_state.align_mode,
                seq_type=st.session_state.seq_type,
                match_score=st.session_state.match_score,
                mismatch_score=st.session_state.mismatch_score,
                open_gap_score=st.session_state.open_gap_score,
                extend_gap_score=st.session_state.extend_gap_score,
            )
        except ValueError as exc:
            st.error(f"Alignment error: {exc}")
            return

        metrics = calculate_similarity_metrics(result, st.session_state.seq_type)
        mutations_df = detect_mutations(result)
        mut_sum = mutation_summary(mutations_df)

        st.session_state.alignment    = result
        st.session_state.metrics      = metrics
        st.session_state.mutations_df = mutations_df
        st.session_state.mut_summary  = mut_sum
        st.session_state.run_done     = True


# ---------------------------------------------------------------------------
# Tab 2 — Alignment View
# ---------------------------------------------------------------------------

def _render_alignment_tab() -> None:
    st.header("🔗 Pairwise Alignment View")

    if not st.session_state.run_done or st.session_state.alignment is None:
        st.info("Run the alignment from **Tab 1** to see results here.")
        return

    result: AlignmentResult = st.session_state.alignment

    col_a, col_b, col_c = st.columns(3)
    col_a.metric("Alignment Score", f"{result.score:.2f}")
    col_b.metric("Alignment Length", f"{result.alignment_length:,}")
    col_c.metric("Mode", result.mode.capitalize())

    st.markdown("---")
    st.subheader("Alignment Map (colour-coded)")
    _render_colored_alignment(result)

    st.markdown("---")
    st.subheader("Formatted Alignment")
    st.code(result.formatted_alignment, language=None)


def _render_colored_alignment(result: AlignmentResult, chunk: int = 80) -> None:
    """
    Render the gapped alignment as colour-coded HTML spans inside Streamlit.
    Green = match, Red = mismatch, Orange = gap.
    """
    a1 = result.aligned_seq1
    a2 = result.aligned_seq2
    label_len = max(len(result.seq1_id), len(result.seq2_id))

    html_parts: list[str] = [
        "<div style='font-family:monospace;font-size:0.85rem;line-height:1.6;overflow-x:auto;'>"
    ]

    for start in range(0, len(a1), chunk):
        s1_chunk = a1[start : start + chunk]
        s2_chunk = a2[start : start + chunk]

        s1_html = ""
        s2_html = ""
        for ch1, ch2 in zip(s1_chunk, s2_chunk):
            if ch1 == "-" or ch2 == "-":
                s1_html += f"<span style='color:#f39c12'>{ch1}</span>"
                s2_html += f"<span style='color:#f39c12'>{ch2}</span>"
            elif ch1 == ch2:
                s1_html += f"<span style='color:#2ecc71'>{ch1}</span>"
                s2_html += f"<span style='color:#2ecc71'>{ch2}</span>"
            else:
                s1_html += f"<span style='color:#e74c3c'>{ch1}</span>"
                s2_html += f"<span style='color:#e74c3c'>{ch2}</span>"

        lbl1 = result.seq1_id.ljust(label_len)
        lbl2 = result.seq2_id.ljust(label_len)
        pos  = f"{start + 1}".rjust(6)

        html_parts.append(
            f"<div><span style='color:#aaa'>{lbl1} {pos} </span>{s1_html}</div>"
            f"<div><span style='color:#aaa'>{lbl2} {'':>6} </span>{s2_html}</div>"
            "<div style='margin-bottom:8px'></div>"
        )

    html_parts.append("</div>")
    st.markdown("".join(html_parts), unsafe_allow_html=True)

    # Legend
    st.markdown(
        "🟢 **Match** &nbsp;&nbsp; 🔴 **Mismatch** &nbsp;&nbsp; 🟠 **Gap**",
        unsafe_allow_html=True,
    )


# ---------------------------------------------------------------------------
# Tab 3 — Metrics & Mutations
# ---------------------------------------------------------------------------

def _render_metrics_tab() -> None:
    st.header("📊 Similarity Metrics & Mutation Table")

    if not st.session_state.run_done:
        st.info("Run the alignment from **Tab 1** to see results here.")
        return

    metrics: dict = st.session_state.metrics
    mutations_df: pd.DataFrame = st.session_state.mutations_df
    mut_sum: dict  = st.session_state.mut_summary

    # KPI row
    k1, k2, k3, k4 = st.columns(4)
    k1.metric("Identity (%)", f"{metrics['percent_identity']:.2f}%")
    k2.metric("Matches", metrics["match_count"])
    k3.metric("Mismatches", metrics["mismatch_count"])
    k4.metric("Gap Positions", metrics["gap_count"])

    st.markdown("---")
    col_left, col_right = st.columns([1, 1], gap="large")

    with col_left:
        st.subheader("Detailed Metrics")
        rows = metrics_to_display_rows(metrics)
        df_metrics = pd.DataFrame(rows)
        st.dataframe(df_metrics, use_container_width=True, hide_index=True)

    with col_right:
        st.subheader("Mutation Summary")
        summary_rows = [
            {"Mutation Type": k, "Count": v} for k, v in mut_sum.items()
        ]
        st.dataframe(pd.DataFrame(summary_rows), use_container_width=True, hide_index=True)

    st.markdown("---")
    st.subheader(f"Mutation Events Table ({len(mutations_df):,} records)")

    if mutations_df.empty:
        st.success("🎉 No mutations detected — sequences are identical within the aligned region.")
    else:
        # Colour-code the Mutation_Type column
        def _style_mut(val: Any) -> str:
            colors_map = {
                "Substitution": "color:#e74c3c;font-weight:bold",
                "Insertion":    "color:#3498db;font-weight:bold",
                "Deletion":     "color:#f39c12;font-weight:bold",
                "Duplication":  "color:#9b59b6;font-weight:bold",
            }
            return colors_map.get(str(val), "")

        styled = (
            mutations_df.style
            .map(_style_mut, subset="Mutation_Type")  # type: ignore
            .format({"Position": "{:,}"})
        )
        st.dataframe(styled, use_container_width=True, hide_index=True, height=400)
        


# ---------------------------------------------------------------------------
# Tab 4 — Visualizations
# ---------------------------------------------------------------------------

def _render_viz_tab() -> None:
    st.header("📈 Visualizations & Interactive Plots")

    if not st.session_state.run_done or st.session_state.alignment is None:
        st.info("Run the alignment from **Tab 1** to see results here.")
        return

    result: AlignmentResult  = st.session_state.alignment
    mut_sum: dict            = st.session_state.mut_summary
    seq_type: str            = st.session_state.seq_type
    rec1: SequenceRecord     = st.session_state.record1
    rec2: SequenceRecord     = st.session_state.record2

    # Chart 1 — Match / Mismatch Map
    st.subheader("1. Alignment Match / Mismatch Map")
    fig_map = plot_alignment_map(result)
    st.plotly_chart(fig_map, use_container_width=True)

    st.markdown("---")

    # Chart 2 — Dot Plot
    st.subheader("2. Dot Plot")
    window_size = st.slider(
        "K-mer window size (bp / aa)",
        min_value=1, max_value=10, value=1, step=1,
        key="dot_window",
    )
    with st.spinner("Generating dot plot…"):
        fig_dot = plot_dot_plot(
            rec1.sequence,
            rec2.sequence,
            seq1_id=rec1.seq_id,
            seq2_id=rec2.seq_id,
            window=window_size,
        )
    st.plotly_chart(fig_dot, use_container_width=True)

    st.markdown("---")

    # Chart 3 — Base Composition
    st.subheader("3. Base / Residue Composition")
    fig_comp = plot_base_composition(
        rec1.sequence,
        rec2.sequence,
        seq1_id=rec1.seq_id,
        seq2_id=rec2.seq_id,
        seq_type=seq_type,
    )
    st.plotly_chart(fig_comp, use_container_width=True)

    st.markdown("---")

    # Chart 4 — Mutation Distribution Pie
    st.subheader("4. Mutation Type Distribution")
    fig_pie = plot_mutation_distribution(mut_sum)
    st.plotly_chart(fig_pie, use_container_width=True)


# ---------------------------------------------------------------------------
# Tab 5 — Downloads
# ---------------------------------------------------------------------------

def _render_download_tab() -> None:
    st.header("⬇ Download & Report Export")

    if not st.session_state.run_done:
        st.info("Run the alignment from **Tab 1** to generate downloadable reports.")
        return

    result: AlignmentResult  = st.session_state.alignment
    metrics: dict            = st.session_state.metrics
    mutations_df: pd.DataFrame = st.session_state.mutations_df
    mut_sum: dict            = st.session_state.mut_summary
    rec1: SequenceRecord     = st.session_state.record1
    rec2: SequenceRecord     = st.session_state.record2

    col_csv, col_pdf = st.columns(2, gap="large")

    with col_csv:
        st.subheader("📄 Mutation Table (CSV)")
        st.markdown(
            "Downloads the complete mutation events table as a comma-separated file "
            "compatible with Excel, R, or any spreadsheet tool."
        )
        if mutations_df is not None and not mutations_df.empty:
            csv_bytes = export_mutations_csv(mutations_df)
            st.download_button(
                label="⬇ Download mutations.csv",
                data=csv_bytes,
                file_name="mutations.csv",
                mime="text/csv",
                use_container_width=True,
            )
            st.caption(f"{len(mutations_df):,} mutation records · UTF-8 CSV")
        else:
            st.info("No mutation events to export (sequences are identical).")

            # Still allow downloading empty table
            empty_csv = "Position,Seq1_Base,Seq2_Base,Mutation_Type,Description\n".encode()
            st.download_button(
                label="⬇ Download empty mutations.csv",
                data=empty_csv,
                file_name="mutations.csv",
                mime="text/csv",
                use_container_width=True,
            )

    with col_pdf:
        st.subheader("📑 Full Analysis Report (PDF)")
        st.markdown(
            "Generates a multi-page PDF containing metadata, similarity metrics, "
            "alignment preview, embedded alignment-map chart (requires **kaleido**), "
            "and full mutation tables."
        )

        score_params = {
            "match_score":      st.session_state.match_score,
            "mismatch_score":   st.session_state.mismatch_score,
            "open_gap_score":   st.session_state.open_gap_score,
            "extend_gap_score": st.session_state.extend_gap_score,
        }

        with st.spinner("Building PDF…"):
            # Build the alignment map chart for embedding
            alignment_fig = plot_alignment_map(result)
            pdf_bytes = generate_pdf_report(
                seq1_id=rec1.seq_id,
                seq2_id=rec2.seq_id,
                seq_type=st.session_state.seq_type,
                alignment_mode=result.mode,
                metrics=metrics,
                alignment_preview=result.formatted_alignment,
                mutations_df=mutations_df,
                mutation_summary=mut_sum,
                score_params=score_params,
                alignment_fig=alignment_fig,
            )

        st.download_button(
            label="⬇ Download report.pdf",
            data=pdf_bytes,
            file_name="sequence_comparison_report.pdf",
            mime="application/pdf",
            use_container_width=True,
        )
        st.caption("Multi-page A4 PDF · ReportLab")

    st.markdown("---")

    # Inline preview of metrics JSON
    with st.expander("📋 View Raw Metrics (JSON)"):
        import json
        st.code(json.dumps(metrics, indent=2, default=str), language="json")

    with st.expander("📋 View Alignment Metadata"):
        meta = {
            "Sequence 1 ID":        rec1.seq_id,
            "Sequence 1 Length":    rec1.length,
            "Sequence 2 ID":        rec2.seq_id,
            "Sequence 2 Length":    rec2.length,
            "Mode":                 result.mode,
            "Score":                result.score,
            "Alignment Length":     result.alignment_length,
            "Seq1 Span":            f"{result.start_seq1}–{result.end_seq1}",
            "Seq2 Span":            f"{result.start_seq2}–{result.end_seq2}",
        }
        st.json(meta)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    _init_state()
    _render_sidebar()

    # Page header
    st.title("🧬 Genome Sequence Comparison Tool")
    st.markdown(
        "Compare DNA, RNA, or protein sequences with pairwise alignment, "
        "mutation detection, interactive visualizations, and downloadable reports."
    )
    st.markdown("---")

    tab1, tab2, tab3, tab4, tab5 = st.tabs(
        [
            "📥 Sequence Input",
            "🔗 Alignment View",
            "📊 Metrics & Mutations",
            "📈 Visualizations",
            "⬇ Downloads",
        ]
    )

    with tab1:
        _render_input_tab()
    with tab2:
        _render_alignment_tab()
    with tab3:
        _render_metrics_tab()
    with tab4:
        _render_viz_tab()
    with tab5:
        _render_download_tab()


if __name__ == "__main__":
    main()