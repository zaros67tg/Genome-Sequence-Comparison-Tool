"""
visualization.py
----------------
Plotly-based interactive charts for the Genome Sequence Comparison Tool.

Three chart types:
  1. Alignment Match/Mismatch Map — horizontal scatter strip by category.
  2. Dot Plot                     — seq1 vs seq2 position scatter (k-mer window).
  3. Base Composition Bar Chart   — comparative nucleotide / amino acid frequencies.
"""

from __future__ import annotations

from collections import Counter

import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import numpy as np

from alignment import AlignmentResult


# ---------------------------------------------------------------------------
# Color palette
# ---------------------------------------------------------------------------

COLORS = {
    "Match":        "#2ecc71",   # Green
    "Mismatch":     "#e74c3c",   # Red
    "Gap":          "#f39c12",   # Orange
    "Seq1":         "#3498db",   # Blue
    "Seq2":         "#9b59b6",   # Purple
    "background":   "#0e1117",
    "grid":         "#2c2c2c",
}


# ---------------------------------------------------------------------------
# Chart 1 — Alignment Match / Mismatch Map
# ---------------------------------------------------------------------------

def plot_alignment_map(result: AlignmentResult) -> go.Figure:
    """
    Horizontal scatter strip showing the alignment category at each column.

    Positions are colour-coded:
      • Green  (Match)    — identical residues
      • Red    (Mismatch) — substitutions
      • Orange (Gap)      — either sequence has a gap

    Parameters
    ----------
    result : AlignmentResult
        Completed pairwise alignment.

    Returns
    -------
    go.Figure
        Interactive Plotly figure.
    """
    a1 = result.aligned_seq1
    a2 = result.aligned_seq2
    n  = len(a1)

    positions:  list[int] = []
    categories: list[str] = []
    labels:     list[str] = []

    for i, (ch1, ch2) in enumerate(zip(a1, a2), start=1):
        if ch1 == "-" or ch2 == "-":
            cat = "Gap"
        elif ch1 == ch2:
            cat = "Match"
        else:
            cat = "Mismatch"

        positions.append(i)
        categories.append(cat)
        labels.append(f"Pos {i}<br>{result.seq1_id}: {ch1}<br>{result.seq2_id}: {ch2}")

    df = pd.DataFrame({"Position": positions, "Category": categories, "Label": labels})

    color_map = {k: COLORS[k] for k in ["Match", "Mismatch", "Gap"]}

    fig = px.scatter(
        df,
        x="Position",
        y="Category",
        color="Category",
        color_discrete_map=color_map,
        hover_name="Label",
        title="Alignment Match / Mismatch Map",
        labels={"Position": "Alignment Position", "Category": ""},
        height=320,
    )

    fig.update_traces(marker=dict(size=6, opacity=0.85, symbol="square"))
    fig.update_layout(
        plot_bgcolor=COLORS["background"],
        paper_bgcolor=COLORS["background"],
        font_color="#ecf0f1",
        legend_title_text="",
        xaxis=dict(showgrid=True, gridcolor=COLORS["grid"]),
        yaxis=dict(showgrid=False),
        margin=dict(l=10, r=10, t=45, b=30),
    )

    return fig


# ---------------------------------------------------------------------------
# Chart 2 — Dot Plot
# ---------------------------------------------------------------------------

def plot_dot_plot(
    seq1: str,
    seq2: str,
    seq1_id: str = "Seq1",
    seq2_id: str = "Seq2",
    window: int = 1,
) -> go.Figure:
    """
    Generate a dot-plot comparing sequence 1 (x-axis) against sequence 2
    (y-axis).  A dot is placed at (i, j) when the k-mer of length *window*
    starting at position i in seq1 matches the k-mer at position j in seq2.

    Parameters
    ----------
    seq1 : str
        Ungapped first sequence.
    seq2 : str
        Ungapped second sequence.
    seq1_id : str
        X-axis label.
    seq2_id : str
        Y-axis label.
    window : int
        K-mer comparison window (1 = single-base comparison).

    Returns
    -------
    go.Figure
    """
    seq1 = seq1.replace("-", "")
    seq2 = seq2.replace("-", "")

    # Cap to avoid performance issues in the browser
    max_len = 2000
    seq1 = seq1[:max_len]
    seq2 = seq2[:max_len]

    xs: list[int] = []
    ys: list[int] = []

    for i in range(len(seq1) - window + 1):
        kmer1 = seq1[i : i + window]
        for j in range(len(seq2) - window + 1):
            if kmer1 == seq2[j : j + window]:
                xs.append(i + 1)
                ys.append(j + 1)

    fig = go.Figure(
        go.Scattergl(
            x=xs,
            y=ys,
            mode="markers",
            marker=dict(size=2, color=COLORS["Match"], opacity=0.6),
            hovertemplate=(
                f"{seq1_id} pos: %{{x}}<br>{seq2_id} pos: %{{y}}<extra></extra>"
            ),
        )
    )

    fig.update_layout(
        title=f"Dot Plot  (window = {window} bp/aa)",
        xaxis_title=f"{seq1_id} Position",
        yaxis_title=f"{seq2_id} Position",
        plot_bgcolor=COLORS["background"],
        paper_bgcolor=COLORS["background"],
        font_color="#ecf0f1",
        xaxis=dict(showgrid=True, gridcolor=COLORS["grid"]),
        yaxis=dict(showgrid=True, gridcolor=COLORS["grid"]),
        height=500,
        margin=dict(l=10, r=10, t=45, b=30),
    )

    return fig


# ---------------------------------------------------------------------------
# Chart 3 — Base / Residue Composition Bar Chart
# ---------------------------------------------------------------------------

def plot_base_composition(
    seq1: str,
    seq2: str,
    seq1_id: str = "Seq1",
    seq2_id: str = "Seq2",
    seq_type: str = "DNA",
) -> go.Figure:
    """
    Grouped bar chart showing the frequency distribution of nucleotides
    (DNA/RNA) or amino acids (Protein) for both sequences side by side.

    Parameters
    ----------
    seq1 : str
        Ungapped first sequence (uppercase).
    seq2 : str
        Ungapped second sequence (uppercase).
    seq1_id : str
        Legend label for sequence 1.
    seq2_id : str
        Legend label for sequence 2.
    seq_type : str
        'DNA', 'RNA', or 'Protein' — controls which residues to track.

    Returns
    -------
    go.Figure
    """
    seq1 = seq1.replace("-", "").upper()
    seq2 = seq2.replace("-", "").upper()

    if seq_type == "RNA":
        residues = list("ACGU")
    elif seq_type == "Protein":
        residues = list("ACDEFGHIKLMNPQRSTVWY")
    else:  # DNA (default)
        residues = list("ACGT")

    cnt1 = Counter(seq1)
    cnt2 = Counter(seq2)

    len1 = len(seq1) or 1
    len2 = len(seq2) or 1

    freq1 = [cnt1.get(r, 0) / len1 * 100 for r in residues]
    freq2 = [cnt2.get(r, 0) / len2 * 100 for r in residues]

    fig = go.Figure(
        [
            go.Bar(
                name=seq1_id,
                x=residues,
                y=freq1,
                marker_color=COLORS["Seq1"],
                opacity=0.85,
                hovertemplate="%{x}: %{y:.2f}%<extra>" + seq1_id + "</extra>",
            ),
            go.Bar(
                name=seq2_id,
                x=residues,
                y=freq2,
                marker_color=COLORS["Seq2"],
                opacity=0.85,
                hovertemplate="%{x}: %{y:.2f}%<extra>" + seq2_id + "</extra>",
            ),
        ]
    )

    fig.update_layout(
        barmode="group",
        title="Base / Residue Composition",
        xaxis_title="Residue",
        yaxis_title="Frequency (%)",
        plot_bgcolor=COLORS["background"],
        paper_bgcolor=COLORS["background"],
        font_color="#ecf0f1",
        legend_title_text="Sequence",
        xaxis=dict(showgrid=False),
        yaxis=dict(showgrid=True, gridcolor=COLORS["grid"]),
        height=420,
        margin=dict(l=10, r=10, t=45, b=30),
    )

    return fig


# ---------------------------------------------------------------------------
# Chart 4 — Mutation Type Distribution Pie
# ---------------------------------------------------------------------------

def plot_mutation_distribution(summary: dict[str, int]) -> go.Figure:
    """
    Pie chart summarising mutation type counts.

    Parameters
    ----------
    summary : dict[str, int]
        Output of :func:`mutation_detection.mutation_summary`.

    Returns
    -------
    go.Figure
    """
    labels = [k for k, v in summary.items() if v > 0]
    values = [v for v in summary.values() if v > 0]

    if not labels:
        fig = go.Figure()
        fig.update_layout(
            title="Mutation Type Distribution (no mutations detected)",
            paper_bgcolor=COLORS["background"],
            font_color="#ecf0f1",
            height=300,
        )
        return fig

    colors_pie = [
        "#e74c3c",   # Substitution — red
        "#3498db",   # Insertion    — blue
        "#f39c12",   # Deletion     — orange
        "#9b59b6",   # Duplication  — purple
    ]

    fig = go.Figure(
        go.Pie(
            labels=labels,
            values=values,
            hole=0.38,
            marker=dict(colors=colors_pie[: len(labels)], line=dict(color="#ffffff", width=1.5)),
            textinfo="label+percent",
            hovertemplate="%{label}: %{value} events<extra></extra>",
        )
    )

    fig.update_layout(
        title="Mutation Type Distribution",
        paper_bgcolor=COLORS["background"],
        font_color="#ecf0f1",
        height=340,
        margin=dict(l=10, r=10, t=45, b=10),
    )

    return fig
