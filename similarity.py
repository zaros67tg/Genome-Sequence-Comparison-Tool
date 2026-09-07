"""
similarity.py
-------------
Computes comprehensive pairwise similarity metrics from aligned sequence strings.
Returns a structured dictionary suitable for display and report generation.
"""

from __future__ import annotations

from typing import Optional

from alignment import AlignmentResult


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def calculate_similarity_metrics(
    result: AlignmentResult,
    seq_type: str = "DNA",
) -> dict[str, object]:
    """
    Derive all similarity and composition metrics from an AlignmentResult.

    Parameters
    ----------
    result : AlignmentResult
        Output of :func:`alignment.run_pairwise_alignment`.
    seq_type : str
        Sequence type ('DNA', 'RNA', or 'Protein').  GC content is only
        meaningful for DNA / RNA sequences.

    Returns
    -------
    dict[str, object]
        Keys (all human-readable strings):
            - ``seq1_length``          : int
            - ``seq2_length``          : int
            - ``alignment_length``     : int
            - ``alignment_score``      : float
            - ``match_count``          : int
            - ``mismatch_count``       : int
            - ``gap_count``            : int
            - ``gap_open_count``       : int
            - ``percent_identity``     : float (0–100)
            - ``gc_content_seq1``      : float | None
            - ``gc_content_seq2``      : float | None
    """
    a1 = result.aligned_seq1
    a2 = result.aligned_seq2

    # Raw lengths (without gaps)
    seq1_length = len(a1.replace("-", ""))
    seq2_length = len(a2.replace("-", ""))
    alignment_length = result.alignment_length

    # Column-wise statistics
    match_count = 0
    mismatch_count = 0
    gap_count = 0
    gap_open_count = 0
    in_gap = False

    for ch1, ch2 in zip(a1, a2):
        if ch1 == "-" or ch2 == "-":
            gap_count += 1
            if not in_gap:
                gap_open_count += 1
                in_gap = True
        else:
            in_gap = False
            if ch1 == ch2:
                match_count += 1
            else:
                mismatch_count += 1

    # Percentage identity (BLAST-style: matches over alignment length)
    percent_identity = (match_count / alignment_length * 100) if alignment_length > 0 else 0.0

    # GC content (DNA / RNA only)
    gc1: Optional[float] = None
    gc2: Optional[float] = None
    if seq_type in ("DNA", "RNA"):
        gc1 = _gc_content(a1.replace("-", ""))
        gc2 = _gc_content(a2.replace("-", ""))

    return {
        "seq1_length": seq1_length,
        "seq2_length": seq2_length,
        "alignment_length": alignment_length,
        "alignment_score": result.score,
        "match_count": match_count,
        "mismatch_count": mismatch_count,
        "gap_count": gap_count,
        "gap_open_count": gap_open_count,
        "percent_identity": round(percent_identity, 4),
        "gc_content_seq1": round(gc1, 4) if gc1 is not None else None,
        "gc_content_seq2": round(gc2, 4) if gc2 is not None else None,
    }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _gc_content(sequence: str) -> float:
    """
    Calculate GC content (%) for a nucleotide sequence.

    Parameters
    ----------
    sequence : str
        Ungapped nucleotide sequence (uppercase, no gaps).

    Returns
    -------
    float
        Percentage of G + C bases (0–100).  Returns 0 for empty input.
    """
    if not sequence:
        return 0.0
    gc = sum(1 for ch in sequence if ch in "GC")
    return (gc / len(sequence)) * 100


def metrics_to_display_rows(metrics: dict[str, object]) -> list[dict[str, str]]:
    """
    Convert the metrics dictionary to a list of display-friendly rows
    for rendering in a Streamlit table or DataFrame.

    Parameters
    ----------
    metrics : dict[str, object]
        Output of :func:`calculate_similarity_metrics`.

    Returns
    -------
    list[dict[str, str]]
        Each element has keys 'Metric' and 'Value'.
    """
    labels = {
        "seq1_length":       "Sequence 1 Length (bp / aa)",
        "seq2_length":       "Sequence 2 Length (bp / aa)",
        "alignment_length":  "Alignment Length",
        "alignment_score":   "Alignment Score",
        "match_count":       "Matches",
        "mismatch_count":    "Mismatches",
        "gap_count":         "Total Gap Positions",
        "gap_open_count":    "Gap Opens (Indel Events)",
        "percent_identity":  "Percent Identity (%)",
        "gc_content_seq1":   "GC Content — Seq 1 (%)",
        "gc_content_seq2":   "GC Content — Seq 2 (%)",
    }

    rows = []
    for key, label in labels.items():
        value = metrics.get(key)
        if value is None:
            continue  # Skip GC for proteins
        if isinstance(value, float):
            display_value = f"{value:.4f}"
        else:
            display_value = str(value)
        rows.append({"Metric": label, "Value": display_value})

    return rows
