"""
alignment.py
------------
Pairwise alignment engine wrapping Bio.Align.PairwiseAligner.
Supports Needleman-Wunsch global and Smith-Waterman local modes
with configurable scoring parameters.

For DNA/RNA sequences, flat match/mismatch scores are used.
For protein sequences, the BLOSUM62 substitution matrix is applied
automatically (Bio.Align.substitution_matrices), which correctly penalises
chemically dissimilar residue swaps (e.g. Leu->Pro) more than conservative
substitutions (e.g. Leu->Ile).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

from Bio.Align import PairwiseAligner, substitution_matrices


# ---------------------------------------------------------------------------
# Data container
# ---------------------------------------------------------------------------

@dataclass
class AlignmentResult:
    """Structured result of a pairwise alignment."""

    seq1_id: str
    seq2_id: str
    mode: str
    score: float
    aligned_seq1: str          # Top aligned string with gap characters
    aligned_seq2: str          # Bottom aligned string with gap characters
    alignment_length: int
    start_seq1: int
    end_seq1: int
    start_seq2: int
    end_seq2: int
    formatted_alignment: str   # Human-readable three-line representation
    scoring_method: str = "flat"  # 'flat' for DNA/RNA, 'BLOSUM62' for Protein


# ---------------------------------------------------------------------------
# Alignment core
# ---------------------------------------------------------------------------

def run_pairwise_alignment(
    seq1: str,
    seq2: str,
    seq1_id: str = "Seq1",
    seq2_id: str = "Seq2",
    mode: Literal["global", "local"] = "global",
    seq_type: str = "DNA",
    match_score: float = 2.0,
    mismatch_score: float = -1.0,
    open_gap_score: float = -2.0,
    extend_gap_score: float = -0.5,
) -> AlignmentResult:
    """
    Perform pairwise sequence alignment using Bio.Align.PairwiseAligner.

    Parameters
    ----------
    seq1 : str
        First (query) sequence string, uppercased, no gaps.
    seq2 : str
        Second (target) sequence string, uppercased, no gaps.
    seq1_id : str
        Label for the first sequence.
    seq2_id : str
        Label for the second sequence.
    mode : {'global', 'local'}
        'global' -> Needleman-Wunsch end-to-end alignment.
        'local'  -> Smith-Waterman local sub-sequence alignment.
    seq_type : str
        'DNA', 'RNA', or 'Protein'.  When 'Protein', the BLOSUM62 substitution
        matrix is used and the match_score / mismatch_score parameters are
        ignored.  Gap penalties always apply regardless of seq_type.
    match_score : float
        Score rewarded for identical residues (DNA/RNA only).
    mismatch_score : float
        Penalty applied to substitution positions (DNA/RNA only).
    open_gap_score : float
        Penalty for opening a new gap (all seq types).
    extend_gap_score : float
        Penalty for extending an existing gap by one position (all seq types).

    Returns
    -------
    AlignmentResult
        Full alignment metadata and aligned strings with gap characters.

    Raises
    ------
    ValueError
        If either sequence is empty or the aligner fails to produce alignments.
    """
    if not seq1 or not seq2:
        raise ValueError("Both sequences must be non-empty to perform alignment.")

    aligner = PairwiseAligner()
    aligner.mode = mode

    # --- Scoring scheme -------------------------------------------------
    if seq_type == "Protein":
        # BLOSUM62 is the standard matrix for protein pairwise alignment.
        # It encodes empirical log-odds substitution frequencies, correctly
        # penalising chemically dissimilar swaps (e.g. L->P) far more than
        # conservative ones (e.g. L->I or L->V).
        blosum62 = substitution_matrices.load("BLOSUM62")
        aligner.substitution_matrix = blosum62
        scoring_method = "BLOSUM62"
    else:
        # Flat integer scoring for nucleotide sequences.
        aligner.match_score = match_score
        aligner.mismatch_score = mismatch_score
        scoring_method = "flat"

    # Gap penalties apply to all sequence types.
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score
    # -----------------------------------------------------------------------

    alignments = aligner.align(seq1, seq2)

    try:
        best = next(iter(alignments))
    except StopIteration:
        raise ValueError(
            "No valid alignment found. Check that the sequences are compatible "
            "with the selected alignment mode and scoring scheme."
        )

    # Extract aligned strings
    aligned_seq1, aligned_seq2 = _extract_aligned_strings(best, seq1, seq2)

    # Build human-readable three-line format
    formatted = _format_alignment(aligned_seq1, aligned_seq2, seq1_id, seq2_id)

    # Coordinates
    coords = best.coordinates  # shape (2, n_blocks+1)
    start_seq1 = int(coords[0, 0])
    end_seq1   = int(coords[0, -1])
    start_seq2 = int(coords[1, 0])
    end_seq2   = int(coords[1, -1])

    return AlignmentResult(
        seq1_id=seq1_id,
        seq2_id=seq2_id,
        mode=mode,
        score=float(best.score),
        aligned_seq1=aligned_seq1,
        aligned_seq2=aligned_seq2,
        alignment_length=len(aligned_seq1),
        start_seq1=start_seq1,
        end_seq1=end_seq1,
        start_seq2=start_seq2,
        end_seq2=end_seq2,
        formatted_alignment=formatted,
        scoring_method=scoring_method,
    )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _extract_aligned_strings(alignment, seq1: str, seq2: str) -> tuple[str, str]:
    """
    Convert a Biopython Alignment object to two gapped strings.

    Parameters
    ----------
    alignment
        A Bio.Align.Alignment object.
    seq1 : str
        Original (ungapped) first sequence.
    seq2 : str
        Original (ungapped) second sequence.

    Returns
    -------
    tuple[str, str]
        (gapped_seq1, gapped_seq2)
    """
    # Biopython ≥1.80 exposes __str__ as a multi-line alignment block.
    # We reconstruct the two gapped rows directly from the coordinate blocks.
    coords = alignment.coordinates  # (2, N)

    gapped1: list[str] = []
    gapped2: list[str] = []

    for col in range(coords.shape[1] - 1):
        s1_start, s1_end = int(coords[0, col]), int(coords[0, col + 1])
        s2_start, s2_end = int(coords[1, col]), int(coords[1, col + 1])

        len1 = s1_end - s1_start
        len2 = s2_end - s2_start

        block_len = max(len1, len2)

        if len1 == 0:
            # Gap in seq1
            gapped1.append("-" * len2)
            gapped2.append(seq2[s2_start:s2_end])
        elif len2 == 0:
            # Gap in seq2
            gapped1.append(seq1[s1_start:s1_end])
            gapped2.append("-" * len1)
        else:
            gapped1.append(seq1[s1_start:s1_end])
            gapped2.append(seq2[s2_start:s2_end])

    return "".join(gapped1), "".join(gapped2)


def _format_alignment(
    aligned_seq1: str,
    aligned_seq2: str,
    seq1_id: str,
    seq2_id: str,
    line_width: int = 60,
) -> str:
    """
    Build a classic three-line pairwise alignment string, wrapped at
    *line_width* columns.

    Parameters
    ----------
    aligned_seq1 : str
        Gapped first sequence.
    aligned_seq2 : str
        Gapped second sequence.
    seq1_id : str
        Label for sequence 1.
    seq2_id : str
        Label for sequence 2.
    line_width : int
        Characters per line (excluding the label prefix).

    Returns
    -------
    str
        Multi-line formatted alignment string.
    """
    label_width = max(len(seq1_id), len(seq2_id), 8) + 2
    lines: list[str] = []

    for i in range(0, len(aligned_seq1), line_width):
        chunk1 = aligned_seq1[i : i + line_width]
        chunk2 = aligned_seq2[i : i + line_width]

        # Middle match line
        middle = "".join(
            "|" if a == b and a != "-" else ("." if a != "-" and b != "-" else " ")
            for a, b in zip(chunk1, chunk2)
        )

        lines.append(f"{seq1_id:<{label_width}}{chunk1}")
        lines.append(f"{'':>{label_width}}{middle}")
        lines.append(f"{seq2_id:<{label_width}}{chunk2}")
        lines.append("")  # blank separator

    return "\n".join(lines)
