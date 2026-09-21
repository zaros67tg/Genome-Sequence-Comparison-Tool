"""
mutation_detection.py
---------------------
Iterates over aligned sequence pairs to locate and categorise coordinate-based
mutation events: substitutions, insertions, deletions, and simple tandem
duplications adjacent to indels.

Output is a structured Pandas DataFrame.
"""

from __future__ import annotations

import pandas as pd

from alignment import AlignmentResult


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MUTATION_TYPES = {
    "SUBSTITUTION": "Substitution",
    "INSERTION":    "Insertion",
    "DELETION":     "Deletion",
    "DUPLICATION":  "Duplication",
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def detect_mutations(result: AlignmentResult) -> pd.DataFrame:
    """
    Walk through two aligned strings (with gaps) and classify every
    non-match position as a mutation event.

    Duplication detection: when an insertion or deletion run is preceded
    or followed by an identical sequence block (simple tandem repeat signal),
    the event is re-labelled as a Duplication.

    Parameters
    ----------
    result : AlignmentResult
        Output of :func:`alignment.run_pairwise_alignment`.  The
        ``aligned_seq1`` and ``aligned_seq2`` attributes must already
        contain the gap characters.

    Returns
    -------
    pd.DataFrame
        Columns:
            - Position       : int   — 1-based alignment column index
            - Seq1_Base      : str   — residue in sequence 1 ('-' = gap)
            - Seq2_Base      : str   — residue in sequence 2 ('-' = gap)
            - Mutation_Type  : str   — Substitution / Insertion / Deletion /
                                       Duplication
            - Description    : str   — Human-readable description
    """
    a1 = result.aligned_seq1
    a2 = result.aligned_seq2

    if len(a1) != len(a2):
        raise ValueError(
            "Aligned sequences must have equal length.  "
            f"Got {len(a1)} vs {len(a2)}."
        )

    rows: list[dict[str, object]] = []

    for pos, (ch1, ch2) in enumerate(zip(a1, a2), start=1):
        if ch1 == ch2:
            continue  # Match — not a mutation

        mut_type, description = _classify_event(ch1, ch2, pos)
        rows.append(
            {
                "Position":      pos,
                "Seq1_Base":     ch1,
                "Seq2_Base":     ch2,
                "Mutation_Type": mut_type,
                "Description":   description,
            }
        )

    df = pd.DataFrame(
        rows,
        columns=["Position", "Seq1_Base", "Seq2_Base", "Mutation_Type", "Description"],
    )

    if df.empty:
        return df

    # Post-process: flag duplications by scanning indel runs
    df = _flag_duplications(df, a1, a2)

    return df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _classify_event(ch1: str, ch2: str, pos: int) -> tuple[str, str]:
    """
    Classify a single column mismatch.

    Parameters
    ----------
    ch1 : str   — Character in aligned sequence 1.
    ch2 : str   — Character in aligned sequence 2.
    pos : int   — 1-based column index.

    Returns
    -------
    tuple[str, str]
        (mutation_type_label, description_string)
    """
    if ch1 == "-" and ch2 != "-":
        return (
            MUTATION_TYPES["INSERTION"],
            f"Position {pos}: Insertion of '{ch2}' in Seq2 (gap in Seq1).",
        )
    if ch1 != "-" and ch2 == "-":
        return (
            MUTATION_TYPES["DELETION"],
            f"Position {pos}: Deletion of '{ch1}' from Seq1 (gap in Seq2).",
        )
    # Both non-gap → substitution
    return (
        MUTATION_TYPES["SUBSTITUTION"],
        f"Position {pos}: Substitution {ch1}->{ch2}.",
    )


def _flag_duplications(
    df: pd.DataFrame,
    a1: str,
    a2: str,
    context_window: int = 6,
    min_repeat_len: int = 2,
) -> pd.DataFrame:
    """
    Re-label Insertion / Deletion rows as Duplication where the inserted
    or deleted base is part of a tandem repeat motif flanking the indel.

    Heuristic
    ---------
    A tandem duplication is inferred when the inserted/deleted residue
    appears as a consecutive run of **at least** *min_repeat_len* identical
    characters in the *context_window* bases immediately before **or**
    immediately after the indel position (in the non-gapped strand).

    This avoids false-positive labelling from coincidental single-base
    repeats (e.g. poly-A background at an unrelated indel).  The threshold
    of 2 is a conservative lower bound; real tandem duplications in clinical
    genomics typically span ≥ 2 bp.

    Parameters
    ----------
    df : pd.DataFrame
        Raw mutation table from :func:`detect_mutations`.
    a1 : str
        Full gapped sequence 1.
    a2 : str
        Full gapped sequence 2.
    context_window : int
        Number of flanking positions to inspect for the repeat motif.
    min_repeat_len : int
        Minimum run of identical consecutive residues required to call
        a duplication (default 2).  Must be >= 2 to avoid false positives.

    Returns
    -------
    pd.DataFrame
        Updated DataFrame with Duplication labels where applicable.
    """
    if min_repeat_len < 2:
        min_repeat_len = 2  # Safety floor — 1 would flag every repeat base

    df = df.copy()

    for idx, row in df.iterrows():
        mut = row["Mutation_Type"]
        if mut not in (MUTATION_TYPES["INSERTION"], MUTATION_TYPES["DELETION"]):
            continue

        pos0 = int(row["Position"]) - 1  # Convert to 0-based alignment index

        if mut == MUTATION_TYPES["INSERTION"]:
            inserted_base = row["Seq2_Base"]
            # Inspect the non-gapped seq2 context flanking this column
            ctx_before = a2[max(0, pos0 - context_window) : pos0].replace("-", "")
            ctx_after  = a2[pos0 + 1 : pos0 + 1 + context_window].replace("-", "")
        else:  # DELETION
            inserted_base = row["Seq1_Base"]
            ctx_before = a1[max(0, pos0 - context_window) : pos0].replace("-", "")
            ctx_after  = a1[pos0 + 1 : pos0 + 1 + context_window].replace("-", "")

        # Count the trailing run of the same base immediately before the gap
        run_before = _trailing_run(ctx_before, inserted_base)
        # Count the leading run of the same base immediately after the gap
        run_after  = _leading_run(ctx_after, inserted_base)

        if run_before >= min_repeat_len or run_after >= min_repeat_len:
            df.at[idx, "Mutation_Type"] = MUTATION_TYPES["DUPLICATION"]
            old_desc = str(row["Description"])
            df.at[idx, "Description"] = (
                old_desc.replace(
                    "Insertion" if mut == MUTATION_TYPES["INSERTION"] else "Deletion",
                    "Duplication",
                )
                + f" (tandem repeat: run of {max(run_before, run_after)} '{inserted_base}')"
            )

    return df


def _trailing_run(text: str, base: str) -> int:
    """Count how many times *base* appears consecutively at the END of *text*."""
    count = 0
    for ch in reversed(text):
        if ch == base:
            count += 1
        else:
            break
    return count


def _leading_run(text: str, base: str) -> int:
    """Count how many times *base* appears consecutively at the START of *text*."""
    count = 0
    for ch in text:
        if ch == base:
            count += 1
        else:
            break
    return count


def mutation_summary(df: pd.DataFrame) -> dict[str, int]:
    """
    Aggregate mutation counts by type.

    Parameters
    ----------
    df : pd.DataFrame
        Output of :func:`detect_mutations`.

    Returns
    -------
    dict[str, int]
        Keys are mutation type labels; values are counts.
    """
    if df.empty:
        return {v: 0 for v in MUTATION_TYPES.values()}
    counts = df["Mutation_Type"].value_counts().to_dict()
    return {v: counts.get(v, 0) for v in MUTATION_TYPES.values()}
