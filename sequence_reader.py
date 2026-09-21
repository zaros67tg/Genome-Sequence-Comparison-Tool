"""
sequence_reader.py
------------------
Handles FASTA file parsing and raw text input, validates character composition,
and returns cleaned, uppercase sequences with structured metadata.
"""

from __future__ import annotations

import io
import re
from dataclasses import dataclass, field
from typing import Optional

from Bio import SeqIO


# ---------------------------------------------------------------------------
# Character-set definitions
# ---------------------------------------------------------------------------

DNA_CHARS: frozenset[str] = frozenset("ACGTN")
RNA_CHARS: frozenset[str] = frozenset("ACGUN")
PROTEIN_CHARS: frozenset[str] = frozenset("ACDEFGHIKLMNPQRSTVWY")

SEQ_TYPE_CHARS: dict[str, frozenset[str]] = {
    "DNA": DNA_CHARS,
    "RNA": RNA_CHARS,
    "Protein": PROTEIN_CHARS,
}


# ---------------------------------------------------------------------------
# Data container
# ---------------------------------------------------------------------------

@dataclass
class SequenceRecord:
    """Structured result of a parsed and validated sequence."""

    seq_id: str
    sequence: str
    seq_type: str
    is_valid: bool
    error_message: str = ""
    invalid_chars: list[str] = field(default_factory=list)

    @property
    def length(self) -> int:
        """Return sequence length."""
        return len(self.sequence)


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def _parse_fasta_bytes(raw_bytes: bytes) -> list[tuple[str, str]]:
    """
    Parse raw FASTA bytes using Bio.SeqIO and return a list of
    (record_id, sequence_string) tuples.

    Parameters
    ----------
    raw_bytes : bytes
        Raw file content.

    Returns
    -------
    list[tuple[str, str]]
        Each tuple is (record_id, sequence).
    """
    handle = io.StringIO(raw_bytes.decode("utf-8", errors="replace"))
    records = list(SeqIO.parse(handle, "fasta"))
    return [(rec.id, str(rec.seq)) for rec in records]


def _parse_raw_text(text: str) -> list[tuple[str, str]]:
    """
    Parse a block of raw text that is either:
      - A FASTA-formatted string (starts with '>'), or
      - A bare sequence string (no header).

    Parameters
    ----------
    text : str
        Raw input text.

    Returns
    -------
    list[tuple[str, str]]
        Each tuple is (record_id, sequence).
    """
    text = text.strip()
    if not text:
        return []

    if text.startswith(">"):
        handle = io.StringIO(text)
        records = list(SeqIO.parse(handle, "fasta"))
        return [(rec.id, str(rec.seq)) for rec in records]

    # Bare sequence — strip whitespace and line breaks
    seq = re.sub(r"\s+", "", text)
    return [("manual_input", seq)]


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def _validate_sequence(
    seq_id: str,
    sequence: str,
    seq_type: str,
) -> SequenceRecord:
    """
    Validate that every character in *sequence* belongs to the allowed
    character set for *seq_type*.

    Parameters
    ----------
    seq_id : str
        Identifier for the sequence (used in error messages).
    sequence : str
        Raw sequence string (will be uppercased internally).
    seq_type : str
        One of 'DNA', 'RNA', 'Protein'.

    Returns
    -------
    SequenceRecord
        Populated record with validation status.
    """
    sequence_upper = sequence.upper().strip()

    allowed = SEQ_TYPE_CHARS.get(seq_type)
    if allowed is None:
        return SequenceRecord(
            seq_id=seq_id,
            sequence=sequence_upper,
            seq_type=seq_type,
            is_valid=False,
            error_message=f"Unknown sequence type '{seq_type}'. Choose DNA, RNA, or Protein.",
        )

    if not sequence_upper:
        return SequenceRecord(
            seq_id=seq_id,
            sequence="",
            seq_type=seq_type,
            is_valid=False,
            error_message="Sequence is empty.",
        )

    invalid = sorted({ch for ch in sequence_upper if ch not in allowed})

    if invalid:
        return SequenceRecord(
            seq_id=seq_id,
            sequence=sequence_upper,
            seq_type=seq_type,
            is_valid=False,
            invalid_chars=invalid,
            error_message=(
                f"Invalid characters for {seq_type} sequence: "
                + ", ".join(f"'{c}'" for c in invalid)
            ),
        )

    return SequenceRecord(
        seq_id=seq_id,
        sequence=sequence_upper,
        seq_type=seq_type,
        is_valid=True,
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def read_sequence_from_upload(
    file_bytes: bytes,
    seq_type: str,
    record_index: int = 0,
) -> SequenceRecord:
    """
    Parse a FASTA file uploaded by the user and return the validated record
    at *record_index* (default: first record).

    Parameters
    ----------
    file_bytes : bytes
        Content of the uploaded file.
    seq_type : str
        Expected sequence type ('DNA', 'RNA', 'Protein').
    record_index : int
        Which FASTA record to use if the file contains multiple.

    Returns
    -------
    SequenceRecord
    """
    try:
        records = _parse_fasta_bytes(file_bytes)
    except Exception as exc:
        return SequenceRecord(
            seq_id="unknown",
            sequence="",
            seq_type=seq_type,
            is_valid=False,
            error_message=f"Failed to parse file: {exc}",
        )

    if not records:
        return SequenceRecord(
            seq_id="unknown",
            sequence="",
            seq_type=seq_type,
            is_valid=False,
            error_message="No FASTA records found in the uploaded file.",
        )

    if record_index >= len(records):
        record_index = 0  # Fallback to first record

    seq_id, sequence = records[record_index]
    return _validate_sequence(seq_id, sequence, seq_type)


def read_sequence_from_text(
    text: str,
    seq_type: str,
    record_index: int = 0,
) -> SequenceRecord:
    """
    Parse and validate a sequence from a raw text string pasted by the user.

    Parameters
    ----------
    text : str
        Raw text (FASTA format or bare sequence).
    seq_type : str
        Expected sequence type ('DNA', 'RNA', 'Protein').
    record_index : int
        Which record to use if multiple FASTA records are present.

    Returns
    -------
    SequenceRecord
    """
    try:
        records = _parse_raw_text(text)
    except Exception as exc:
        return SequenceRecord(
            seq_id="unknown",
            sequence="",
            seq_type=seq_type,
            is_valid=False,
            error_message=f"Failed to parse text input: {exc}",
        )

    if not records:
        return SequenceRecord(
            seq_id="unknown",
            sequence="",
            seq_type=seq_type,
            is_valid=False,
            error_message="No sequence content found in the text input.",
        )

    if record_index >= len(records):
        record_index = 0

    seq_id, sequence = records[record_index]
    return _validate_sequence(seq_id, sequence, seq_type)


def get_all_records_from_upload(
    file_bytes: bytes,
    seq_type: str,
) -> list[SequenceRecord]:
    """
    Return validated SequenceRecords for every FASTA entry in an uploaded file.

    Parameters
    ----------
    file_bytes : bytes
        Content of the uploaded file.
    seq_type : str
        Expected sequence type.

    Returns
    -------
    list[SequenceRecord]
    """
    try:
        records = _parse_fasta_bytes(file_bytes)
    except Exception as exc:
        return [
            SequenceRecord(
                seq_id="unknown",
                sequence="",
                seq_type=seq_type,
                is_valid=False,
                error_message=f"Failed to parse file: {exc}",
            )
        ]

    return [_validate_sequence(sid, seq, seq_type) for sid, seq in records]
