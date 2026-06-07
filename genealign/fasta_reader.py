"""FASTA file reader module."""

from dataclasses import dataclass
from typing import List


@dataclass
class Sequence:
    """Represents a FASTA sequence record."""
    header: str
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)

    def __repr__(self) -> str:
        return f"Sequence(header='{self.header[:30]}...', length={self.length})"


def read_fasta(filepath: str) -> List[Sequence]:
    """Read a FASTA format file and return list of Sequence objects.

    Args:
        filepath: Path to the FASTA file.

    Returns:
        List of Sequence objects.
    """
    sequences: List[Sequence] = []
    current_header = None
    current_seq_parts: List[str] = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    sequences.append(Sequence(
                        header=current_header,
                        sequence="".join(current_seq_parts)
                    ))
                current_header = line[1:].strip()
                current_seq_parts = []
            else:
                current_seq_parts.append(line.upper())

        if current_header is not None:
            sequences.append(Sequence(
                header=current_header,
                sequence="".join(current_seq_parts)
            ))

    if not sequences:
        raise ValueError(f"No sequences found in FASTA file: {filepath}")

    return sequences


def gc_content(sequence: str) -> float:
    """Calculate GC content of a DNA sequence.

    Args:
        sequence: DNA sequence string.

    Returns:
        GC content as a fraction (0.0 to 1.0).
    """
    if not sequence:
        return 0.0
    gc = sum(1 for base in sequence if base in ("G", "C"))
    return gc / len(sequence)
