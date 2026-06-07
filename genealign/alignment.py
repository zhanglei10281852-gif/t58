"""Smith-Waterman local alignment algorithm implementation."""

from dataclasses import dataclass
from typing import Tuple, List


@dataclass
class AlignmentResult:
    """Result of a Smith-Waterman local alignment."""
    score: int
    ref_aligned: str
    query_aligned: str
    ref_start: int
    ref_end: int
    query_start: int
    query_end: int

    @property
    def alignment_length(self) -> int:
        return len(self.ref_aligned)

    @property
    def matches(self) -> int:
        return sum(1 for a, b in zip(self.ref_aligned, self.query_aligned)
                   if a == b and a != "-" and b != "-")

    @property
    def mismatches(self) -> int:
        return sum(1 for a, b in zip(self.ref_aligned, self.query_aligned)
                   if a != b and a != "-" and b != "-")

    @property
    def insertions(self) -> int:
        return sum(1 for a, b in zip(self.ref_aligned, self.query_aligned)
                   if a == "-" and b != "-")

    @property
    def deletions(self) -> int:
        return sum(1 for a, b in zip(self.ref_aligned, self.query_aligned)
                   if a != "-" and b == "-")


class SmithWaterman:
    """Smith-Waterman local alignment algorithm.

    Scoring:
        Match: +2
        Mismatch: -1
        Gap penalty: -2 (linear)
    """

    def __init__(self, match: int = 2, mismatch: int = -1, gap: int = -2):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap

    def align(self, ref_seq: str, query_seq: str) -> AlignmentResult:
        """Perform local alignment between two sequences.

        Args:
            ref_seq: Reference sequence.
            query_seq: Query sequence.

        Returns:
            AlignmentResult with the optimal local alignment.
        """
        m = len(ref_seq)
        n = len(query_seq)

        score_matrix = [[0] * (n + 1) for _ in range(m + 1)]
        traceback = [[None] * (n + 1) for _ in range(m + 1)]

        max_score = 0
        max_i, max_j = 0, 0

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if ref_seq[i - 1] == query_seq[j - 1]:
                    diag_score = score_matrix[i - 1][j - 1] + self.match
                else:
                    diag_score = score_matrix[i - 1][j - 1] + self.mismatch

                up_score = score_matrix[i - 1][j] + self.gap
                left_score = score_matrix[i][j - 1] + self.gap

                best = max(0, diag_score, up_score, left_score)
                score_matrix[i][j] = best

                if best == 0:
                    traceback[i][j] = None
                elif best == diag_score:
                    traceback[i][j] = "diag"
                elif best == up_score:
                    traceback[i][j] = "up"
                else:
                    traceback[i][j] = "left"

                if best > max_score:
                    max_score = best
                    max_i, max_j = i, j

        ref_aligned: List[str] = []
        query_aligned: List[str] = []
        i, j = max_i, max_j

        while traceback[i][j] is not None:
            if traceback[i][j] == "diag":
                ref_aligned.append(ref_seq[i - 1])
                query_aligned.append(query_seq[j - 1])
                i -= 1
                j -= 1
            elif traceback[i][j] == "up":
                ref_aligned.append(ref_seq[i - 1])
                query_aligned.append("-")
                i -= 1
            else:
                ref_aligned.append("-")
                query_aligned.append(query_seq[j - 1])
                j -= 1

        ref_aligned_str = "".join(reversed(ref_aligned))
        query_aligned_str = "".join(reversed(query_aligned))

        return AlignmentResult(
            score=max_score,
            ref_aligned=ref_aligned_str,
            query_aligned=query_aligned_str,
            ref_start=i,
            ref_end=max_i - 1,
            query_start=j,
            query_end=max_j - 1,
        )


def format_alignment(
    ref_aligned: str,
    query_aligned: str,
    ref_start: int,
    query_start: int,
    line_width: int = 60,
) -> str:
    """Format alignment output in BLAST-like style.

    Args:
        ref_aligned: Aligned reference sequence (with gaps).
        query_aligned: Aligned query sequence (with gaps).
        ref_start: Start position in reference (0-indexed).
        query_start: Start position in query (0-indexed).
        line_width: Number of characters per line.

    Returns:
        Formatted alignment string.
    """
    lines = []
    total_len = len(ref_aligned)

    ref_pos = ref_start
    query_pos = query_start

    for start in range(0, total_len, line_width):
        end = min(start + line_width, total_len)
        ref_segment = ref_aligned[start:end]
        query_segment = query_aligned[start:end]

        match_line = ""
        for a, b in zip(ref_segment, query_segment):
            if a == b and a != "-":
                match_line += "|"
            elif a == "-" or b == "-":
                match_line += " "
            else:
                match_line += "x"

        ref_bases = sum(1 for c in ref_segment if c != "-")
        query_bases = sum(1 for c in query_segment if c != "-")

        ref_end_pos = ref_pos + ref_bases - 1
        query_end_pos = query_pos + query_bases - 1

        lines.append(f"Ref   {ref_pos + 1:4d}  {ref_segment}  {ref_end_pos + 1:4d}")
        lines.append(f"          {match_line}")
        lines.append(f"Query {query_pos + 1:4d}  {query_segment}  {query_end_pos + 1:4d}")
        lines.append("")

        ref_pos += ref_bases
        query_pos += query_bases

    return "\n".join(lines)
