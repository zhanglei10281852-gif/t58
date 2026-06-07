"""Alignment statistics module."""

from dataclasses import dataclass
from typing import List

from .alignment import AlignmentResult
from .mutations import Mutation, MutationType
from .fasta_reader import gc_content


@dataclass
class AlignmentStats:
    """Statistics for an alignment."""
    ref_length: int
    query_length: int
    ref_gc: float
    query_gc: float
    alignment_length: int
    alignment_score: int
    coverage: float
    identity: float
    num_mutations: int
    num_snps: int
    num_insertions: int
    num_deletions: int
    mutation_density: float

    def summary_dict(self) -> dict:
        return {
            "Reference length (bp)": self.ref_length,
            "Query length (bp)": self.query_length,
            "Reference GC content": f"{self.ref_gc * 100:.2f}%",
            "Query GC content": f"{self.query_gc * 100:.2f}%",
            "Alignment length (bp)": self.alignment_length,
            "Alignment score": self.alignment_score,
            "Coverage": f"{self.coverage * 100:.2f}%",
            "Sequence identity": f"{self.identity * 100:.2f}%",
            "Total mutations": self.num_mutations,
            "SNPs": self.num_snps,
            "Insertions": self.num_insertions,
            "Deletions": self.num_deletions,
            "Mutation density (per kb)": f"{self.mutation_density:.2f}",
        }


def compute_stats(
    ref_seq: str,
    query_seq: str,
    alignment: AlignmentResult,
    mutations: List[Mutation],
) -> AlignmentStats:
    """Compute alignment statistics.

    Args:
        ref_seq: Full reference sequence.
        query_seq: Full query sequence.
        alignment: Alignment result.
        mutations: List of detected mutations.

    Returns:
        AlignmentStats object.
    """
    ref_len = len(ref_seq)
    query_len = len(query_seq)
    ref_gc = gc_content(ref_seq)
    query_gc = gc_content(query_seq)

    align_len = alignment.alignment_length
    matches = alignment.matches
    score = alignment.score

    ref_covered_bases = alignment.ref_end - alignment.ref_start + 1
    coverage = ref_covered_bases / ref_len if ref_len > 0 else 0.0
    identity = matches / align_len if align_len > 0 else 0.0

    num_snps = sum(1 for m in mutations if m.mutation_type == MutationType.SNP)
    num_insertions = sum(1 for m in mutations if m.mutation_type == MutationType.INSERTION)
    num_deletions = sum(1 for m in mutations if m.mutation_type == MutationType.DELETION)
    num_mutations = len(mutations)

    mutation_density = (num_mutations / align_len * 1000) if align_len > 0 else 0.0

    return AlignmentStats(
        ref_length=ref_len,
        query_length=query_len,
        ref_gc=ref_gc,
        query_gc=query_gc,
        alignment_length=align_len,
        alignment_score=score,
        coverage=coverage,
        identity=identity,
        num_mutations=num_mutations,
        num_snps=num_snps,
        num_insertions=num_insertions,
        num_deletions=num_deletions,
        mutation_density=mutation_density,
    )
