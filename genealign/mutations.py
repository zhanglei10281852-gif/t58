"""Mutation detection and annotation module."""

from dataclasses import dataclass
from typing import List
from enum import Enum


class MutationType(str, Enum):
    SNP = "SNP"
    INSERTION = "insertion"
    DELETION = "deletion"


@dataclass
class Mutation:
    """Represents a mutation event."""
    position: int
    mutation_type: MutationType
    ref_base: str
    query_base: str
    context: str

    @property
    def label(self) -> str:
        if self.mutation_type == MutationType.SNP:
            return f"{self.ref_base}{self.position + 1}{self.query_base}"
        elif self.mutation_type == MutationType.INSERTION:
            return f"{self.position + 1}ins{self.query_base}"
        else:
            return f"{self.position + 1}del{self.ref_base}"

    def __repr__(self) -> str:
        return f"Mutation({self.label}, type={self.mutation_type.value})"


def detect_mutations(
    ref_aligned: str,
    query_aligned: str,
    ref_start: int,
    context_size: int = 5,
) -> List[Mutation]:
    """Detect mutations from aligned sequences.

    Args:
        ref_aligned: Aligned reference sequence (with gaps).
        query_aligned: Aligned query sequence (with gaps).
        ref_start: Start position in reference (0-indexed).
        context_size: Number of flanking bases to include in context.

    Returns:
        List of Mutation objects sorted by position.
    """
    raw_events: List[dict] = []
    ref_seq_no_gaps = ref_aligned.replace("-", "")

    i = 0
    ref_pos = ref_start

    while i < len(ref_aligned):
        ref_base = ref_aligned[i]
        query_base = query_aligned[i]

        if ref_base == "-":
            inserted_bases = []
            while i < len(ref_aligned) and ref_aligned[i] == "-":
                inserted_bases.append(query_aligned[i])
                i += 1
            insertion_seq = "".join(inserted_bases)
            raw_events.append({
                "type": MutationType.INSERTION,
                "position": ref_pos - 1,
                "ref": "-",
                "alt": insertion_seq,
            })
            continue

        elif query_base == "-":
            deleted_bases = []
            start_ref_pos = ref_pos
            while i < len(ref_aligned) and query_aligned[i] == "-":
                deleted_bases.append(ref_aligned[i])
                ref_pos += 1
                i += 1
            deletion_seq = "".join(deleted_bases)
            raw_events.append({
                "type": MutationType.DELETION,
                "position": start_ref_pos,
                "ref": deletion_seq,
                "alt": "-",
            })
            continue

        elif ref_base != query_base:
            raw_events.append({
                "type": MutationType.SNP,
                "position": ref_pos,
                "ref": ref_base,
                "alt": query_base,
            })
            ref_pos += 1
            i += 1

        else:
            ref_pos += 1
            i += 1

    merged = _merge_adjacent_events(raw_events)

    mutations = []
    for event in merged:
        pos_in_ref = event["position"] - ref_start
        if event["type"] == MutationType.INSERTION:
            context = _get_context(
                ref_seq_no_gaps,
                pos_in_ref,
                pos_in_ref,
                context_size,
                insertion=True,
                inserted=event["alt"],
            )
            mutations.append(Mutation(
                position=event["position"],
                mutation_type=event["type"],
                ref_base=event["ref"],
                query_base=event["alt"],
                context=context,
            ))
        elif event["type"] == MutationType.DELETION:
            end_pos_in_ref = pos_in_ref + len(event["ref"]) - 1
            context = _get_context(
                ref_seq_no_gaps,
                pos_in_ref,
                end_pos_in_ref,
                context_size,
                deletion=True,
                deleted=event["ref"],
            )
            mutations.append(Mutation(
                position=event["position"],
                mutation_type=event["type"],
                ref_base=event["ref"],
                query_base=event["alt"],
                context=context,
            ))
        else:
            context = _get_context(
                ref_seq_no_gaps,
                pos_in_ref,
                pos_in_ref,
                context_size,
            )
            mutations.append(Mutation(
                position=event["position"],
                mutation_type=event["type"],
                ref_base=event["ref"],
                query_base=event["alt"],
                context=context,
            ))

    return mutations


def _merge_adjacent_events(events: List[dict]) -> List[dict]:
    """Merge adjacent insertion or deletion events.

    When the alignment algorithm splits a multi-base indel into multiple
    parts (due to ambiguous alignment), merge them into a single event.

    Args:
        events: List of raw mutation events.

    Returns:
        List of merged mutation events.
    """
    if not events:
        return []

    merged = []
    current = dict(events[0])

    for event in events[1:]:
        if current["type"] == event["type"] and current["type"] in (
            MutationType.INSERTION, MutationType.DELETION
        ):
            distance = event["position"] - current["position"]
            if current["type"] == MutationType.DELETION:
                distance -= len(current["ref"])

            if distance <= 1:
                if current["type"] == MutationType.INSERTION:
                    current["alt"] += event["alt"]
                else:
                    current["ref"] += event["ref"]
                continue

        merged.append(current)
        current = dict(event)

    merged.append(current)
    return merged


def _get_context(
    ref_seq: str,
    start_idx: int,
    end_idx: int,
    context_size: int,
    insertion: bool = False,
    deletion: bool = False,
    inserted: str = "",
    deleted: str = "",
) -> str:
    """Get flanking sequence context around a mutation.

    Args:
        ref_seq: Reference sequence without gaps.
        start_idx: Start index in ref_seq (0-based).
        end_idx: End index in ref_seq (0-based, inclusive).
        context_size: Number of bases on each side.
        insertion: Whether this is an insertion.
        deletion: Whether this is a deletion.
        inserted: Inserted bases (for insertion context).
        deleted: Deleted bases (for deletion display).

    Returns:
        Context string with the mutation in brackets.
    """
    seq_len = len(ref_seq)

    left_start = max(0, start_idx - context_size)
    left = ref_seq[left_start:start_idx] if start_idx > 0 else ""

    right_start = end_idx + 1
    right = ref_seq[right_start:right_start + context_size] if right_start < seq_len else ""

    if insertion:
        middle = f"[{inserted}]"
        return f"{left}{middle}{right}"
    elif deletion:
        middle = f"[{deleted}del]"
        return f"{left}{middle}{right}"
    else:
        middle = ref_seq[start_idx:end_idx + 1]
        return f"{left}[{middle}]{right}"
