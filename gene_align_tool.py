#!/usr/bin/env python3
"""Gene Alignment Tool - CLI entry point.

A command-line tool for gene sequence alignment and mutation annotation.
Supports Smith-Waterman local alignment, mutation detection, statistics,
and HTML report generation.
"""

import argparse
import sys
import os

from genealign.fasta_reader import read_fasta
from genealign.alignment import SmithWaterman, format_alignment
from genealign.mutations import detect_mutations
from genealign.stats import compute_stats
from genealign.report import generate_html_report


def _load_sequences(ref_path: str, query_path: str):
    """Load reference and query sequences from FASTA files."""
    ref_seqs = read_fasta(ref_path)
    query_seqs = read_fasta(query_path)
    ref_seq = ref_seqs[0]
    query_seq = query_seqs[0]
    return ref_seq, query_seq


def _do_align(ref_seq, query_seq):
    """Run Smith-Waterman alignment and return result."""
    sw = SmithWaterman(match=2, mismatch=-1, gap=-2)
    alignment = sw.align(ref_seq.sequence, query_seq.sequence)
    return alignment


def cmd_align(args):
    """Handle 'align' subcommand."""
    ref_seq, query_seq = _load_sequences(args.reference, args.query)
    alignment = _do_align(ref_seq, query_seq)

    print("=" * 60)
    print("  SMITH-WATERMAN LOCAL ALIGNMENT")
    print("=" * 60)
    print()
    print(f"Reference: {ref_seq.header}")
    print(f"Query:     {query_seq.header}")
    print()
    print(f"Alignment score: {alignment.score}")
    print(f"Alignment length: {alignment.alignment_length} bp")
    print(f"Reference range: {alignment.ref_start + 1} - {alignment.ref_end + 1}")
    print(f"Query range:     {alignment.query_start + 1} - {alignment.query_end + 1}")
    print()
    print("-" * 60)
    print("Alignment:")
    print("-" * 60)
    print(format_alignment(
        alignment.ref_aligned,
        alignment.query_aligned,
        alignment.ref_start,
        alignment.query_start,
        line_width=60,
    ))
    print("=" * 60)


def cmd_mutations(args):
    """Handle 'mutations' subcommand."""
    ref_seq, query_seq = _load_sequences(args.reference, args.query)
    alignment = _do_align(ref_seq, query_seq)
    mutations = detect_mutations(
        alignment.ref_aligned,
        alignment.query_aligned,
        alignment.ref_start,
        context_size=args.context,
    )

    print("=" * 70)
    print("  MUTATION ANALYSIS")
    print("=" * 70)
    print()
    print(f"Reference: {ref_seq.header}")
    print(f"Query:     {query_seq.header}")
    print()
    print(f"Total mutations found: {len(mutations)}")
    print()
    print("-" * 70)
    print(f"{'Pos':>6}  {'Type':<10} {'Ref':<6} {'Query':<8} {'Label':<12} Context")
    print("-" * 70)

    for mut in mutations:
        print(
            f"{mut.position + 1:6d}  "
            f"{mut.mutation_type.value:<10} "
            f"{mut.ref_base:<6} "
            f"{mut.query_base:<8} "
            f"{mut.label:<12} "
            f"{mut.context}"
        )

    print("-" * 70)
    print("=" * 70)


def cmd_stats(args):
    """Handle 'stats' subcommand."""
    ref_seq, query_seq = _load_sequences(args.reference, args.query)
    alignment = _do_align(ref_seq, query_seq)
    mutations = detect_mutations(
        alignment.ref_aligned,
        alignment.query_aligned,
        alignment.ref_start,
        context_size=5,
    )
    stats = compute_stats(ref_seq.sequence, query_seq.sequence, alignment, mutations)

    print("=" * 60)
    print("  ALIGNMENT STATISTICS")
    print("=" * 60)
    print()
    print(f"Reference: {ref_seq.header}")
    print(f"Query:     {query_seq.header}")
    print()

    summary = stats.summary_dict()
    max_key_len = max(len(k) for k in summary.keys())

    for key, value in summary.items():
        print(f"  {key:<{max_key_len}} : {value}")

    print()
    print("=" * 60)


def cmd_report(args):
    """Handle 'report' subcommand."""
    ref_seq, query_seq = _load_sequences(args.reference, args.query)
    alignment = _do_align(ref_seq, query_seq)
    mutations = detect_mutations(
        alignment.ref_aligned,
        alignment.query_aligned,
        alignment.ref_start,
        context_size=5,
    )
    stats = compute_stats(ref_seq.sequence, query_seq.sequence, alignment, mutations)

    output_path = args.output
    generate_html_report(
        ref_header=ref_seq.header,
        query_header=query_seq.header,
        alignment=alignment,
        mutations=mutations,
        stats=stats,
        output_path=output_path,
    )

    print(f"HTML report generated: {os.path.abspath(output_path)}")
    print(f"  - Alignment score: {alignment.score}")
    print(f"  - Mutations found: {len(mutations)}")
    print(f"  - Sequence identity: {stats.identity * 100:.2f}%")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        prog="gene-align",
        description="Gene sequence alignment and mutation annotation tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  gene-align align -r examples/reference.fasta -q examples/query.fasta
  gene-align mutations -r examples/reference.fasta -q examples/query.fasta
  gene-align stats -r examples/reference.fasta -q examples/query.fasta
  gene-align report -r examples/reference.fasta -q examples/query.fasta -o report.html
        """,
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument(
        "-r", "--reference",
        required=True,
        help="Reference sequence FASTA file",
    )
    common_parser.add_argument(
        "-q", "--query",
        required=True,
        help="Query sequence FASTA file",
    )

    align_parser = subparsers.add_parser(
        "align",
        parents=[common_parser],
        help="Smith-Waterman local alignment",
        description="Perform Smith-Waterman local alignment and display results",
    )
    align_parser.set_defaults(func=cmd_align)

    mut_parser = subparsers.add_parser(
        "mutations",
        parents=[common_parser],
        help="Detect mutations from alignment",
        description="Identify SNPs, insertions, and deletions from alignment",
    )
    mut_parser.add_argument(
        "-c", "--context",
        type=int,
        default=5,
        help="Number of flanking bases in context sequence (default: 5)",
    )
    mut_parser.set_defaults(func=cmd_mutations)

    stats_parser = subparsers.add_parser(
        "stats",
        parents=[common_parser],
        help="Alignment statistics",
        description="Compute alignment statistics and metrics",
    )
    stats_parser.set_defaults(func=cmd_stats)

    report_parser = subparsers.add_parser(
        "report",
        parents=[common_parser],
        help="Generate HTML report",
        description="Generate a visual HTML report of alignment results",
    )
    report_parser.add_argument(
        "-o", "--output",
        default="alignment_report.html",
        help="Output HTML file path (default: alignment_report.html)",
    )
    report_parser.set_defaults(func=cmd_report)

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    args.func(args)


if __name__ == "__main__":
    main()
