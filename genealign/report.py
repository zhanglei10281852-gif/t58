"""HTML report generation module."""

from typing import List

from .alignment import AlignmentResult
from .mutations import Mutation, MutationType
from .stats import AlignmentStats


def generate_html_report(
    ref_header: str,
    query_header: str,
    alignment: AlignmentResult,
    mutations: List[Mutation],
    stats: AlignmentStats,
    output_path: str,
) -> None:
    """Generate an HTML visualization report.

    Args:
        ref_header: Reference sequence header.
        query_header: Query sequence header.
        alignment: Alignment result.
        mutations: List of mutations.
        stats: Alignment statistics.
        output_path: Path to save the HTML report.
    """
    sorted_mutations = sorted(mutations, key=lambda m: m.position)

    summary_rows = "".join(
        f"<tr><td>{key}</td><td><strong>{value}</strong></td></tr>"
        for key, value in stats.summary_dict().items()
    )

    mutation_rows = []
    for mut in sorted_mutations:
        if mut.mutation_type == MutationType.SNP:
            bg_class = "snp"
        elif mut.mutation_type == MutationType.INSERTION:
            bg_class = "insertion"
        else:
            bg_class = "deletion"
        mutation_rows.append(
            f"<tr class='{bg_class}'>"
            f"<td>{mut.position + 1}</td>"
            f"<td>{mut.mutation_type.value}</td>"
            f"<td>{mut.ref_base}</td>"
            f"<td>{mut.query_base}</td>"
            f"<td>{mut.label}</td>"
            f"<td class='context'>{mut.context}</td>"
            f"</tr>"
        )
    mutation_rows_html = "".join(mutation_rows)

    similarity_svg = _generate_similarity_svg(alignment, sorted_mutations)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Gene Alignment Report</title>
<style>
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        background-color: #f5f7fa;
        color: #2c3e50;
        line-height: 1.6;
        padding: 30px;
    }}
    .container {{
        max-width: 1100px;
        margin: 0 auto;
        background: white;
        border-radius: 12px;
        box-shadow: 0 4px 20px rgba(0,0,0,0.08);
        padding: 40px;
    }}
    h1 {{
        color: #2c3e50;
        font-size: 28px;
        margin-bottom: 10px;
        border-bottom: 3px solid #3498db;
        padding-bottom: 10px;
    }}
    h2 {{
        color: #2c3e50;
        font-size: 20px;
        margin-top: 30px;
        margin-bottom: 15px;
        padding-left: 12px;
        border-left: 4px solid #3498db;
    }}
    .seq-info {{
        background: #ecf0f1;
        padding: 12px 16px;
        border-radius: 6px;
        margin-bottom: 10px;
        font-size: 13px;
        word-break: break-all;
    }}
    .seq-info strong {{ color: #34495e; }}
    table {{
        width: 100%;
        border-collapse: collapse;
        margin-bottom: 20px;
        font-size: 14px;
    }}
    th, td {{
        padding: 10px 14px;
        text-align: left;
        border-bottom: 1px solid #e0e6ed;
    }}
    th {{
        background-color: #34495e;
        color: white;
        font-weight: 600;
    }}
    tr:nth-child(even) {{ background-color: #f8fafc; }}
    tr.snp {{ background-color: #fff5f5; }}
    tr.snp:hover {{ background-color: #ffe0e0; }}
    tr.insertion {{ background-color: #f0fff4; }}
    tr.insertion:hover {{ background-color: #d4f5dc; }}
    tr.deletion {{ background-color: #fffaf0; }}
    tr.deletion:hover {{ background-color: #fceccf; }}
    .context {{
        font-family: 'Courier New', Courier, monospace;
        font-size: 12px;
        color: #555;
    }}
    .legend {{
        display: flex;
        gap: 20px;
        margin-bottom: 15px;
        flex-wrap: wrap;
    }}
    .legend-item {{
        display: flex;
        align-items: center;
        gap: 8px;
        font-size: 13px;
    }}
    .legend-box {{
        width: 16px;
        height: 16px;
        border-radius: 3px;
    }}
    .snp-color {{ background-color: #e74c3c; }}
    .insertion-color {{ background-color: #27ae60; }}
    .deletion-color {{ background-color: #f39c12; }}
    .svg-container {{
        background: #fafbfc;
        border: 1px solid #e0e6ed;
        border-radius: 6px;
        padding: 16px;
        overflow-x: auto;
    }}
    svg {{ display: block; }}
    .footer {{
        margin-top: 40px;
        padding-top: 20px;
        border-top: 1px solid #e0e6ed;
        color: #95a5a6;
        font-size: 12px;
        text-align: center;
    }}
</style>
</head>
<body>
<div class="container">
    <h1>🧬 Gene Alignment Report</h1>

    <h2>Sequence Information</h2>
    <div class="seq-info"><strong>Reference:</strong> {ref_header}</div>
    <div class="seq-info"><strong>Query:</strong> {query_header}</div>

    <h2>Alignment Summary</h2>
    <table>
        <tr><th>Metric</th><th>Value</th></tr>
        {summary_rows}
    </table>

    <h2>Sequence Similarity Map</h2>
    <div class="legend">
        <div class="legend-item"><div class="legend-box snp-color"></div>SNP</div>
        <div class="legend-item"><div class="legend-box insertion-color"></div>Insertion</div>
        <div class="legend-item"><div class="legend-box deletion-color"></div>Deletion</div>
    </div>
    <div class="svg-container">
        {similarity_svg}
    </div>

    <h2>Mutation Sites ({len(sorted_mutations)})</h2>
    <table>
        <tr>
            <th>Position</th>
            <th>Type</th>
            <th>Ref</th>
            <th>Query</th>
            <th>Label</th>
            <th>Context</th>
        </tr>
        {mutation_rows_html}
    </table>

    <div class="footer">
        Generated by Gene Alignment Tool | Smith-Waterman Local Alignment
    </div>
</div>
</body>
</html>
"""
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)


def _generate_similarity_svg(
    alignment: AlignmentResult,
    mutations: List[Mutation],
    width: int = 900,
    height: int = 80,
) -> str:
    """Generate an inline SVG visualization of alignment similarity.

    Args:
        alignment: Alignment result.
        mutations: List of mutations.
        width: SVG width.
        height: SVG height.

    Returns:
        SVG string.
    """
    align_len = alignment.alignment_length
    if align_len == 0:
        return "<svg></svg>"

    bar_height = 24
    bar_y = (height - bar_height) // 2

    mutation_mark_lines = []
    for mut in mutations:
        rel_pos = (mut.position - alignment.ref_start) / align_len
        x = rel_pos * width

        if mut.mutation_type == MutationType.SNP:
            color = "#e74c3c"
        elif mut.mutation_type == MutationType.INSERTION:
            color = "#27ae60"
        else:
            color = "#f39c12"

        mutation_mark_lines.append(
            f"<line x1='{x:.1f}' y1='{bar_y - 4}' x2='{x:.1f}' y2='{bar_y + bar_height + 4}' "
            f"stroke='{color}' stroke-width='2' />"
        )

    mutation_marks = "".join(mutation_mark_lines)

    svg = f"""<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg">
        <rect x="0" y="{bar_y}" width="{width}" height="{bar_height}" fill="#2ecc71" rx="4" />
        {mutation_marks}
        <text x="0" y="{height - 5}" font-size="11" fill="#7f8c8d">Position {alignment.ref_start + 1}</text>
        <text x="{width}" y="{height - 5}" font-size="11" fill="#7f8c8d" text-anchor="end">Position {alignment.ref_end + 1}</text>
    </svg>"""

    return svg
