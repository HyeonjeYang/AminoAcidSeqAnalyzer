import os
from datetime import datetime

from Bio import SeqIO

from .properties import calculate_aa_properties


def print_quick_analysis(sequence):
    """Prints repeat motifs and physicochemical properties for a single sequence."""
    print(f"\nLength: {len(sequence)}")

    properties = calculate_aa_properties(sequence)
    print("\nPhysicochemical properties:")
    for key, value in properties.items():
        print(f"  {key}: {value}")


def find_record(fasta_file, record_id=None):
    """Returns the SeqRecord matching record_id, or the first record if record_id is None."""
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record_id is None or record.id == record_id:
            return record
    return None


def generate_report(output_dir, source, matrix, figures, mut_matrix=None, tables=None):
    """Writes an HTML summary report combining the feature matrix and saved figures."""
    os.makedirs(output_dir, exist_ok=True)
    tables = tables or []

    html_parts = [
        "<html><head><meta charset='utf-8'><title>Amino Acid Sequence Analysis Report</title></head><body>",
        "<h1>Amino Acid Sequence Analysis Report</h1>",
        f"<p>Source: {source}</p>",
        f"<p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>",
        "<h2>Feature Matrix Summary</h2>",
        matrix.head(20).to_html(),
    ]

    for title, img_path in figures:
        html_parts.append(f"<h2>{title}</h2><img src='{os.path.basename(img_path)}' style='max-width:100%;'>")

    for title, table in tables:
        html_parts.append(f"<h2>{title}</h2>")
        html_parts.append(table.to_html(index=False) if hasattr(table, "to_html") else str(table))

    if mut_matrix is not None:
        html_parts.append("<h2>In Silico Mutagenesis Matrix</h2>")
        html_parts.append(mut_matrix.to_html())

    html_parts.append("</body></html>")

    report_path = os.path.join(output_dir, "summary.html")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(html_parts))
    print(f"\nReport generated at '{report_path}'")
