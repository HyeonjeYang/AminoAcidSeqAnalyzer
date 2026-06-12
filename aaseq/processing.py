import concurrent.futures
from itertools import islice

import pandas as pd
from tqdm import tqdm
from Bio import SeqIO

from .motifs import getrep
from .properties import calculate_aa_properties


def process_record(record_id, seq, min_len, max_len, min_reps):
    """Worker function for single sequence analysis (for parallel processing)"""
    repeat_counts = getrep(seq, min_len, max_len, min_reps)
    properties = calculate_aa_properties(seq)

    # Merge motif search results and chemical property data
    result_dict = {**repeat_counts, **properties}
    return pd.DataFrame([result_dict], index=[record_id])


def process_fasta_to_matrix(infile, min_len=3, max_len=10, min_reps=2, workers=4, chunk_size=5000):
    all_data = []
    record_iter = SeqIO.parse(infile, "fasta")

    # Implement parallel processing to maximize CPU core utilization using ProcessPoolExecutor
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        chunk_idx = 1
        while True:
            # Load only chunk_size items into memory using islice to prevent overflow of large files
            batch = list(islice(record_iter, chunk_size))
            if not batch:
                break

            futures = []
            for record in batch:
                # Pass record.id and string sequence instead of SeqRecord object for serialization speed
                futures.append(
                    executor.submit(process_record, record.id, str(record.seq), min_len, max_len, min_reps)
                )

            # Apply tqdm to track parallel processing progress by chunk
            for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc=f"Analyzing chunk {chunk_idx}"):
                res = future.result()
                if res is not None:
                    all_data.append(res)
            chunk_idx += 1

    if not all_data:
        print("No sequences processed.")
        return None

    # Do not cast to int as property data contains floats
    matrix = pd.concat(all_data).fillna(0)
    return matrix
