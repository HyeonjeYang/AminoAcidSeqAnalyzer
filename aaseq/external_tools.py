import os
import shutil
import subprocess


class ExternalToolError(RuntimeError):
    pass


def require_executable(executable):
    path = shutil.which(executable)
    if path is None:
        raise ExternalToolError(
            f"Required executable '{executable}' was not found on PATH. "
            "Install the tool separately and rerun this command."
        )
    return path


def run_mmseqs_easy_cluster(
    fasta_file,
    output_prefix,
    tmp_dir,
    min_seq_id=None,
    coverage=None,
    sensitivity=None,
    mmseqs_bin="mmseqs",
    extra_args=None,
):
    """Runs MMseqs2 easy-cluster through an external executable.

    This function does not bundle or reimplement MMseqs2. It only calls a
    user-installed binary and returns the expected output paths.
    """
    mmseqs_path = require_executable(mmseqs_bin)
    os.makedirs(tmp_dir, exist_ok=True)

    command = [mmseqs_path, "easy-cluster", fasta_file, output_prefix, tmp_dir]
    if min_seq_id is not None:
        command.extend(["--min-seq-id", str(min_seq_id)])
    if coverage is not None:
        command.extend(["-c", str(coverage)])
    if sensitivity is not None:
        command.extend(["-s", str(sensitivity)])
    if extra_args:
        command.extend(extra_args)

    subprocess.run(command, check=True)
    return {
        "cluster_tsv": f"{output_prefix}_cluster.tsv",
        "representatives_fasta": f"{output_prefix}_rep_seq.fasta",
        "all_sequences_fasta": f"{output_prefix}_all_seqs.fasta",
    }


def run_foldseek_easy_search(
    query_structures,
    target_database,
    output_tsv,
    tmp_dir,
    foldseek_bin="foldseek",
    extra_args=None,
):
    """Runs Foldseek easy-search through an external executable.

    Foldseek is GPL-licensed upstream, so this project intentionally does not
    vendor or copy its source code. Users who install Foldseek separately can
    call it from this wrapper.
    """
    foldseek_path = require_executable(foldseek_bin)
    os.makedirs(tmp_dir, exist_ok=True)

    command = [foldseek_path, "easy-search", query_structures, target_database, output_tsv, tmp_dir]
    if extra_args:
        command.extend(extra_args)
    subprocess.run(command, check=True)
    return output_tsv
