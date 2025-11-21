#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import gzip
import importlib.util
from Bio.Seq import Seq
import math
import re

def print_logo():
    logo = r"""
🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬
██╗   ██╗██╗██████╗   ██████╗         ██████╗ ███████╗████████╗██████║ ██║  ██║
██║   ██║██║██╔══██╗ ██╔═══██╗        ██╔═══╝ ██╔══██╗╚══██╔══╝██╔═══  ██║  ██║
██║   ██║██║███████║ ██║   ██║🧬🦠🧬 ██║     ███████║   ██║   ██║     ███████║   
╚██╗ ██╔╝██║██║  ██║ ██║   ██║        ██║     ██║  ██║   ██║   ██║     ██║  ██║ 
 ╚████╔╝ ██║██║  ██║ ████████║        ██████  ██║  ██║   ██║   ██████║ ██║  ██║
  ╚═══╝  ╚═╝╚═    ═╝ ╚═══════╝       ╚══════╝╚══╝ ╚══╝   ╚═╝  ╚══════╝ ══╝  ══╝
🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬 V.1.1 by Morán et al. 2025 🦠🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬
"""
    print(logo)
    print(" 🧬 Detection and identification of DNA viruses from ONT, RCA-R-ONT, and Illumina data 🧬\n")

def check_dependencies():
    print("🔍 Checking dependencies...")
    required_python_modules = ["Bio"]
    required_commands = ["fastp", "fastplong", "kraken2", "flye", "python", "blastn", "awk", "spades.py"]

    for module in required_python_modules:
        if importlib.util.find_spec(module) is None:
            print(f"❌ Missing module: {module}. Install it with pip.")
            sys.exit(1)

    for cmd in required_commands:
        if not shutil.which(cmd):
            print(f"❌ Missing executable in PATH: {cmd}")
            sys.exit(1)

    global extract_kraken_script_path
    extract_script = Path("extract_kraken_reads.py")
    if not extract_script.exists():
        extract_script = Path(__file__).parent / "extract_kraken_reads.py"
        if not extract_script.exists():
            print("❌ Missing extract_kraken_reads.py")
            sys.exit(1)

    extract_kraken_script_path = str(extract_script.resolve())
    print("✅ Dependencies verified.\n")

def run_command(command, description, suppress_output=False, check_existing_output=None):
    print(f"🔧 Running: {description}")
    if check_existing_output and Path(check_existing_output).exists() and Path(check_existing_output).stat().st_size > 0:
        print(f"⏭️ {description}: already exists, skipping.")
        return
    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        if not suppress_output:
            print(result.stdout)
            print(result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"❌ Error while running: {e.cmd}\n{e.stderr}")
        sys.exit(1)

########################################
# Kraken2 runner
########################################

def run_kraken2(input_fastq, output_report, output_output, db_path, threads):
    global memory_mapping

    if output_report.exists() and output_report.stat().st_size > 0 and \
       output_output.exists() and output_output.stat().st_size > 0:
        print(f"⏭️ Kraken2 already done: {input_fastq.name}")
        return

    if not db_path.exists():
        print(f"❌ Kraken DB path does not exist: {db_path}")
        sys.exit(1)

    if not any(db_path.glob("*.k2d")):
        print(f"❌ Kraken DB incomplete: no *.k2d files found in {db_path}")
        sys.exit(1)

    cmd = [
        "kraken2",
        "--db", str(db_path),
        "--report", str(output_report),
        "--output", str(output_output),
        "--threads", str(threads),
        str(input_fastq)
    ]

    if memory_mapping:
        cmd.insert(1, "--memory-mapping")

    run_command(cmd, f"Kraken2 on {input_fastq.name}")

########################################
# Extract Kraken reads
########################################

def extract_kraken_reads(input_output, input_fastq, output_fastq):
    cmd = [
        "python", extract_kraken_script_path,
        "-k", str(input_output),
        "-s", str(input_fastq),
        "-o", str(output_fastq),
        "-t", "2", "2157", "2759",
        "--exclude",
        "--fastq-output"
    ]
    run_command(cmd, "Extracting non-host non-bacterial reads", suppress_output=True)

########################################
# FASTQ opening util
########################################

def open_fastq(filepath):
    return gzip.open(filepath, "rt") if str(filepath).endswith(".gz") else open(filepath, "r")

########################################
# Low complexity filter
########################################

def is_low_complexity(seq, entropy_threshold=1.8):
    freqs = [seq.count(n)/len(seq) for n in "ATCG"]
    entropy = -sum(f * math.log2(f) for f in freqs if f > 0)
    return entropy < entropy_threshold

def filter_low_complexity_reads(input_fastq, output_fastq):
    total = kept = discarded = 0
    with open_fastq(input_fastq) as infile, open(output_fastq, "w") as outfile:
        for record in SeqIO.parse(infile, "fastq"):
            total += 1
            if not is_low_complexity(str(record.seq)):
                SeqIO.write(record, outfile, "fastq")
                kept += 1
            else:
                discarded += 1
    print(f"🧼 Complexity filter: {kept} kept, {discarded} removed, total {total}")

########################################
# CLAMP QUALITIES (new)
########################################

def clamp_qualities(input_fastq, output_fastq, max_phred=40):
    """
    SPAdes only accepts qualities ≤ Q40 (ASCII 73).
    This truncates everything above Q40.
    """
    with open(input_fastq) as fin, open(output_fastq, "w") as fout:
        lines = fin.read().splitlines()
        for i in range(0, len(lines), 4):
            fout.write(lines[i] + "\n")
            fout.write(lines[i+1] + "\n")
            fout.write(lines[i+2] + "\n")
            qual = lines[i+3]
            clamped = ''.join(chr(min(ord(q), 33 + max_phred)) for q in qual)
            fout.write(clamped + "\n")
########################################
# Fragmentation for RCA-R
########################################

def fragment_by_motif(input_fastq, output_fastq, motif="TAATATTA"):
    print(f"✂️  Fragmenting by motif '{motif}' and RC...")

    motif = motif.upper()
    motif_rc = str(Seq(motif).reverse_complement())
    motif_len = len(motif)

    total_reads = 0
    reads_with_motif = 0
    total_fragments = []

    with open_fastq(input_fastq) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            total_reads += 1
            seq = str(record.seq).upper()

            positions = []
            for i in range(len(seq)):
                if seq.startswith(motif, i) or seq.startswith(motif_rc, i):
                    positions.append(i)

            if positions:
                reads_with_motif += 1
                last_end = 0

                for i, start in enumerate(positions):
                    frag_seq = seq[last_end:start + motif_len]
                    qual = record.letter_annotations["phred_quality"][last_end:start + motif_len]

                    frag = SeqRecord(
                        Seq(frag_seq),
                        id=f"{record.id}_frag{i+1}",
                        description="",
                        letter_annotations={"phred_quality": qual}
                    )
                    total_fragments.append(frag)
                    last_end = start + motif_len

                # Add final fragment if needed
                if last_end < len(seq):
                    frag_seq = seq[last_end:]
                    qual = record.letter_annotations["phred_quality"][last_end:]

                    frag = SeqRecord(
                        Seq(frag_seq),
                        id=f"{record.id}_frag{len(positions)+1}",
                        description="",
                        letter_annotations={"phred_quality": qual}
                    )
                    total_fragments.append(frag)

            else:
                # Keep whole read
                total_fragments.append(record)

    with gzip.open(output_fastq, "wt") if str(output_fastq).endswith(".gz") else open(output_fastq, "w") as out:
        SeqIO.write(total_fragments, out, "fastq")

    print(f"📊 Fragmentation done: {total_reads} reads processed")
    print(f"🔍 Reads with motif: {reads_with_motif}")
    print(f"📄 Total fragments generated: {len(total_fragments)}\n")


########################################
# Size filtering
########################################

def filter_viral_reads_by_size(input_fastq, output_fastq, min_len=100, max_len=3000):
    print(f"🧪 Filtering viral reads by length ({min_len}–{max_len} nt)...")

    total = 0
    kept = 0
    selected = []

    with open_fastq(input_fastq) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            total += 1
            if min_len <= len(record.seq) <= max_len:
                selected.append(record)
                kept += 1

    with open(output_fastq, "w") as out:
        SeqIO.write(selected, out, "fastq")

    print(f"📊 Total reads: {total}")
    print(f"✅ Kept for assembly: {kept}")
    print(f"🗑️ Discarded: {total - kept}\n")


########################################
# Normalize IDs (required by Flye)
########################################

def normalize_read_ids(input_fastq, output_fastq):
    from itertools import count
    c = count(1)

    with open_fastq(input_fastq) as infile, open(output_fastq, "w") as outfile:
        for rec in SeqIO.parse(infile, "fastq"):
            new_id = f"read{next(c):06d}"
            rec.id = new_id
            rec.name = new_id
            rec.description = ""
            SeqIO.write(rec, outfile, "fastq")


########################################
# Recentrifuge
########################################

def run_recentrifuge(kraken_output_files, kraken_db_dir, output_html):
    if output_html.exists() and output_html.stat().st_size > 0:
        print(f"⏭️ Recentrifuge already exists: {output_html.name}")
        return

    print("🧪 Running Recentrifuge...")

    k_args = []
    for f in kraken_output_files:
        k_args.extend(["-k", str(f)])

    cmd = [
        "rcf"
    ] + k_args + [
        "-n", str(kraken_db_dir),
        "-o", str(output_html),
        "-m", "50"
    ]

    run_command(cmd, "Recentrifuge", suppress_output=True)


########################################
# MAIN
########################################

def main():
    print_logo()
    check_dependencies()

    from argparse import RawTextHelpFormatter

    parser = argparse.ArgumentParser(
        description="📜 ViroCatch v1.1 - Pipeline for viral DNA detection from HTS data",
        add_help=True,
        usage="",
        formatter_class=lambda prog: RawTextHelpFormatter(prog, max_help_position=40)
    )

    group = parser.add_argument_group("🔧 Arguments")
    group.add_argument("--read_type", choices=["long", "short"], required=True)
    group.add_argument("--input_ONT", type=Path, help="ONT FASTQ")
    group.add_argument("--input_Illumina_R1", type=Path)
    group.add_argument("--input_Illumina_R2", type=Path)
    group.add_argument("--outdir", required=True, type=Path)
    group.add_argument("--threads", type=int, default=8)
    group.add_argument("--db", required=True, type=Path)
    group.add_argument("--kraken_db", required=True, type=Path)
    group.add_argument("--memory_mapping", action="store_true")
    group.add_argument("--quality_threshold", type=int, default=20)
    group.add_argument("--motif", type=str, default="TAATATTA")
    group.add_argument("--rca_r_product", action="store_true")
    group.add_argument("--force_spades", action="store_true")
    group.add_argument("--genome_size", type=str, default="2.7k")

    args = parser.parse_args()

    global memory_mapping
    memory_mapping = args.memory_mapping

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    qc_and_trimmed_dir = outdir / "1_QC_and_Trimmed"
    kraken_dir = outdir / "2_Kraken_Classification"
    assembly_dir = outdir / "3_Assembly"
    blast_dir = outdir / "4_Blast_Results"

    for d in [qc_and_trimmed_dir, kraken_dir, assembly_dir, blast_dir]:
        d.mkdir(exist_ok=True)

    assembly_file = None
    ########################################################
    # SHORT READS (ILLUMINA)
    ########################################################

    if args.read_type == "short":
        r1 = args.input_Illumina_R1
        r2 = args.input_Illumina_R2
        base_name = r1.stem.split(".")[0]

        trimmed_r1 = qc_and_trimmed_dir / f"{base_name}_R1_trimmed.fastq.gz"
        trimmed_r2 = qc_and_trimmed_dir / f"{base_name}_R2_trimmed.fastq.gz"

        run_command([
            "fastp",
            "-i", str(r1),
            "-I", str(r2),
            "-o", str(trimmed_r1),
            "-O", str(trimmed_r2),
            "-q", str(args.quality_threshold),
            "-w", str(args.threads),
            "--html", str(qc_and_trimmed_dir / f"{base_name}_fastp_report.html"),
            "--json", str(qc_and_trimmed_dir / f"{base_name}_fastp_report.json")
        ], f"Quality control with fastp for {base_name}",
            check_existing_output=trimmed_r1)

        # Kraken2
        kraken_report_r1 = kraken_dir / f"{base_name}_R1_kraken2_report.txt"
        kraken_output_r1 = kraken_dir / f"{base_name}_R1_kraken2_output.txt"
        kraken_report_r2 = kraken_dir / f"{base_name}_R2_kraken2_report.txt"
        kraken_output_r2 = kraken_dir / f"{base_name}_R2_kraken2_output.txt"

        run_kraken2(trimmed_r1, kraken_report_r1, kraken_output_r1, args.kraken_db, args.threads)
        run_kraken2(trimmed_r2, kraken_report_r2, kraken_output_r2, args.kraken_db, args.threads)

        # Extract viral reads
        viral_r1 = kraken_dir / f"{base_name}_viral_R1.fastq"
        viral_r2 = kraken_dir / f"{base_name}_viral_R2.fastq"

        extract_kraken_reads(kraken_output_r1, trimmed_r1, viral_r1)
        extract_kraken_reads(kraken_output_r2, trimmed_r2, viral_r2)

        # Remove low complexity
        cleaned_r1 = kraken_dir / f"{base_name}_viral_R1_cleaned.fastq"
        cleaned_r2 = kraken_dir / f"{base_name}_viral_R2_cleaned.fastq"

        filter_low_complexity_reads(viral_r1, cleaned_r1)
        filter_low_complexity_reads(viral_r2, cleaned_r2)

        # Read sync for SPAdes
        synced_r1 = assembly_dir / f"{base_name}_R1_synced.fastq"
        synced_r2 = assembly_dir / f"{base_name}_R2_synced.fastq"

        run_command([
            "python", extract_kraken_script_path,
            "-r1", str(cleaned_r1),
            "-r2", str(cleaned_r2),
            "-o1", str(synced_r1),
            "-o2", str(synced_r2),
            "--sync-pairs"
        ], "Synchronize paired-end reads",
            check_existing_output=synced_r1)

        # SPAdes
        assembly_fasta = assembly_dir / "assembly.fasta"

        if assembly_fasta.exists() and assembly_fasta.stat().st_size > 0:
            print(f"⏭️ SPAdes assembly already exists: {assembly_fasta.name}")
        else:
            spades_cmd = [
                "python", str(shutil.which("spades.py")),
                "-1", str(synced_r1),
                "-2", str(synced_r2),
                "-o", str(assembly_dir),
                "--meta",
                "--threads", str(args.threads)
            ]
            run_command(spades_cmd, "Assembly with SPAdes (Illumina)")

            scaffolds_file = assembly_dir / "scaffolds.fasta"
            contigs_file = assembly_dir / "contigs.fasta"

            if scaffolds_file.exists() and scaffolds_file.stat().st_size > 0:
                scaffolds_file.rename(assembly_fasta)
            elif contigs_file.exists() and contigs_file.stat().st_size > 0:
                contigs_file.rename(assembly_fasta)
            else:
                print("❌ SPAdes did not generate an assembly.")

        # Recentrifuge
        rec_html = kraken_dir / f"{base_name}_recentrifuge.html"
        run_recentrifuge([kraken_output_r1, kraken_output_r2], args.kraken_db, rec_html)

    ########################################################
    # LONG READS (ONT)
    ########################################################

    elif args.read_type == "long":

        input_ont = args.input_ONT
        base_name = input_ont.stem.split(".")[0]

        trimmed_fastq = qc_and_trimmed_dir / f"{base_name}_trimmed.fastq.gz"
        fastplong_html = qc_and_trimmed_dir / f"{base_name}_fastplong_report.html"
        fastplong_json = qc_and_trimmed_dir / f"{base_name}_fastplong_report.json"

        run_command([
            "fastplong",
            "-i", str(input_ont),
            "-o", str(trimmed_fastq),
            "-q", str(args.quality_threshold),
            "-w", str(args.threads),
            "--html", str(fastplong_html),
            "--json", str(fastplong_json)
        ], f"Quality control with fastplong for {base_name}",
            check_existing_output=trimmed_fastq)

        # RCA-R fragmentation
        if args.rca_r_product:
            fragments_fastq = qc_and_trimmed_dir / f"{base_name}_fragments.fastq"
            if not fragments_fastq.exists():
                fragment_by_motif(trimmed_fastq, fragments_fastq, args.motif)
            else:
                print(f"⏭️ Fragmentation already exists: {fragments_fastq}")
            input_for_kraken = fragments_fastq
        else:
            input_for_kraken = trimmed_fastq

        # Kraken
        kraken_report = kraken_dir / f"{base_name}_kraken2_report.txt"
        kraken_output = kraken_dir / f"{base_name}_kraken2_output.txt"

        run_kraken2(input_for_kraken, kraken_report, kraken_output, args.kraken_db, args.threads)

        # Extract viral reads
        viral_fastq = kraken_dir / f"{base_name}_viral.fastq"
        extract_kraken_reads(kraken_output, input_for_kraken, viral_fastq)

        # Low complexity removal
        cleaned_viral_fastq = kraken_dir / f"{base_name}_viral_cleaned.fastq"
        filter_low_complexity_reads(viral_fastq, cleaned_viral_fastq)

        # Size filtering (not for RCA-R)
        filtered_for_flye = cleaned_viral_fastq
        if not args.rca_r_product:
            filtered_for_flye = qc_and_trimmed_dir / f"{base_name}_viral_size_filtered.fastq"
            filter_viral_reads_by_size(cleaned_viral_fastq, filtered_for_flye, 500, 3000)

        # Normalize IDs (required by Flye)
        normalized_for_flye = filtered_for_flye.parent / f"{filtered_for_flye.stem}_normalized.fastq"
        normalize_read_ids(filtered_for_flye, normalized_for_flye)

        ###############################
        # SPAdes forced mode
        ###############################
        if args.force_spades:
            print("⚠️ SPAdes forced by user. Flye skipped.")

            spades_ready = assembly_dir / "viral_reads_clamped.fastq"
            print("🔧 Clamping qualities to Q40 for SPAdes...")
            clamp_qualities(cleaned_viral_fastq, spades_ready)

            spades_cmd = [
                "python", str(shutil.which("spades.py")),
                "--nanopore", str(spades_ready),
                "--threads", str(args.threads),
                "-o", str(assembly_dir)
            ]
            run_command(spades_cmd, "Assembly with SPAdes (ONT, clamped)")

            assembly_candidate = assembly_dir / "contigs.fasta"
            if assembly_candidate.exists():
                assembly_file = assembly_candidate
            else:
                print("❌ SPAdes did not produce contigs.")

        ###############################
        # Flye mode (default)
        ###############################
        else:
            flye_output_dir = assembly_dir / "flye_output"
            flye_output_dir.mkdir(exist_ok=True)

            flye_cmd = [
                "flye",
                "--nano-raw", str(normalized_for_flye),
                "--out-dir", str(flye_output_dir),
                "--threads", str(args.threads),
                "--genome-size", str(args.genome_size),
                "--meta"
            ]
            run_command(flye_cmd, "Assembly with Flye",
                        check_existing_output=flye_output_dir / "assembly.fasta")

            flye_assembly = flye_output_dir / "assembly.fasta"
            if flye_assembly.exists():
                assembly_file = assembly_dir / "assembly.fasta"
                shutil.copy(flye_assembly, assembly_file)
                print(f"🧬 Flye assembly copied to {assembly_file}")
            else:
                print("❌ Flye failed to assemble.")

        # Recentrifuge
        rec_html = kraken_dir / f"{base_name}_recentrifuge.html"
        run_recentrifuge([kraken_output], args.kraken_db, rec_html)

    ########################################################
    # BLASTn FINAL
    ########################################################

    if assembly_file and assembly_file.exists():
        blast_out = blast_dir / "blastn_results.tsv"

        blast_cmd = [
            "blastn",
            "-query", str(assembly_file),
            "-db", str(args.db),
            "-out", str(blast_out),
            "-num_threads", str(args.threads),
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                       "qstart qend sstart send evalue bitscore stitle",
            "-max_target_seqs", "1",
            "-max_hsps", "1"
        ]

        run_command(blast_cmd, "BLASTn analysis", check_existing_output=blast_out)

        sorted_blast = blast_dir / "blastn_results.sorted.tsv"
        with open(sorted_blast, "w") as f:
            subprocess.run(["sort", "-k3,3nr", str(blast_out)], stdout=f, check=True)

        shutil.move(str(sorted_blast), str(blast_out))

        print(f"🔬 BLASTn results saved: {blast_out}")


########################################################
# ENTRY POINT
########################################################

if __name__ == "__main__":
    extract_kraken_script_path = "extract_kraken_reads.py"
    main()