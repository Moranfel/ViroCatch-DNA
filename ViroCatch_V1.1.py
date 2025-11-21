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
🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬 V.1.0 by Morán et al. 2025 🦠🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬🦠🦠🦠🧬
"""
    print(logo)
    print(" 🧬 Detection and identification of DNA viruses from ONT, RCA-R-ONT, and Illumina data 🧬\n")

def check_dependencies():
    print("🔍 Checking dependencies...")
    required_python_modules = ["Bio"]
    required_commands = ["fastp", "fastplong", "kraken2", "flye", "python", "blastn", "awk", "spades.py"]
    for module in required_python_modules:
        if importlib.util.find_spec(module) is None:
            print(f"❌ Missing module: {module}. Instálalo con pip.")
            sys.exit(1)
    for cmd in required_commands:
        if not shutil.which(cmd):
            print(f"❌ Missing command in PATH: {cmd}")
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
    if check_existing_output and Path(check_existing_output).exists():
        print(f"✅ {description} Already exists, skipping.")
        return
    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        if not suppress_output:
            print(result.stdout)
            print(result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"❌ Error while running: {e.cmd}\n{e.stderr}")
        sys.exit(1)

def run_kraken2(input_fastq, output_report, output_output, db_path, threads):
    global memory_mapping
    if output_report.exists() and output_report.stat().st_size > 0 and \
       output_output.exists() and output_output.stat().st_size > 0:
        print(f"✅ Kraken2 classification already exists for {input_fastq.name}, Skipping.")
        return

    kraken2_command = [
        "kraken2", "--db", str(db_path),
        "--report", str(output_report),
        "--output", str(output_output),
        "--threads", str(threads),
        str(input_fastq)
    ]

    if memory_mapping:
        kraken2_command.insert(1, "--memory-mapping")

    run_command(kraken2_command, f"Kraken2 para {input_fastq.name}", check_existing_output=output_report)

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
    run_command(cmd, f"Excluding Eukaryotic 🌱 and Prokaryotic 🧫 reads", suppress_output=True, check_existing_output=output_fastq)

def filter_viral_fragments_by_motif_and_size(input_fastq, output_fastq, motif="TAATATTA"):
    print(f"🧪 Filtering viral reads by size (500-3000bp) and motif '{motif}' or RC...")
    motif = motif.upper()
    motif_rc = str(Seq(motif).reverse_complement())
    kept = 0

    selected = []
    with open_fastq(input_fastq) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            seq = str(record.seq).upper()
            if (motif in seq or motif_rc in seq) and (500 <= len(seq) <= 3000):
                selected.append(record)
                kept += 1

    with open(output_fastq, "w") as out:
        SeqIO.write(selected, out, "fastq")

    print(f"✅ Selected viral reads: {kept}")

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
    print(f"✅ Conserved reads for Flye: {kept}")
    print(f"🗑️ Discarded reads: {total - kept}\n")

def open_fastq(filepath):
    return gzip.open(filepath, "rt") if str(filepath).endswith(".gz") else open(filepath, "r")

def is_low_complexity(seq, entropy_threshold=1.8):
    freqs = [seq.count(n)/len(seq) for n in "ATCG"]
    entropy = -sum(f * math.log2(f) for f in freqs if f > 0)
    return entropy < entropy_threshold

    # Código muerto mantenido sin modificar para respetar la opción A
    freqs = [seq.count(n)/len(seq) for n in "ATCG"]
    entropy = -sum(f * math.log2(f) for f in freqs if f > 0)
    return repetitive or entropy < entropy_threshold

def filter_low_complexity_reads(input_fastq, output_fastq):
    total = 0
    kept = 0
    discarded = 0

    with open_fastq(input_fastq) as infile, open(output_fastq, "w") as outfile:
        for record in SeqIO.parse(infile, "fastq"):
            total += 1
            seq = str(record.seq).upper()
            if not is_low_complexity(seq):
                SeqIO.write(record, outfile, "fastq")
                kept += 1
            else:
                discarded += 1

    print(f"🧼 Filtered by low complexity: {kept} conserved reads, {discarded} discarded from {total} total.")

def fragment_by_motif(input_fastq, output_fastq, motif="TAATATTA"):
    print(f"✂️  Fragmenting by motif '{motif}' and its reverse complement...")
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
                    frag = SeqRecord(
                        Seq(frag_seq),
                        id=f"{record.id}_frag{i + 1}",
                        description="",
                        letter_annotations={"phred_quality": record.letter_annotations["phred_quality"][last_end:start + motif_len]}
                    )
                    total_fragments.append(frag)
                    last_end = start + motif_len

                if last_end < len(seq):
                    frag_seq = seq[last_end:]
                    frag = SeqRecord(
                        Seq(frag_seq),
                        id=f"{record.id}_frag{len(positions)+1}",
                        description="",
                        letter_annotations={"phred_quality": record.letter_annotations["phred_quality"][last_end:]}
                    )
                    total_fragments.append(frag)

            else:
                total_fragments.append(record)

    with gzip.open(output_fastq, "wt") if str(output_fastq).endswith(".gz") else open(output_fastq, "w") as out:
        SeqIO.write(total_fragments, out, "fastq")

    print("📊 Fragmentation completed without size filtering.")
    print(f"   🧬 Processed reads: {total_reads}")
    print(f"   🔍 Reads with motif: {reads_with_motif}")
    print(f"   📄 Total fragments generated: {len(total_fragments)}\n")

def normalize_read_ids(input_fastq, output_fastq):
    from itertools import count
    counter = count(1)
    with open_fastq(input_fastq) as infile, open(output_fastq, "w") as outfile:
        for record in SeqIO.parse(infile, "fastq"):
            record.id = f"read{next(counter):06d}"
            record.name = record.id
            record.description = ""
            SeqIO.write(record, outfile, "fastq")

def synchronize_paired_reads(r1_path, r2_path, synced_r1_path, synced_r2_path):
    print("🔁 Synchronizing paired reads between R1 and R2...")

    r1_path_str = str(r1_path)
    r2_path_str = str(r2_path)
    synced_r1_path_str = str(synced_r1_path)
    synced_r2_path_str = str(synced_r2_path)

    if not Path(r1_path_str).is_file() or not Path(r2_path_str).is_file():
        print("❌ Input files R1 or R2 do not exist.")
        sys.exit(1)

    try:
        r1_reads = SeqIO.index(r1_path_str, "fastq")
        r2_reads = SeqIO.index(r2_path_str, "fastq")
    except Exception as e:
        print(f"❌ Error indexing reads: {e}")
        sys.exit(1)

    def normalize_keys(index):
        return {key.split()[0]: key for key in index.keys()}

    r1_keymap = normalize_keys(r1_reads)
    r2_keymap = normalize_keys(r2_reads)
    shared_ids = set(r1_keymap.keys()) & set(r2_keymap.keys())

    print(f"🧬 Shared reads found: {len(shared_ids)}")

    written = 0
    try:
        with open(synced_r1_path_str, "w") as r1_out, open(synced_r2_path_str, "w") as r2_out:
            for base_id in shared_ids:
                r1_record = r1_reads.get(r1_keymap[base_id])
                r2_record = r2_reads.get(r2_keymap[base_id])
                if r1_record and r2_record:
                    SeqIO.write(r1_record, r1_out, "fastq")
                    SeqIO.write(r2_record, r2_out, "fastq")
                    written += 1

        print(f"✅ {written} Synchronized pairs written.")
    except Exception as e:
        print(f"❌ Error writing synchronized files: {e}")
        sys.exit(1)
    finally:
        r1_reads.close()
        r2_reads.close()

def run_recentrifuge(kraken_output_files, kraken_db_dir, output_html):
    if output_html.exists() and output_html.stat().st_size > 0:
        print(f"✅ Recentrifuge already executed, existing file: {output_html.name}")
        return

    print("🧪 Running Recentrifuge on Kraken2 output files...")

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

def main():
    print_logo()
    check_dependencies()

    from argparse import RawTextHelpFormatter

    parser = argparse.ArgumentParser(
        description="📜 ViroCatch v1.0 - Pipeline for viral DNA detection from HTS data",
        add_help=True,
        usage="",
        formatter_class=lambda prog: RawTextHelpFormatter(prog, max_help_position=40)
    )

    group = parser.add_argument_group("🔧 Arguments")
    group.add_argument("--read_type", choices=["long", "short"], required=True, help="Tipo de lectura: 'short' (Illumina) o 'long' (ONT)")
    group.add_argument("--input_ONT", type=Path, metavar="", help="Archivo FASTQ de ONT (lecturas largas)")
    group.add_argument("--input_Illumina_R1", type=Path, metavar="", help="FASTQ R1 de Illumina")
    group.add_argument("--input_Illumina_R2", type=Path, metavar="", help="FASTQ R2 de Illumina")
    group.add_argument("--outdir", required=True, type=Path, metavar="", help="Kraken2 database directory")
    group.add_argument("--threads", type=int, default=8, help="Número de hilos de procesamiento")
    group.add_argument("--db", required=True, type=Path, metavar="", help="Base de datos BLASTn (preformateada con makeblastdb)")
    group.add_argument("--kraken_db", required=True, type=Path, metavar="", help="Directorio de base de datos Kraken2- Download from https://benlangmead.github.io/aws-indexes/k2")
    group.add_argument("--memory_mapping", action="store_true", help="Activar --memory-mapping en Kraken2")
    group.add_argument("--quality_threshold", type=int, metavar="", default=20, help="Umbral de calidad mínima para trimming")
    group.add_argument("--motif", type=str, default="TAATATTA", help="Motivo para fragmentación RCA-R y filtrado - Por defecto TAATATTA ")
    group.add_argument("--rca_r_product", action="store_true", help="Activar procesamiento RCA-R para ONT")
    group.add_argument("--force_spades", action="store_true", help="Forzar ensamblaje con SPAdes incluso con lecturas largas")
    group.add_argument("--genome_size", type=str, default="2.7k", help="Tamaño estimado del genoma para Flye (ej. 2k, 10k, 150k, 1m)")

    args = parser.parse_args()
    global memory_mapping
    memory_mapping = args.memory_mapping

    if args.read_type == "short":
        if not args.input_Illumina_R1 or not args.input_Illumina_R2:
            print("❌ --input_Illumina_R1 y R2 are required with --read_type short.")
            sys.exit(1)
        if args.rca_r_product:
            print("❌ --rca_r_product can only be used with --read_type long.")
            sys.exit(1)
    elif args.read_type == "long" and not args.input_ONT:
        print("❌ --input_ONT is required with --read_type long.")
        sys.exit(1)

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    qc_and_trimmed_dir = outdir / "1_QC_and_Trimmed"
    kraken_dir = outdir / "2_Kraken_Classification"
    assembly_dir = outdir / "3_Assembly"
    blast_dir = outdir / "4_Blast_Results"

    for d in [qc_and_trimmed_dir, kraken_dir, assembly_dir, blast_dir]:
        d.mkdir(exist_ok=True)

    assembly_file = None

    if args.read_type == "short":

        input_r1 = args.input_Illumina_R1
        input_r2 = args.input_Illumina_R2
        base_name = input_r1.name.split('_R1')[0]

        trimmed_r1 = qc_and_trimmed_dir / f"{base_name}_R1_trimmed.fastq.gz"
        trimmed_r2 = qc_and_trimmed_dir / f"{base_name}_R2_trimmed.fastq.gz"
        fastp_report_html = qc_and_trimmed_dir / f"{base_name}_fastp_report.html"
        fastp_report_json = qc_and_trimmed_dir / f"{base_name}_fastp_report.json"

        fastp_cmd = [
            "fastp", "-i", str(input_r1), "-o", str(trimmed_r1),
            "-I", str(input_r2), "-O", str(trimmed_r2),
            "-q", str(args.quality_threshold), "-w", str(args.threads),
            "-h", str(fastp_report_html), "-j", str(fastp_report_json),
            "--complexity_threshold", "30"
        ]

        run_command(fastp_cmd, f"Control de calidad y recorte con fastp para {base_name}", check_existing_output=trimmed_r1)

        kraken_report_r1 = kraken_dir / f"{base_name}_R1_kraken2_report.txt"
        kraken_output_r1 = kraken_dir / f"{base_name}_R1_kraken2_output.txt"
        run_kraken2(trimmed_r1, kraken_report_r1, kraken_output_r1, args.kraken_db, args.threads)

        kraken_report_r2 = kraken_dir / f"{base_name}_R2_kraken2_report.txt"
        kraken_output_r2 = kraken_dir / f"{base_name}_R2_kraken2_output.txt"
        run_kraken2(trimmed_r2, kraken_report_r2, kraken_output_r2, args.kraken_db, args.threads)

        viral_fastq_r1 = kraken_dir / f"{base_name}_R1_viral_reads.fastq"
        extract_kraken_reads(kraken_output_r1, trimmed_r1, viral_fastq_r1)

        cleaned_fastq_r1 = kraken_dir / f"{base_name}_R1_viral_cleaned.fastq"
        filter_low_complexity_reads(viral_fastq_r1, cleaned_fastq_r1)

        viral_fastq_r2 = kraken_dir / f"{base_name}_R2_viral_reads.fastq"
        extract_kraken_reads(kraken_output_r2, trimmed_r2, viral_fastq_r2)

        cleaned_fastq_r2 = kraken_dir / f"{base_name}_R2_viral_cleaned.fastq"
        filter_low_complexity_reads(viral_fastq_r2, cleaned_fastq_r2)

        synced_r1 = kraken_dir / f"{base_name}_R1_viral_reads_synced.fastq"
        synced_r2 = kraken_dir / f"{base_name}_R2_viral_reads_synced.fastq"
        synchronize_paired_reads(cleaned_fastq_r1, cleaned_fastq_r2, synced_r1, synced_r2)

        if not synced_r1.is_file() or synced_r1.stat().st_size == 0 or \
           not synced_r2.is_file() or synced_r2.stat().st_size == 0:
            print(f"❌ Error: Synchronized files not found or empty. Skipping assembly.")
        else:
            assembly_fasta = assembly_dir / "assembly.fasta"
            if assembly_fasta.exists() and assembly_fasta.stat().st_size > 0:
                print(f"✅ SPAdes assembly already exists: {assembly_fasta.name}, Skipping.")
                assembly_file = assembly_fasta
            else:
                spades_cmd = [
                    "python", str(shutil.which("spades.py")),
                    "-1", str(synced_r1), "-2", str(synced_r2),
                    "-o", str(assembly_dir),
                    "--meta", "--threads", str(args.threads)
                ]
                run_command(spades_cmd, f"Assembly with SPAdes metaviral")

                scaffolds_file = assembly_dir / "scaffolds.fasta"
                contigs_file = assembly_dir / "contigs.fasta"

                if scaffolds_file.exists() and scaffolds_file.stat().st_size > 0:
                    scaffolds_file.rename(assembly_fasta)
                    assembly_file = assembly_fasta
                elif contigs_file.exists() and contigs_file.stat().st_size > 0:
                    contigs_file.rename(assembly_fasta)
                    assembly_file = assembly_fasta
                else:
                    print("⚠️ No valid assembly file found.")

        recentrifuge_html = kraken_dir / f"{base_name}_recentrifuge.html"
        run_recentrifuge([kraken_output_r1, kraken_output_r2], args.kraken_db, recentrifuge_html)

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
        ], f"Quality control with fastplong for {base_name}", check_existing_output=trimmed_fastq)

        if args.rca_r_product:
            fragments_fastq = qc_and_trimmed_dir / f"{base_name}_fragments.fastq"
            if not fragments_fastq.exists() or fragments_fastq.stat().st_size == 0:
                fragment_by_motif(trimmed_fastq, fragments_fastq, args.motif)
            else:
                print(f"✅ Fragmentation already exists: {fragments_fastq.name}, skipping.")

            input_for_kraken = fragments_fastq
        else:
            input_for_kraken = trimmed_fastq

        kraken_report = kraken_dir / f"{base_name}_kraken2_report.txt"
        kraken_output = kraken_dir / f"{base_name}_kraken2_output.txt"
        run_kraken2(input_for_kraken, kraken_report, kraken_output, args.kraken_db, args.threads)

        viral_fastq = kraken_dir / f"{base_name}_viral_reads.fastq"
        extract_kraken_reads(kraken_output, input_for_kraken, viral_fastq)

        cleaned_viral_fastq = kraken_dir / f"{base_name}_viral_cleaned.fastq"
        filter_low_complexity_reads(viral_fastq, cleaned_viral_fastq)
        filtered_for_flye = cleaned_viral_fastq

        print(f"✅ Will be used directly '{filtered_for_flye.name}' for assembly with Flye.")

        if not args.rca_r_product and args.motif == "TAATATTA":
            filtered_for_flye = qc_and_trimmed_dir / f"{base_name}_viral_filtered_size_and_complexity.fastq"
            filter_viral_reads_by_size(cleaned_viral_fastq, filtered_for_flye, min_len=500, max_len=3000)

        flye_output_dir = assembly_dir / "flye_output"
        flye_output_dir.mkdir(exist_ok=True)

        normalized_for_flye = filtered_for_flye.parent / (filtered_for_flye.stem + "_normalized.fastq")
        normalize_read_ids(filtered_for_flye, normalized_for_flye)

        if args.force_spades:
            print("⚠️ SPAdes forzado por el usuario. Flye será omitido.")
        else:
            flye_cmd = [
                "flye",
                "--nano-raw", str(normalized_for_flye),
                "--out-dir", str(flye_output_dir),
                "--threads", str(args.threads),
                "--genome-size", str(args.genome_size),
                "--min-overlap", "1000",
                "--meta"
            ]
            run_command(flye_cmd, "Ensamblaje con Flye", check_existing_output=flye_output_dir / "assembly.fasta")

        if args.force_spades:
            print("🔧 Ejecutando SPAdes (modo ONT) forzado por el usuario")
            print("⚠️  Flye será omitido automáticamente")

            spades_cmd = [
                "python", str(shutil.which("spades.py")),
                "--nanopore", str(cleaned_viral_fastq),
                "--threads", str(args.threads),
                "-o", str(assembly_dir)
            ]

            run_command(spades_cmd, "Assembly with SPAdes (ONT, --nanopore)")

            assembly_candidate = assembly_dir / "contigs.fasta"
            if assembly_candidate.exists() and assembly_candidate.stat().st_size > 0:
                assembly_file = assembly_candidate
            else:
                print("❌ SPAdes no generó contigs. Revisa warnings.log y spades.log")

        else:
            flye_cmd = [
                "flye",
                "--nano-raw", str(normalized_for_flye),
                "--out-dir", str(flye_output_dir),
                "--threads", str(args.threads),
                "--genome-size", str(args.genome_size),
                "--min-overlap", "1000",
                "--meta"
            ]
            run_command(flye_cmd, "Ensamblaje con Flye",
                        check_existing_output=flye_output_dir / "assembly.fasta")

            flye_assembly = flye_output_dir / "assembly.fasta"
            if flye_assembly.exists() and flye_assembly.stat().st_size > 0:
                assembly_file = assembly_dir / "assembly.fasta"
                shutil.copy(flye_assembly, assembly_file)
                print(f"🧬 Ensamblaje Flye copiado a {assembly_file}")
            else:
                print("❌ Flye no generó ensamblaje.")
        if not assembly_file or not assembly_file.exists():
            flye_assembly = flye_output_dir / "assembly.fasta"
            if flye_assembly.exists() and flye_assembly.stat().st_size > 0:
                final_assembly = assembly_dir / "assembly.fasta"
                shutil.copy(flye_assembly, final_assembly)
                assembly_file = final_assembly
                print(f"⚠️ Unfiltered assembly copied as fallback to {final_assembly}")
            else:
                print("⚠️ Flye assembly not generated.")

        recentrifuge_html = kraken_dir / f"{base_name}_recentrifuge.html"
        run_recentrifuge([kraken_output], args.kraken_db, recentrifuge_html)

    if assembly_file and assembly_file.exists():
        blast_out = blast_dir / "blastn_results.tsv"
        blast_cmd = [
            "blastn", "-query", str(assembly_file),
            "-db", str(args.db),
            "-out", str(blast_out),
            "-num_threads", str(args.threads),
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle",
            "-max_target_seqs", "1",
            "-max_hsps", "1"
        ]
        run_command(blast_cmd, "Análisis BLASTn", check_existing_output=blast_out)

        sorted_blast = blast_dir / "blastn_results.sorted.tsv"
        print("🔧 Running: Sorting BLASTN results by % identity")

        with open(sorted_blast, "w") as f:
            subprocess.run(["sort", "-k3,3nr", str(blast_out)], stdout=f, check=True)

        shutil.move(str(sorted_blast), str(blast_out))

if __name__ == "__main__":
    extract_kraken_script_path = "extract_kraken_reads.py"
    main()