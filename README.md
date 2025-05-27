# 🦠 BegomoHunter - HTS Begomovirus Detection Pipeline

**BegomoHunter is a multi‐stage bioinformatics Python pipeline for detecting and assembling Begomovirus genomes from Illumina, Oxford Nanopore (ONT), and RCA-R Amplicon sequencing data**.

> 🧬 Versión 1.0 - Developed by Félix Morán

![BegomoHunter Logo](https://img.shields.io/badge/status-STABLE-green.svg)
![Python](https://img.shields.io/badge/python-3.6%2B-blue.svg)
![License](https://img.shields.io/badge/license-MIT-lightgrey.svg)

---

## ⚙️ COMMANDS

python Begomo_Hunter.py \
  --read_type <long|short> \
  --input_ONT <ONT.fastq.gz> \
  --input_Illumina_R1 <R1.fastq.gz> \
  --input_Illumina_R2 <R2.fastq.gz> \
  --outdir <output_directory> \
  --threads \ By default: 8 
  --db <blast_db_path> \
  --kraken_db <kraken2_db_path> \
  [--quality_threshold] \ By default: 20
  [--motif] \ By default: TAATATTA
  [--rca_r_product]

Note: To avoid issues, do not use spaces in the input file names (--input_ONT, --input_Illumina_R1, --input_Illumina_R2) or in --outdir

---
## 🧬 DATABASES 🦠 

Parameter --db: Path to a pre‐built BLAST nucleotide database. You can create one from a reference FASTA file using makeblastdb.
Example:
  makeblastdb \
  -in begomovirus_ref.fasta \
  -dbtype nucl \
  -out blast_db/begomo_db \
  -parse_seqids

Parameter --kraken_db: Path to a Kraken2 database, built as described in [Kraken2](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown). It's recommended to use the collection [PlusPFP](https://benlangmead.github.io/aws-indexes/k2)

## 🚀 FEATURES

- Quality control with [`fastp`](https://github.com/OpenGene/fastp) for short reads and [`fastplong`](https://github.com/OpenGene/fastplong) for long reads
- Taxonomic classification with [`Kraken2`](https://ccb.jhu.edu/software/kraken2/)
- Complexity and size filtering for begomovirus
- Viral‐motif‐based fragmentation  (just for RCA-R sequences  [--rca_r_product])
- Assembly of short reads with [`SPAdes`](https://cab.spbu.ru/software/spades/) or long reads with [`Flye`](https://github.com/fenderglass/Flye)
- Final analysis with [`BLASTn`](https://blast.ncbi.nlm.nih.gov/) and  [`Recentrifuge`](https://github.com/khyox/recentrifuge)

---

## 🧰 EXTERNAL REQUIREMENTS & DEPENDENCIES

| Software       | Description                         |
|----------------|-------------------------------------|
| `fastp`        | QC - Illumina                       |
| `fastplong`    | QC - ONT                            |
| `kraken2`      | Taxonomy Clasification              |
| `blastn`       | Final analysis of contigs           |
| `flye`         | Assembler long reads                |
| `spades.py`    | Assembler short reads               |
| `rcf`          | View taxonomy Clasification         |
| `awk`, `python`, `gzip`, `sort` | Standard utilities |

### 📦 PYTHON

- Python 3.6 o superior
- Required modules:
  - `biopython`

Instalación:
```bash
pip install biopython
