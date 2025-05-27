# 🦠 BegomoHunter v1.0 - HTS Begomovirus Detection Pipeline

**BegomoHunter** es un pipeline bioinformático multietapa para la detección y ensamblaje de virus **Begomovirus** a partir de datos de secuenciación **Illumina**, **Oxford Nanopore (ONT)** y **Amplicones de RCA-R**.

> 🧬 Versión 1.0 - Desarrollado por Félix Morán

![BegomoHunter Logo](https://img.shields.io/badge/status-STABLE-green.svg)
![Python](https://img.shields.io/badge/python-3.6%2B-blue.svg)
![License](https://img.shields.io/badge/license-MIT-lightgrey.svg)

---

## ⚙️ COMANDOS

python Begomo_Hunter.py \
  --read_type <long|short> \
  --input_ONT <ONT.fastq.gz> \
  --input_Illumina_R1 <R1.fastq.gz> \
  --input_Illumina_R2 <R2.fastq.gz> \
  --outdir <output_directory> \
  --threads 8 \
  --db <blast_db_path> \
  --kraken_db <kraken2_db_path> \
  [--quality_threshold 20] \
  [--motif TAATATTA] \
  [--rca_r_product]

---
## 🧬 BASES DE DATOS 

Parametro --db: Este parámetro requiere una base de datos BLAST previamente construida. Puedes crear una a partir de un archivo FASTA de referencia usando makeblastdb.
ejemplo: 
  makeblastdb \
  -in begomovirus_ref.fasta \
  -dbtype nucl \
  -out blast_db/begomo_db \
  -parse_seqids

Parametro --kraken_db: Kraken2 database, se debe contruir como se menciona en ^*[Kraken2](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown) y se recomuenda usar la collection [PlusPFP](https://benlangmead.github.io/aws-indexes/k2)

## 🚀 CARACTERÍSTICAS

- Control de calidad con [`fastp`](https://github.com/OpenGene/fastp) y `fastplong`
- Clasificación taxonómica con [`Kraken2`](https://ccb.jhu.edu/software/kraken2/)
- Filtro por complejidad y tamaño
- Fragmentación por motivo viral (RCA-R compatible)
- Ensamblaje short reads [`SPAdes`](https://cab.spbu.ru/software/spades/) o para long reads [`Flye`](https://github.com/fenderglass/Flye)
- Análisis final con [`BLASTn`](https://blast.ncbi.nlm.nih.gov/) y [`Recentrifuge`](https://github.com/khyox/recentrifuge)

---

## 🧰 REQUISITOS-DEPENDENCIAS

### 🔧 Dependencias externas

| Software       | Descripción                |
|----------------|----------------------------|
| `fastp`        | QC de Illumina             |
| `fastplong`    | QC para ONT                |
| `kraken2`      | Clasificación taxonómica   |
| `blastn`       | Análisis final de ensamblaje |
| `flye`         | Ensamblador ONT            |
| `spades.py`    | Ensamblador Illumina       |
| `rcf`          | CLI de Recentrifuge        |
| `awk`, `python`, `gzip`, `sort` | Utilidades estándar |

### 📦 PYTHON

- Python 3.6 o superior
- Módulos requeridos:
  - `biopython`

Instalación:
```bash
pip install biopython
