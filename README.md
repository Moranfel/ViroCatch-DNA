# 🦠 BegomoHunter v1.0 - HTS Begomovirus Detection Pipeline

**BegomoHunter** es un pipeline bioinformático multietapa para la detección y ensamblaje de virus **Begomovirus** a partir de datos de secuenciación **Illumina**, **Oxford Nanopore (ONT)** y **RCA-R**.

> 🧬 Versión 1.0 - Desarrollado por Moran et al., 2025

![BegomoHunter Logo](https://img.shields.io/badge/status-STABLE-green.svg)
![Python](https://img.shields.io/badge/python-3.6%2B-blue.svg)
![License](https://img.shields.io/badge/license-MIT-lightgrey.svg)

---

## 🚀 Características

- Control de calidad con [`fastp`](https://github.com/OpenGene/fastp) y `fastplong`
- Clasificación taxonómica con [`Kraken2`](https://ccb.jhu.edu/software/kraken2/)
- Filtro por complejidad y tamaño
- Fragmentación por motivo viral (RCA-R compatible)
- Ensamblaje con [`SPAdes`](https://cab.spbu.ru/software/spades/) o [`Flye`](https://github.com/fenderglass/Flye)
- Análisis final con [`BLASTn`](https://blast.ncbi.nlm.nih.gov/) y [`Recentrifuge`](https://github.com/khyox/recentrifuge)

---

## 🧰 Requisitos

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

### 📦 Python

- Python 3.6 o superior
- Módulos requeridos:
  - `biopython`

Instalación:
```bash
pip install biopython
