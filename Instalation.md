# ViroCatch-DNA v1.0

🦠 **ViroCatch-DNA** is a bioinformatics pipeline for detecting and assembling **Begomovirus** genomes from Illumina, Oxford Nanopore (ONT), and RCA-R amplicon sequencing data.

---

## 📋 Repository Contents

- `environment.yml` – Conda environment definition with all required dependencies.  
- `Begomo_Hunter.py` – Main pipeline script.

---

## ⚙️ Installation

1. **Clone the repository**  
   ```bash
   git clone https://github.com/Moranfel/Begomo_Hunter.git
   cd ViroCatch-DNA
2. conda env create -f ViroCatch-DNA_environment.yml
3. conda activate ViroCatch-DNA
4. python Begomo_Hunter.py --help

