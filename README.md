# Building biocontrol backwards
This repository contains scripts and supplementary data for the manuscript titled: Repository for manuscript 'Building biocontrol backwards for forward-thinking plant health management based on efficient mechanisms"

## Contents

- `data/` - Supplementary data

## Overview

This study employed a mechanism-first approach to biocontrol discovery, leveraging functional metagenomics to identify Actinobacterial stearoyl-CoA desaturase genes as ecological indicators of *Xylella fastidiosa* risk. The analysis involved two primary stages: (1) searching metagenomes for stearoyl-CoA desaturase genes, and (2) validating this association through qPCR on plant samples.

---

## Programs and Dependencies

The following bioinformatics tools and databases were used for this analysis:

### Taxonomic Classification
- **Kraken2** (v2.1.2): Taxonomic classification of metagenomic reads  
  GitHub: https://github.com/DerrickWood/kraken2  
  Databases: https://benlangmead.github.io/aws-indexes/k2

### Functional Gene Search
- **Kaiju**: Fast metagenomic classification for protein-coding genes  
  GitHub: https://github.com/bioinformatics-centre/kaiju

### Reference Databases
- **UniProt**: Protein sequence database  
  Website: https://www.uniprot.org/

---

## Analysis Workflow

### 1. Searching for Potential Quorum Signal Surrogates

#### 1.1 Data Download
SRA data were downloaded for all 59 samples listed in `data/MetagenomeList-Xf.csv`. Each sample entry contains plant host and NCBI SRA run accession (SRR).

#### 1.2 Taxonomic Screening with Kraken2

To detect *Xylella fastidiosa* in metagenomic samples, we classified all reads using Kraken2 with the PlusPF database:

```bash
cd /kraken2/kraken2-2.1.2/

for i in $(ls <input_dir>*fastq.gz | grep "_R1" | cut -f 1 -d "_"); do 
    N=$(basename $i .fastq.gz)
    ./kraken2 --db PlusPF \
        --gzip-compressed \
        --paired \
        --report <output_dir>/$N.kreport \
        --output <output_dir>/$N.kraken \
        ${i}_R1.fastq.gz ${i}_R2.fastq.gz
done
```

**Parameters:**
- `--db PlusPF`: Kraken2 PlusPF database (bacteria, archaea, viruses, fungi, plant)
- `--paired`: Paired-end read mode
- `--report`: Classification report with abundance estimates
- `--output`: Detailed classification per read

#### 1.3 Functional Gene Search with Kaiju

To identify stearoyl-CoA desaturase genes from bacteria and fungi, we created a custom Kaiju database from UniProt sequences (22,000 sequences matching query: "stearoyl-CoA desaturase").

**Database Construction:**

```bash
cd kaiju/bin

# Create index
./kaiju-mkbwt -n 5 \
    -a ACDEFGHIKLMNPQRSTVWY \
    -o <output_dir> \
    <fasta_path>

# Create FM-index for fast search
./kaiju-mkfmi <output_dir>
```

**Functional Classification:**

```bash
cd kaiju/bin

for i in $(ls <raw_data_dir>*.fastq | grep "_R1" | cut -f 1 -d "_"); do
    N=$(basename $i .fastq)
    ./kaiju \
        -t /path/to/nodes.dmp \
        -f <output_dir>.fmi \
        -i ${i}_R1.fastq \
        -j ${i}_R2.fastq \
        -o <kaiju_output_dir>/$N.out \
        -v \
        -s 60 \
        -e 10 \
        -m 20
done
```

**Parameters:**
- `-t`: NCBI nodes.dmp file for taxonomy
- `-f`: FM-index path
- `-i, -j`: Forward and reverse read files
- `-o`: Output classification file
- `-v`: Verbose output
- `-s`: Sensitivity (60 = high sensitivity)
- `-e`: Maximum number of E-values to report
- `-m`: Minimum match length (amino acids)

---

### 2. Validating Actinobacteria Stearoyl-CoA Desaturase as an Xf Risk Indicator

#### 2.1 qPCR Assay Design and Validation

The following primers and probes were used for quantitative real-time PCR (qPCR): `data/qPCR_primers.csv`

**Final reaction volume:** 10 µL

| Component | Concentration | Volume |
|-----------|---------------|--------|
| SsoAdvanced Universal SYBR Green Supermix | 5 µL | 5.0 µL |
| Forward Primer | 300 nM | variable |
| Reverse Primer | 300 nM | variable |
| Template DNA | 0.04 ng | 2.0 µL |
| Nuclease-free Water | - | 2.4 µL |

**Cycling Conditions:**
- Denaturation: 95°C, 30 s
- Annealing/Extension: 60°C, 30 s
- Plate read and melt curve
- 40 cycles total

**Technical Replicates:** Three per sample

#### 2.3 Data Analysis

Cq values and melt curve temperatures were extracted from qPCR runs. Absolute quantification was performed using:
- Standard curves generated from serial dilutions
- Conversion to log₁₀ absolute copy number
- Correlation analysis between *Xf* and Actinobacterial stearoyl-CoA desaturase abundance

---

## File Structure

```
xf-microbiome-biocontrol/
├── README.md                          # This file
├── LICENSE                            # Repository license
├── data/
│   ├── MetagenomeList-Xf.xlsx        # Sample metadata (59 metagenomes)
│   ├── srr_sample_mapping.csv        # SRA accession to sample mapping
│   └── supplemental_samples.csv       # Additional qPCR samples
├── scripts/
│   ├── 01_kraken2_classification.sh   # Kraken2 taxonomic profiling
│   ├── 02_kaiju_functional_search.sh  # Kaiju gene search
│   └── 03_qpcr_analysis.R             # qPCR data analysis
├── docs/
│   ├── METHODS.md                     # Detailed methods
│   ├── DATABASE_CONSTRUCTION.md       # Custom database creation
│   └── PRIMER_DESIGN.md               # Primer development notes
└── results/
    ├── kraken2/                       # Kraken2 classification reports
    ├── kaiju/                         # Kaiju functional annotations
    └── qpcr/                          # qPCR quantification data
```

---

## Quick Start

### Prerequisites
- Kraken2 v2.1.2 or later
- Kaiju
- R (v3.6+) with packages: ggplot2, tidyverse, Hmisc
- SRA Toolkit for downloading fastq files from NCBI

### Running the Analysis

1. **Download metagenomes** from NCBI SRA using accessions in `data/MetagenomeList-Xf.xlsx`

2. **Run taxonomic classification:**
   ```bash
   bash scripts/01_kraken2_classification.sh
   ```

3. **Search for stearoyl-CoA desaturase genes:**
   ```bash
   bash scripts/02_kaiju_functional_search.sh
   ```

4. **Analyze qPCR data:**
   ```R
   source("scripts/03_qpcr_analysis.R")
   ```

---

## Metadata

### Sample Information
- **Total metagenomes analyzed:** 59 (publicly available)
- **Supplemental qPCR samples:** 61
- **qPCR targets detected in both assays:** 37 samples
- **Undetectable *Xf*:** 14 samples
- **Undetectable stearoyl-CoA desaturase:** 16 samples
