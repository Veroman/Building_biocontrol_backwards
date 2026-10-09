# Building biocontrol backwards

Scripts, processed data, and supplementary tables for the manuscript:

> **Building biocontrol backwards for forward-thinking plant health management based on efficient mechanisms**
> Cindy E. Morris, Sara Campigli, Guido Marchi, and Veronica Roman-Reyna.

Archived release: scholarsphere

## Overview

This study used functional metagenomics to identify microbial genes associated with plant health, using *Xylella fastidiosa* (Xf) as a case study. We searched olive metagenomes for stearoyl-CoA (Δ9) desaturase genes, which we hypothesised would be negatively associated with Xf. The analysis had two stages:

1. screening Xf-related metagenomes for stearoyl-CoA desaturase genes, and
2. validating the relationship between Actinomycetota Δ9-desaturase abundance and Xf using qPCR.

## Repository contents
---

## Programs and Dependencies

The following bioinformatics tools and databases were used for this analysis:

### Taxonomic Classification
- **Kraken2:** Taxonomic classification of metagenomic reads  
  GitHub: https://github.com/DerrickWood/kraken2  
  Databases: https://benlangmead.github.io/aws-indexes/k2

### Functional Gene Search
- **Kaiju:** Metagenomic classification for protein-coding genes  
  GitHub: https://github.com/bioinformatics-centre/kaiju

### Reference Databases and samples
- **UniProt:** Protein sequence database  
  Website: https://www.uniprot.org/
- **NCBI SRA Toolkit:** extract FASTQ-files from SRA-accessions  
  GitHub: https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit

---

## Analysis Workflow

### 1. Searching for Potential Quorum Signal Surrogates

#### 1.1 Data Download
SRA data were downloaded, and the FASTQ files extracted using NCBI SRA Toolkit prefetch and fasterq-dump.
All metagenomic samples (raw data) is listed in `data/MetagenomeList-Xf.csv`. Each sample entry contains the plant host and the NCBI SRA run accession (SRR).

#### 1.2 Taxonomic Screening with Kraken2

To detect *Xylella fastidiosa* in metagenomic samples, we classified all reads using Kraken2 with the PlusPF database:

```bash
cd /kraken2/kraken2-2.1.2/

for i in $(ls <input_dir>/*fastq.gz | grep "_R1" | cut -f 1 -d "_"); do 
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
- `--db PlusPF`: Kraken2 PlusPF database (bacteria, archaea, viruses, protozoa, and fungi)
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
    -o <output_dir>/proteins \
    <fasta_path>

# Create FM-index for fast search
./kaiju-mkfmi <output_dir>/proteins
```

**Functional Classification:**

```bash
cd kaiju/bin

for i in $(ls <raw_data_dir>/*.fastq | grep "_R1" | cut -f 1 -d "_"); do
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
- initial denaturation: 98°C, 3 min
- Denaturation: 95°C, 15 s
- Annealing/Extension: 60°C, 40 s
- 40 cycles total
- Melt curve:
  -  65°C, 5 s
  -  95°C, 0.5°C/cy

**Technical Replicates:** Three per sample

#### 2.3 Data Analysis

Cq values and melt curve temperatures were extracted from qPCR runs. Absolute quantification was performed using:
- Standard curves generated from serial dilutions
- Conversion to log₁₀ absolute copy number
- Correlation analysis between *Xf* and Actinobacterial stearoyl-CoA desaturase abundance

---
