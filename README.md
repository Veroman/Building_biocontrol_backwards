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

```
config.sh                          paths, database locations, and every analysis threshold
environment.yml                    conda environment with pinned software versions
delta9.faa                         File with Delta9 proteins for the database.

scripts/
  00_download.sh                   SRA download (prefetch + fasterq-dump)
  01_qc_fastp.sh                   adapter / quality trimming
  02_kraken2_xylella.sh            Kraken2
  03_build_desaturase_db.sh        custom UniProt desaturase database + Kaiju index
  04_kaiju_desaturase.sh           Kaiju search of host-free reads; lineage of each hit
  05_summarize.py                  builds every processed table in results/
```
## Software and databases

| Step | Software (version) | Database (version / date) |
|---|---|---|
| Download | SRA Toolkit (3.1.1) | NCBI SRA |
| Read QC | fastp (0.23.4) | – |
| Taxonomic classification | Kraken2 (2.1.2), Bracken (2.9), KrakenTools (1.2) | Kraken2 PlusPF, build `k2_pluspf_2022` |
| Delta9 search | Kaiju (1.10.1) | Custom UniProtKB database, release `20XX_XX` (downloaded YYYY-MM-DD); NCBI taxonomy downloaded YYYY-MM-DD |
| Statistics | Python 3.11, SciPy (1.13); R (4.x) | – |

*Versions above are placeholders taken from `environment.yml`; replace them with the versions in `work/logs/versions_*.txt`.*

## Running the pipeline

```bash
conda env create -f environment.yml
conda activate biocontrol
# edit paths in config.sh, then:
source config.sh
bash scripts/00_download.sh
bash scripts/01_qc_fastp.sh
bash scripts/03_kraken2_xylella.sh
bash scripts/04_build_desaturase_db.sh
bash scripts/05_kaiju_desaturase.sh
python3 scripts/06_summarize.py
Rscript scripts/07_figures.R
`
## Analysis workflow

### 1. Data

The 59 paired-end shotgun metagenomes listed in the repository were downloaded with `prefetch` and `fasterq-dump`. 

### 2. Read quality control

Reads were trimmed with fastp: automatic adapter detection for paired-end data, sliding-window trimming from the 3′ end (window 4, mean Q ≥ 20), removal of reads < 50 bp, with > 5 ambiguous bases, or of low complexity.

### 3. Taxonomic classification and Xylella detection

Reads were classified with Kraken2 against the PlusPF database (archaea, bacteria, viruses, plasmids, human, UniVec_Core, protozoa and fungi) with `--confidence 0.1`. Species abundances were re-estimated with Bracken (read length 150, species level, threshold 10 reads). Bacterial and fungal fractions were taken from the Kraken2 report as reads assigned to the domain Bacteria (taxid 2) and kingdom Fungi (taxid 4751).

### 4. Custom stearoyl-CoA desaturase database

UniProtKB was queried for proteins named "stearoyl-CoA desaturase" or with EC 1.14.19.1, restricted to Bacteria (taxid 2) and Fungi (taxid 4751):

```
((protein_name:"stearoyl-CoA desaturase") OR (ec:1.14.19.1)) AND ((taxonomy_id:2) OR (taxonomy_id:4751))
```

Plant and animal sequences were excluded so that residual host reads cannot match the database. Sequences shorter than 250 aa or longer than 500 aa were removed as likely fragments or fusion proteins. FASTA headers were rewritten as `<UniProt accession>_<NCBI taxid>` as required by Kaiju, and the index was built with `kaiju-mkbwt` and `kaiju-mkfmi`. The database is in `data/delta9_db.faa`.

### 6. Desaturase gene search

Host-free read pairs were searched against the database with Kaiju in greedy mode with: at most 3 mismatches (`-e 3`), minimum match score 75 (`-s 75`), minimum match length 20 aa (`-m 20`), maximum E-value 0.001 (`-E 0.001`), and SEG low-complexity filtering (default). Reads matching several database sequences equally well were assigned to the lowest common ancestor (LCA) of those sequences (Kaiju default). The superkingdom, phylum, class and genus of each hit were added with `kaiju-addTaxonNames`. A read was counted as an Actinomycetota Δ9-desaturase read when its LCA was within the phylum Actinomycetota (taxid 201174); reads assigned above phylum level were counted in the total but not attributed to a phylum.




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
