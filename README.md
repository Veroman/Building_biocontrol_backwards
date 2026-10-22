# Building biocontrol backwards

Scripts, processed data, and supplementary tables for the manuscript:

> **Building biocontrol backwards for forward-thinking plant health management based on efficient mechanisms**
> Cindy E. Morris, Sara Campigli, Guido Marchi, and Veronica Roman-Reyna.

Archived release: <scholarsphere>

## Overview

This study used functional metagenomics to identify microbial genes associated with plant health, using *Xylella fastidiosa* (Xf) as a case study. We searched olive metagenomes for stearoyl-CoA (Δ9) desaturase genes, which we hypothesised would be negatively associated with Xf. 

## Software and databases

| Step | Software (version) | Database (version / date) | Link
|---|---|---|---|
| Download | SRA Toolkit (3.1.1) | NCBI SRA |https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump
| Read QC | Trimmomatic (0.41) | – |https://github.com/usadellab/trimmomatic
| Taxonomic classification | Kraken2 (2.1.2) | Kraken2 PlusPF, build `k2_pluspf_2022` | https://github.com/DerrickWood/kraken2/wiki
| Delta9 search | Kaiju (1.10.1) | Custom UniProtKB database, (accessed 2022) |https://github.com/bioinformatics-centre/kaiju

## Running the pipeline

# Step 1 - Download raw reads from NCBI SRA

```bash
#!/usr/bin/env bash
# 1.1) prefetch  1.2) fasterq-dump

SRATOOLKIT=/PROGRAMS/sratoolkit.3.1.1-ubuntu64/bin
OUT=/path/to/raw_data      # change to your output folder
LIST=path/to/accessions.txt
mkdir -p ${OUT}

cd ${SRATOOLKIT}

# 1.1 Download .sra files
./prefetch --option-file ${LIST} --output-directory ${OUT} --max-size u

# 1.2 Convert to paired FASTQ and compress
for i in $(cat ${LIST}); do
    ./fasterq-dump ${OUT}/${i} --split-files --threads 8 --outdir ${OUT}
    gzip ${OUT}/${i}_1.fastq ${OUT}/${i}_2.fastq
done
```

# Step 2 - Kraken2 

```bash
#!/usr/bin/env bash

KRAKEN2=/path/to/kraken2/
KRAKEN_DB= /path/PlusPF

for R1 in "${HOSTFREE_DIR}"/*_1.fastq.gz; do
    S=$(basename "$R1" _1.fastq.gz)
    R2="${HOSTFREE_DIR}/${S}_2.fastq.gz"

    # 3a Kraken2
    kraken2 --db "$KRAKEN_DB" --threads "$THREADS" \
        --paired --gzip-compressed \
        --report "${KRAKEN_DIR}/${S}.kreport" \
        --output "${KRAKEN_DIR}/${S}.kraken" \
        "$R1" "$R2"
```
# Step 3 - Kraken2 


### 2. Read quality control

Reads were trimmed with Trimmomatic.

### 3. Taxonomic classification and Xylella detection

The reads were classified with Kraken2 against the PlusPF database (archaea, bacteria, viruses, plasmids, human, UniVec_Core, protozoa and fungi) downloaded in 2024. Bacterial and fungal fractions were taken from the Kraken2 report as reads assigned to the domain Bacteria (taxid 2) and kingdom Fungi (taxid 4751).

### 4. Custom stearoyl-CoA desaturase database

UniProtKB was queried for proteins named "stearoyl-CoA desaturase" or with EC 1.14.19.1, restricted to Bacteria (taxid 2) and Fungi (taxid 4751):

```
((protein_name:"stearoyl-CoA desaturase") OR (ec:1.14.19.1)) AND ((taxonomy_id:2) OR (taxonomy_id:4751))
```

Plant and animal sequences were excluded so that residual host reads cannot match the database. Sequences shorter than 250 aa or longer than 500 aa were removed as likely fragments or fusion proteins. FASTA headers were rewritten as `<UniProt accession>_<NCBI taxid>` as required by Kaiju, and the index was built with `kaiju-mkbwt` and `kaiju-mkfmi`. The database is in `data/delta9_db.faa`.

### 6. Desaturase gene search

Host-free read pairs were searched against the database with Kaiju in greedy mode with: at most 3 mismatches (`-e 3`), minimum match score 75 (`-s 75`), minimum match length 20 aa (`-m 20`), maximum E-value 0.001 (`-E 0.001`), and SEG low-complexity filtering (default). Reads matching several database sequences equally well were assigned to the lowest common ancestor (LCA) of those sequences (Kaiju default). The superkingdom, phylum, class and genus of each hit were added with `kaiju-addTaxonNames`. A read was counted as an Actinomycetota Δ9-desaturase read when its LCA was within the phylum Actinomycetota (taxid 201174); reads assigned above phylum level were counted in the total but not attributed to a phylum.

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
