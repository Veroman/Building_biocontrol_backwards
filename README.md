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

# Step 2 - Clean raw reads
Clean reads with low Phred scores.

```bash
TRIMMOMATIC=/path/to/Trimmomatic-0.41
THREADS=8
OUT2=/path/to/qc_reads 
mkdir -p ${OUT2}/unpaired

for R1 in ${OUT}/*_1.fastq.gz; do
    S=$(basename $R1 _1.fastq.gz)
    java -jar ${TRIMMOMATIC}/trimmomatic-0.39.jar PE \
        -threads ${THREADS} -phred33 \
        -summary ${OUT}/${S}.trimmomatic.summary \
        ${OUT}/${S}_1.fastq.gz ${OUT}/${S}_2.fastq.gz \
        ${OUT2}/${S}_1.fastq.gz ${OUT2}/unpaired/${S}_1U.fastq.gz \
        ${OUT2}/${S}_2.fastq.gz ${OUT2}/unpaired/${S}_2U.fastq.gz \
        ILLUMINACLIP:${TRIMMOMATIC}/adapters/TruSeq3-PE-2.fa:2:30:10:2:True \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 \
        2> ${OUT2}/${S}.trimmomatic.log
done
```

# Step 3 - Classified reads using Kraken2 and PlusPF
The reads were classified with Kraken2 against the PlusPF database (archaea, bacteria, viruses, plasmids, human, UniVec_Core, protozoa and fungi) downloaded in 2022.

```bash
#!/usr/bin/env bash

KRAKEN2=/path/to/kraken2/
KRAKEN_DB= /path/PlusPF
OUT3=/path/to/kraken2
cd ${KRAKEN2}

for R1 in "${OUT}"/*_1.fastq.gz; do
    S=$(basename "$R1" _1.fastq.gz)

    ./kraken2 --db "$KRAKEN_DB" --threads "$THREADS" \
        --paired --gzip-compressed \
        --report "${OUT3}/${S}.kreport" \
        --output "${OUT3}/${S}.kraken" \
       ${IN}/${S}_1.fastq.gz ${IN}/${S}_2.fastq.gz
done
```

# Step 4 - Search for reads associated with stearoyl-CoA desaturase using Kaiju and Uniprot
UniProtKB was queried for proteins named "stearoyl-CoA desaturase" or with EC 1.14.19.1, restricted to Bacteria (taxid 2) and Fungi (taxid 4751):
((protein_name:"stearoyl-CoA desaturase") OR (ec:1.14.19.1)) AND ((taxonomy_id:2) OR (taxonomy_id:4751))

The sequences shorter than 250 aa or longer than 500 aa were removed as likely fragments or fusion proteins. FASTA headers were rewritten as `<UniProt accession>_<NCBI taxid>` as required by Kaiju, and the index was built with `kaiju-mkbwt` and `kaiju-mkfmi`. The database is in `data/delta9_db.faa`.

```bash
KAIJU=/path/to/kaiju/bin
DB_FASTA=/data/delta9_db.faa
TAXONOMY=/path/to/taxonomy     # folder with nodes.dmp and names.dmp
DB_DIR=/path/to/kaiju_db
OUT4=/path/to/kaiju_output
THREADS=8

cd ${KAIJU}

# 1. Create index
./kaiju-mkbwt -n ${THREADS} \
    -a ACDEFGHIKLMNPQRSTVWY \
    -o ${DB_DIR}/proteins \
    ${DB_FASTA}
./kaiju-mkfmi ${DB_DIR}/proteins

# 2. Search for delta9 reads
for R1 in ${OUT2}/*_1.fastq.gz; do
    S=$(basename $R1 _1.fastq.gz)

    ./kaiju -z ${THREADS} \
        -t ${TAXONOMY}/nodes.dmp \
        -f ${DB_DIR}/proteins.fmi \
        -i ${IN}/${S}_1.fastq.gz \
        -j ${IN}/${S}_2.fastq.gz \
        -o ${OUT4}/${S}.kaiju.out \
        -a greedy \
        -e 3 \
        -s 75 \
        -m 20 \
        -E 0.001 \
        -v
# 3. Classified reads + lineage

    grep '^C' ${OUT4}/${S}.kaiju.out > ${OUT4}/${S}.kaiju.classified.out
    ./kaiju-addTaxonNames \
        -t ${TAXONOMY}/nodes.dmp -n ${TAXONOMY}/names.dmp \
        -r superkingdom,phylum,class,genus \
        -i ${OUT4}/${S}.kaiju.classified.out \
        -o ${OUT4}/${S}.kaiju.names.out

# 4. Read counts per phylum
    ./kaiju2table -t ${TAXONOMY}/nodes.dmp -n ${TAXONOMY}/names.dmp \
        -r phylum -o ${OUT4}/${S}.phylum.tsv ${OUT4}/${S}.kaiju.out
done
```

**Parameters:**
- `-f`: FM-index path
- `-i, -j`: Forward and reverse sequencing files
- `-o`: Output classification file
- `-v`: Verbose output
- `-s`: Sensitivity (60 = high sensitivity)
- `-e`: Maximum number of E-values to report
- `-m`: Minimum match length (amino acids)
- `-E`: Maximum E-value
