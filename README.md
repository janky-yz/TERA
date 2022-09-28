# TEA
A developed pipeline for Transposable element analysis.

This repository contains TEA script (script/) and TE reference files (hg19, hg38, chm13v2.0)
## Software requirements

Required: python3 (pandas, numpy), samtools, gffread v0.11.6, BEDTools
TEA detect (TE identification): Trinity, StringTie, GMAP
TEA quant (TE quantification): RSEM or kallisto

## Obtaining TEA

```
git clone https://github.com/janky-yz/TEA.git
chmod u+x TEA/script/TEA
```

## Usage

TEA is composed of TE identification (TEA detect) and quantification (TEA quant).

### Input files

```
RNA sequencing files: paired end, in FASTQ format (fastq or fastq.gz)
TE reference bed: TE annotation in BED format (8 fields: chrom, chromStart, chromEnd, ID, name, strand, family, class). See Reference/
Reference genome: reference genome sequnce (FASTA) and annotation (GTF) is required. You can create genome index manually before running TEA.
```

### TEA

```bash

usage: TEA [-h] {detect,quant} ...

TEA: pipeline for Transposable element analysis

positional arguments:
  {detect,quant}  sub-command help
    detect        detect help
    quant         quant help

optional arguments:
  -h, --help      show this help message and exit
```

### TEA detect

TEA detect is designed for TE transcript and exon detection

```bash

Example:
TEA detect -fq1 ${fastq1} -fq2 ${fastq2} -e ${ref_TE_bed} -p ${prefix} -r ${ref_genome_fasta} -a ${ref_genome_gtf} -t ${nthread}
```

```bash

usage: TEA detect [-h] -fq1 FASTQ1 -fq2 FASTQ2 -e TE -r REF_GENOME -a
                      ANNOTATION [-s {RF,FR}] [-o OUTPUT_DIR] [-p PREFIX]
                      [-m MERGE] [-S STAR_INDEX] [-G GMAP_INDEX]
                      [-g GMAP_INDEX_NAME] [-t NTHREAD]
                      [--genomeSAindexNbases GENOMESAINDEXNBASES]
                      [--nthreadsort NTHREADSORT] [--nRAMsort NRAMSORT]
                      [--nRAMassem NRAMASSEM] [--max_intron MAX_INTRON]
                      [--min_identity MIN_IDENTITY]
                      [--min_coverage MIN_COVERAGE]

optional arguments:
  -h, --help            show this help message and exit
  -fq1 FASTQ1, --fastq1 FASTQ1
                        Read1 in FASTQ format (required)
  -fq2 FASTQ2, --fastq2 FASTQ2
                        Read1 in FASTQ format (required)
  -e TE, --TE TE        TE position in BED format (required)
  -r REF_GENOME, --ref_genome REF_GENOME
                        Reference genome in FASTA format (required)
  -a ANNOTATION, --annotation ANNOTATION
                        Genome annotation in GTF format (required)
  -s {RF,FR}, --stranded_type {RF,FR}
                        Strand-specific RNA-seq read orientation: RF or FR
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory (default: .)
  -p PREFIX, --prefix PREFIX
                        Prefix for output file name (default: TEA)
  -m MERGE, --merge MERGE
                        Merge pattern. 1: local, 2: global (default: 1)
  -S STAR_INDEX, --STAR_index STAR_INDEX
                        Path to the directory where STAR index generated
                        (default: STAR_index)
  -G GMAP_INDEX, --GMAP_index GMAP_INDEX
                        Path to the directory where GMAP index generated
                        (default: GMAP_index)
  -g GMAP_INDEX_NAME, --GMAP_index_name GMAP_INDEX_NAME
                        GMAP index name (default: GRCh38)
  -t NTHREAD, --nthread NTHREAD
                        Number of threads to run TEA (default: 1)
  --genomeSAindexNbases GENOMESAINDEXNBASES
                        length (bases) of the SA pre-indexing string for
                        creating STAR index. Typically between 10 and 15. For
                        small genomes, this parameter must be scaled down to
                        min(14, log2(GenomeLength)/2-1)
  --nthreadsort NTHREADSORT
                        Number of threads for BAM sorting
  --nRAMsort NRAMSORT   Maximum available RAM (bytes) for sorting BAM
                        (default: 10000000000).
  --nRAMassem NRAMASSEM
                        Maximum available RAM (Gb) for assembly (default: 10G)
  --max_intron MAX_INTRON
                        Maximum intron length of transcripts (default: 200000)
  --min_identity MIN_IDENTITY
                        Minimum identity of assembled transcripts (default: 0.95)
  --min_coverage MIN_COVERAGE
                        Minimum coverage of assembled transcripts (default: 0.95)
```

In the output directory, you will see output files:
1. ${prefix}.gtf (all transcripts including TE and nonTE)
2. ${prefix}_TE_exon.bed (TE exons): 11 columns for chromosome, start, end, transcript_id, gene_id, strand, exon_type, TE_ID, TE_name, TE_family, TE_class

### TEA quant

TEA quant is designed for TE transcript and exon quantification

```bash

Example:
TEA quant --kallisto/--rsem -fq1 ${fastq1} -fq2 ${fastq2} -e ${ref_TE_bed} -p ${prefix} -r ${ref_genome_fasta} -a ${ref_TEA_gtf} -t ${nthread}
```

```bash

usage: TEA quant [-h] -fq1 FASTQ1 -fq2 FASTQ2 -e TE -r REF_GENOME -a
                     ANNOTATION [-s {RF,FR}] [-o OUTPUT_DIR] [-p PREFIX]
                     [-t NTHREAD] [--kallisto] [--rsem] [--TE_exon TE_EXON]

optional arguments:
  -h, --help            show this help message and exit
  -fq1 FASTQ1, --fastq1 FASTQ1
                        Read1 in FASTQ format (required)
  -fq2 FASTQ2, --fastq2 FASTQ2
                        Read1 in FASTQ format (required)
  -e TE, --TE TE        TE position in BED format (required)
  -r REF_GENOME, --ref_genome REF_GENOME
                        Reference genome in FASTA format (required)
  -a ANNOTATION, --annotation ANNOTATION
                        Genome annotation in GTF format (required)
  -s {RF,FR}, --stranded_type {RF,FR}
                        Strand-specific RNA-seq read orientation: RF or FR
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory (default: .)
  -p PREFIX, --prefix PREFIX
                        Prefix for output file name (default: TEA)
  -t NTHREAD, --nthread NTHREAD
                        Number of threads to run TEA (default: 1)
  --kallisto            Specific Quantification by kallisto
  --rsem                Specific Quantification by RSEM
  --TE_exon TE_EXON     TE exon annotation in BED format (Only for
                        quantification analysis)
```

In the output directory, you will see output files:
1. ${prefix}_transcript_quant.out (quantification results for all transcripts including TE and nonTE): 5 columns for transcript_id, length, eff_length, count, TPM
2. ${prefix}_TE_exon_quant.out (quantification results for TE exon clusters): 3 columns for exon_cluster, count, TPM
3. ${prefix}_TE_exon_anno.bed (information of TE exons): 12 columns for chromosome, start, end, transcript_id, gene_id, strand, exon_type, exon_cluster, TE_ID, TE_name, TE_family, TE_class

### Example pipeline

```bash

Example:
1. TE identification (single sample)
TEA detect -fq1 ${fastq1} -fq2 ${fastq2} -e ${ref_TE_bed} -p test -r ${ref_genome_fasta} -a ${ref_genome_gtf} -t ${nthread}

2. Optional: if you have two or more samples, you can merge all TEA GTF files (generated by step 1) by meta-assembly tools (Cuffmerge is recommended)
ls *gtf >gtf.list
cuffmerge -p ${nthread} gtf.list
gffread -T --sort-alpha -o ref_TEA.gtf merged_asm/merged.gtf

3. TE quantification (single sample)
TEA quant --kallisto -fq1 ${fastq1} -fq2 ${fastq2} -e ${ref_TE_bed} -p test -r ${ref_genome_fasta} -a ${ref_TEA_gtf} -t ${nthread}
```

## Copyright and License Information

Copyright (C) 2021 Jianqi She (janky666@bjmu.edu.cn). See the LICENSE file for license rights and limitations.
