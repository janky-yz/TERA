# TERA
A developed pipeline for Transposable element-derived RNA analysis.

This repository contains TERA script (script/) and TERA reference files (hg19, hg38, chm13v2.0)
## Software requirements

1. Required: Python3 (pandas, numpy), R
2. TERA detect (TE identification): SERVE (v1.0), StringTie (v2.1.1), Telescope (v1.0.3)
3. TERA anno (TE annotation): BEDTools (v2.27.1)
4. TERA quant (TE quantification): Telescope (v1.0.3), RSEM (v1.2.28) or Kallisto (0.44.0)

The software required can be installed via conda:
```
conda env create -n tera -f environment.yml
conda activate tera
conda install -c bioconda gffread=0.12.1 stringtie bedtools
```

Or via singularity
```
singularity build TERA.sif TERA.def
```

## Obtaining TERA

```
git clone https://github.com/janky-yz/TERA.git
chmod u+x TERA/script/TERA
```

## Usage

TERA is composed of TE identification (TERA detect), annotation (TERA anno) and quantification (TERA quant).

### Input files

```
RNA sequencing files: paired end, in FASTQ format (fastq or fastq.gz)
TE reference: TE annotation in BED format (8 fields: chrom, chromStart, chromEnd, ID, name, strand, family, class) and GTF format (use script/TEbedtogtf.R).
Reference genome: reference genome sequnce (FASTA) and annotation (GTF) is required. You can create genome index manually before running TERA.
```

### TERA

```bash

usage: TERA [-h] {detect,anno,quant} ...

TERA: pipeline for Transposable Element-derived RNA Analysis

positional arguments:
  {detect,anno,quant}  sub-command help
    detect        detect help
    anno          anno help
    quant         quant help

optional arguments:
  -h, --help      show this help message and exit
```

### TERA detect

TERA detect is designed for teRNA detection

```bash

Example:
TERA detect -fq1 ${fastq1} -fq2 ${fastq2} --TE_bed ${ref_TE_bed} --TE_gtf ${ref_TE_gtf} -p ${prefix} -r ${ref_genome_fasta} -a ${ref_genome_gtf} -t ${nthread}
```

```bash

usage: tera.py detect [-h] -fq1 FASTQ1 -fq2 FASTQ2 --TE_bed TE_BED --TE_gtf
                      TE_GTF -r REF_GENOME -a ANNOTATION [-s {RF,FR}]
                      [-o OUTPUT_DIR] [-p PREFIX] [-m {1,2}] [-S STAR_INDEX]
                      [-G GMAP_INDEX] [-g GMAP_INDEX_NAME] [-t NTHREAD]
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
  --TE_bed TE_BED       TE position in BED format (required)
  --TE_gtf TE_GTF       TE position in GTF format (required)
  -r REF_GENOME, --ref_genome REF_GENOME
                        Reference genome in FASTA format (required)
  -a ANNOTATION, --annotation ANNOTATION
                        Genome annotation in GTF format (required)
  -s {RF,FR}, --stranded_type {RF,FR}
                        Strand-specific RNA-seq read orientation: RF or FR
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory (default: .)
  -p PREFIX, --prefix PREFIX
                        Prefix for output file name (default: TERA)
  -m {1,2}, --merge {1,2}
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
                        Number of threads to run TERA (default: 1)
  --genomeSAindexNbases GENOMESAINDEXNBASES
                        length (bases) of the SA pre-indexing string for
                        creating STAR index. Typically between 10 and 15. For
                        small genomes, this parameter must be scaled down to
                        min(14, log2(GenomeLength)/2-1)
  --nthreadsort NTHREADSORT
                        Number of threads for BAM sorting
  --nRAMsort NRAMSORT   Maximum available RAM (bytes) for sorting BAM (default: 10000000000).
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
1. ${prefix}.gtf (TE transcripts)
2. ${prefix}_TE_exon.bed (TE exons): 11 columns for chromosome, start, end, transcript_id, gene_id, strand, exon_type, TE_ID, TE_name, TE_family, TE_class

### TERA anno

TERA anno is designed for teRNA annotation

```bash

Example:
TERA anno -i ${TERA_detect_gtf} --TE_exon ${TERA_detect_bed} --TE_bed ${ref_TE_bed} -p ${prefix} -a ${ref_genome_gtf}
```

```bash
usage: tera.py anno [-h] [-i INPUT] [--TE_bed TE_BED] [-o OUTPUT_DIR]
                    [-p PREFIX] [-a ANNOTATION] [--TE_exon TE_EXON]
                    [-d EXON_DIFF]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input file of TE transcripts in GTF format (required)
  --TE_bed TE_BED       TE position in BED format (required)
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory (default: .)
  -p PREFIX, --prefix PREFIX
                        Prefix for output file name (default: TEA)
  -a ANNOTATION, --annotation ANNOTATION
                        Genome annotation in GTF format (required)
  --TE_exon TE_EXON     TE exon annotation in BED format
  -d EXON_DIFF, --exon_diff EXON_DIFF
                        Maximum difference (bp) of exon ends (default: 5)
```

In the output directory, you will see output files:
1. ${prefix}_TE_exon_anno.bed (TE exons): 8 columns for chromosome, start, end, transcript id, gene id, strand, exon type, exon class
2. ${prefix}.TE.exon.anno.txt (annotation for TE exons): 10 columns for chromosome, start, end, transcript id, gene id, strand, exon type, exon class, exon ID, type
3. ${prefix}.TE.unit.anno.txt (annotation for TE units): 10 columns for chromosome, start, end, Dfam id, family id, strand, superfamily id, class id, TE ID, type

### TERA quant

TERA quant is designed for teRNA quantification

```bash

Example:
1. Exon-level: TERA quant -fq1 ${fastq1} -fq2 ${fastq2} -l 1 -q rsem --TE_gtf ${TERA_detect_gtf} --TE_exon ${TERA_anno_bed} -p ${prefix} -r ${ref_genome_fasta} -a ${ref_genome_gtf} -t ${nthread}
2. Unit-level and Family-level: TERA quant -fq1 ${fastq1} -fq2 ${fastq2} -l 2 --TE_gtf ${ref_TE_gtf} -p ${prefix} -r ${ref_genome_fasta} -a ${ref_genome_gtf} -t ${nthread}
```

```bash

usage: tera.py quant [-h] [-fq1 FASTQ1] [-fq2 FASTQ2] [--TE_gtf TE_GTF]
                     [-l LEVEL] [-s STRANDED_TYPE] [-o OUTPUT_DIR] [-p PREFIX]
                     [-r REF_GENOME] [-a ANNOTATION] [--TE_exon TE_EXON]
                     [-q QUANT] [-i INDEX] [-t NTHREAD]

optional arguments:
  -h, --help            show this help message and exit
  -fq1 FASTQ1, --fastq1 FASTQ1
                        Read1 in FASTQ format (required)
  -fq2 FASTQ2, --fastq2 FASTQ2
                        Read1 in FASTQ format (required)
  --TE_gtf TE_GTF       TE position in GTF format (required)
  -l LEVEL, --level LEVEL
                        TE level: 1 for exon-level, 2 for unit-level and
                        family-level (default: 1)
  -s STRANDED_TYPE, --stranded_type STRANDED_TYPE
                        Strand-specific RNA-seq read orientation: RF or FR
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory (default: .)
  -p PREFIX, --prefix PREFIX
                        Prefix for output file name (default: TEEN)
  -r REF_GENOME, --ref_genome REF_GENOME
                        Reference genome in FASTA format (required)
  -a ANNOTATION, --annotation ANNOTATION
                        Genome annotation in GTF format (required)
  --TE_exon TE_EXON     TE exon annotation in BED format
  -q QUANT, --quant QUANT
                        Quantification pattern: kallisto or rsem (required for
                        exon-level quantification)
  -i INDEX, --index INDEX
                        Index name (default: ./TEEN_index/TEEN)
  -t NTHREAD, --nthread NTHREAD
                        Number of threads to run TEEN (default: 1)
```

In the output directory, you will see output files:
1. ${prefix}.transcript.quant.out (quantification results for all transcripts including TE and nonTE): 5 columns for transcript id, length, eff length, count, TPM
2. ${prefix}.TE.exon.quant.out (quantification results for TE exon clusters): 3 columns for exon cluster, count, TPM
3. ${prefix}.TE.unit.quant.out (quantification results for TE units): 2 columns for TE ID, count
4. ${prefix}.TE.family.quant.out (quantification results for TE families): 2 columns for family ID, count

## Copyright and License Information

Copyright (C) 2022 Jianqi She (janky666@bjmu.edu.cn). See the [LICENSE](https://github.com/janky-yz/TERA/blob/main/LICENSE) file for license rights and limitations.
