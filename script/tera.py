#!/usr/bin/env python3
# Author: Janky

import argparse
import os
import subprocess
import csv
import gzip
import multiprocessing as mp
import pandas as pd
import numpy as np
import contextlib
import shutil
from datetime import datetime
from collections import defaultdict

@contextlib.contextmanager
def cd(cd_path):
	saved_path = os.getcwd()
	os.chdir(cd_path)
	yield
	os.chdir(saved_path)


def create_STAR_index(args):
	STAR_index = os.path.abspath(args.STAR_index)
	cmd = 'STAR'+' --runMode genomeGenerate' \
		+' --runThreadN '+str(args.nthread) \
		+' --genomeDir '+STAR_index \
		+' --genomeFastaFiles '+os.path.abspath(args.ref_genome) \
		+' --sjdbGTFfile '+os.path.abspath(args.annotation) \
		+' --sjdbOverhang 100' \
		+' --genomeSAindexNbases '+str(args.genomeSAindexNbases)

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def STAR(args):
	STAR_index = os.path.abspath(args.STAR_index)
	fastq1 = os.path.abspath(args.fastq1)
	fastq2 = os.path.abspath(args.fastq2)
	cmd = 'STAR'+' --runThreadN '+str(args.nthread) \
		+' --genomeDir '+STAR_index \
		+' --outFileNamePrefix '+out_dir+'/1_align/'+args.prefix+'_' \
		+' --readFilesIn '+fastq1+' '+fastq2 \
		+' --outFilterType BySJout' \
		+' --outFilterIntronMotifs RemoveNoncanonical' \
		+' --outSAMtype BAM SortedByCoordinate' \
        	+' --outSAMattributes NH HI AS nM NM' \
		+' --twopassMode Basic' \
		+' --outSAMstrandField intronMotif'

	if fastq1.endswith('.gz'):
		cmd += ' --readFilesCommand zcat'
	if args.nthreadsort:
		cmd += ' --outBAMsortingThreadN '+str(args.nthreadsort)
	if args.nRAMsort:
		cmd += ' --limitBAMsortRAM '+str(args.nRAMsort)

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def create_BAM_index(args):
	cmd = 'samtools index'+' -@ '+str(args.nthread)+' '+align_bam

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def filtTE(args):
	STAR_index = os.path.abspath(args.STAR_index)
	ref_TE_gtf = os.path.abspath(args.TE_gtf)
	cmd = 'telescope assign --reassign_mode average --overlap_threshold 0.1 --outdir ./' \
		+' --attribute transcript_id '+align_bam+' '+ref_TE_gtf

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	TE_dict = defaultdict()
	with open('telescope-telescope_report.tsv','r') as fh:
		for line in fh.readlines():
			line = line.strip().split('\t')
			if line[0][0]=='#': continue
			if line[0]=='transcript' or line[0]=='__no_feature': continue

			ID = line[0]
			count = float(line[2])

			TE_dict[ID] = count

	fho = open('Telescope_filt_TE.bed','w')
	with open(ref_TE_gtf,'r') as fh:
		for line in fh.readlines():
			line = line.strip().split('\t')

			if line[0][0]=='#': continue

			chrn = line[0]
			entry_type = line[2]
			start = int(line[3])-1
			end = int(line[4])
			strand = line[6]

			attributes = defaultdict()
			for a in line[8].replace('"', '').split(';')[:-1]:
				kv = a.strip().split(' ')
				attributes[kv[0]] = kv[1]

			TE_id = attributes['transcript_id']
			TE_name = attributes['gene_id']
			if TE_id in TE_dict and TE_dict[TE_id]>0:
				print("%s\t%d\t%d\t%s\t%s\t%s" % (chrn, start, end, TE_id, TE_name, strand), file=fho, flush=True)
	fho.close()

	cmd = 'bedtools slop -i Telescope_filt_TE.bed -g '+STAR_index+'/chrNameLength.txt -b 5000 >Telescope_filt_TE_10K.bed'
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'rm -rf telescope-checkpoint.npz telescope-telescope_report.tsv Telescope_filt_TE.bed'
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def SERVE(args):
	STAR_index = os.path.abspath(args.STAR_index)
	GMAP_index = os.path.abspath(args.GMAP_index)
	ref_fasta = os.path.abspath(args.ref_genome)
	ref_gtf = os.path.abspath(args.annotation)
	ref_TE_bed = os.path.abspath(args.TE_bed)
	ref_TE_gtf = os.path.abspath(args.TE_gtf)
	fastq1 = os.path.abspath(args.fastq1)
	fastq2 = os.path.abspath(args.fastq2)


	cmd = 'python3 '+os.path.join(script_dir, 'SERVE.py') \
		+' -fq1 '+args.fastq1+' -fq2 '+args.fastq2 \
		+' -e Telescope_filt_TE_10K.bed' \
		+' -r '+ref_fasta \
		+' -a '+ref_gtf \
		+' -p '+args.prefix \
		+' -S '+STAR_index \
		+' -G '+GMAP_index \
		+' -g '+args.GMAP_index_name \
		+' -t '+str(args.nthread) \
		+' -m '+args.nRAMassem \
		+' --count 1'
	if args.stranded_type:
		cmd += ' -s '+args.stranded_type

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def StringTie(args):
	cmd = 'stringtie'+' -p '+str(args.nthread) \
		+' -G '+os.path.abspath(args.annotation) \
		+' -o '+STRG_gtf \
		+' '+align_bam

	if args.stranded_type:
		cmd += ' --'+args.stranded_type.lower()

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def TEAM(args):
	cmd = 'python3 '+os.path.join(script_dir, 'team.py') \
		+' -m '+str(args.merge) \
		+' --S1 '+STRG_gtf+' --S2 '+SERVE_gtf \
		+' -r '+os.path.abspath(args.TE_bed) \
		+' -b '+TE_bam \
		+' -o '+out_dir \
		+' -p '+args.prefix

	if args.stranded_type:
		cmd += ' -s'

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def TEA(args):
	cmd = 'python3 '+os.path.join(script_dir, 'tea.py') \
		+' -a '+args.annotation \
		+' -i '+args.input \
		+' --TE_bed '+args.TE_bed \
		+' -o '+out_dir \
		+' -p '+args.prefix \
		+' -d '+str(args.exon_diff)

	if args.TE_exon:
		cmd += ' --TE_exon '+args.TE_exon

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')


def TEEN(args):
	TEEN_index = os.path.abspath(args.index)
	ref_fasta = os.path.abspath(args.ref_genome)
	ref_gtf = os.path.abspath(args.annotation)
	ref_TE_gtf = os.path.abspath(args.TE_gtf)
	fastq1 = os.path.abspath(args.fastq1)
	fastq2 = os.path.abspath(args.fastq2)

	cmd = 'python3 '+os.path.join(script_dir, 'teen.py') \
		+' -r '+ref_fasta \
		+' -a '+ref_gtf \
		+' --TE_gtf '+ref_TE_gtf \
		+' -i '+TEEN_index \
		+' -t '+str(args.nthread) \
		+' -fq1 '+fastq1 \
		+' -fq2 '+fastq2 \
		+' -o '+out_dir \
		+' -p '+args.prefix \

	if args.level:
		cmd += ' --level '+args.level
	if args.quant:
		cmd += ' --quant '+args.quant
	if args.stranded_type:
		cmd += ' -s '+args.stranded_type

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def detect(args):
	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] TERA detection start.', flush=True)

	STAR_index = os.path.abspath(args.STAR_index)
	GMAP_index = os.path.abspath(args.GMAP_index)

	with cd(out_dir):
		if not os.path.exists(STAR_index):
			os.makedirs(STAR_index)
			if not args.annotation:
				print('ERROR: Lack annotation file (--annotation)')
			if not args.ref_genome:
				print('ERROR: Lack reference genome (--ref_genome)')

			print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to create STAR index.', flush=True)

			create_STAR_index(args)

			print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish (In ./STAR_index/).', flush=True)

		if not os.path.exists(GMAP_index):
			os.makedirs(GMAP_index)
			if not args.ref_genome:
				print('ERROR: Lack reference genome (--ref_genome)')
			print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to create GMAP index.', flush=True)

			create_GMAP_index(args)

			print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish (In ./GMAP_index/).', flush=True)


		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to align RNA-seq reads to reference genome.', flush=True)

		if not os.path.exists(align_dir):
			os.makedirs(align_dir)

		STAR(args)

		create_BAM_index(args)

		filtTE(args)
		
		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish.', flush=True)


		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to de novo assemble transcripts by SERVE.', flush=True)

		SERVE(args)

		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish.', flush=True)


		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to assemble transcripts by StringTie.', flush=True)

		StringTie(args)

		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish.', flush=True)


		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to run TEAM.', flush=True)

		TEAM(args)

		cmd = 'rm -rf Telescope_filt_TE_10K.bed '+align_dir+' '+assem_dir+' '+qc_dir
		subprocess.check_call(cmd, shell=True, executable='/bin/bash')

		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish.', flush=True)

	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] TERA identification is done.', flush=True)

def anno(args):
	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] TERA annotation start.', flush=True)

	TEA(args)

	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] TERA annotation is done.', flush=True)

def quant(args):
	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] TERA quantification start.', flush=True)

	TEEN(args)

	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] TERA quantification is done.', flush=True)


parser = argparse.ArgumentParser(description='TERA: Transposable Element-derived RNA Analysis')
subparsers = parser.add_subparsers(help='sub-command help')

parser_detect = subparsers.add_parser('detect', help='detect help')
parser_detect.add_argument('-fq1', '--fastq1', help='Read1 in FASTQ format (required)', required=True)
parser_detect.add_argument('-fq2', '--fastq2', help='Read1 in FASTQ format (required)', required=True)
parser_detect.add_argument('--TE_bed', help='TE position in BED format (required)', required=True)
parser_detect.add_argument('--TE_gtf', help='TE position in GTF format (required)', required=True)
parser_detect.add_argument('-r', '--ref_genome', help='Reference genome in FASTA format (required)', required=True)
parser_detect.add_argument('-a', '--annotation', help='Genome annotation in GTF format (required)', required=True)
parser_detect.add_argument('-s', '--stranded_type', help='Strand-specific RNA-seq read orientation: RF or FR', choices=['RF', 'FR'])
parser_detect.add_argument('-o', '--output_dir', default='.', help='Output directory (default: .)')
parser_detect.add_argument('-p', '--prefix', default='TERA', help='Prefix for output file name (default: TERA)')
parser_detect.add_argument('-m', '--merge', default=1, type=int, help='Merge pattern. 1: local, 2: global (default: 1)', choices=[1, 2])
parser_detect.add_argument('-S', '--STAR_index', default='./STAR_index', help='Path to the directory where STAR index generated (default: STAR_index)')
parser_detect.add_argument('-G', '--GMAP_index', default='./GMAP_index', help='Path to the directory where GMAP index generated (default: GMAP_index)')
parser_detect.add_argument('-g', '--GMAP_index_name', default='GRCh38', help='GMAP index name (default: GRCh38)')
parser_detect.add_argument('-t', '--nthread', type=int, default=1, help='Number of threads to run TERA (default: 1)')
parser_detect.add_argument('--genomeSAindexNbases', default=14, type=int, help='length (bases) of the SA pre-indexing string for creating STAR index. Typically between 10 and 15. For small genomes, this parameter must be scaled down to min(14, log2(GenomeLength)/2-1)')
parser_detect.add_argument('--nthreadsort', type=int, help='Number of threads for BAM sorting')
parser_detect.add_argument('--nRAMsort', type=int, default=10000000000, help='Maximum available RAM (bytes) for sorting BAM (default: 10000000000).')
parser_detect.add_argument('--nRAMassem', default='10G', help='Maximum available RAM (Gb) for assembly (default: 10G)')
parser_detect.add_argument('--max_intron', default=200000, type=int, help='Maximum intron length of transcripts (default: 200000)')
parser_detect.add_argument('--min_identity', default=0.95, help='Minimum identity of assembled transcripts (default: 0.95)')
parser_detect.add_argument('--min_coverage', default=0.95, help='Minimum coverage of assembled transcripts (default: 0.95)')

parser_detect.set_defaults(func=detect)

parser_anno = subparsers.add_parser('anno', help='anno help')
parser_anno.add_argument('-i', '--input', help='Input file of TE transcripts in GTF format (required)')
parser_anno.add_argument('--TE_bed', help='TE position in BED format (required)')
parser_anno.add_argument('-o', '--output_dir', default='.', help='Output directory (default: .)')
parser_anno.add_argument('-p', '--prefix', default='TEA', help='Prefix for output file name (default: TEA)')
parser_anno.add_argument('-a', '--annotation', help='Genome annotation in GTF format (required)')
parser_anno.add_argument('--TE_exon', help='TE exon annotation in BED format')
parser_anno.add_argument('-d', '--exon_diff', type=int, default=5, help='Maximum difference (bp) of exon ends (default: 5)')

parser_anno.set_defaults(func=anno)

parser_quant = subparsers.add_parser('quant', help='quant help')
parser_quant.add_argument('-fq1', '--fastq1', help='Read1 in FASTQ format (required)')
parser_quant.add_argument('-fq2', '--fastq2', help='Read1 in FASTQ format (required)')
parser_quant.add_argument('--TE_gtf', help='TE position in GTF format (required)')
parser_quant.add_argument('-l', '--level', default=1, help='TE level: 1 for exon-level, 2 for unit-level and family-level (default: 1)')
parser_quant.add_argument('-s', '--stranded_type', help='Strand-specific RNA-seq read orientation: RF or FR')
parser_quant.add_argument('-o', '--output_dir', default='.', help='Output directory (default: .)')
parser_quant.add_argument('-p', '--prefix', default='TEEN', help='Prefix for output file name (default: TEEN)')
parser_quant.add_argument('-r', '--ref_genome', help='Reference genome in FASTA format (required)')
parser_quant.add_argument('-a', '--annotation', help='Genome annotation in GTF format (required)')
parser_quant.add_argument('--TE_exon', help='TE exon annotation in BED format')
parser_quant.add_argument('-q', '--quant', help='Quantification pattern: kallisto or rsem (required for exon-level quantification)')
parser_quant.add_argument('-i', '--index', default='./TEEN_index/TEEN', help='Index name (default: ./TEEN_index/TEEN)')
parser_quant.add_argument('-t', '--nthread', type=int, default=1, help='Number of threads to run TEEN (default: 1)')

parser_quant.set_defaults(func=quant)

args = parser.parse_args()
script_dir = os.path.abspath(os.path.dirname(__file__))
out_dir = os.path.abspath(args.output_dir)
align_dir = out_dir+'/1_align'
assem_dir = out_dir+'/2_assem'
qc_dir = out_dir+'/3_qc'


align_bam = align_dir+'/'+args.prefix+'_Aligned.sortedByCoord.out.bam'
TE_bam = align_dir+'/'+args.prefix+'_TE.bam'

SERVE_gtf = qc_dir+'/'+args.prefix+'_SERVE.gtf'
STRG_gtf = qc_dir+'/'+args.prefix+'_STRG.gtf'
TEAM_gtf = out_dir+'/'+args.prefix+'_TEAM.gtf'

if not os.path.exists(out_dir):
    os.makedirs(out_dir)


print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running TEA.', flush=True)

args.func(args)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Congratulations!!! TEA Finished.', flush=True)
