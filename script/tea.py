#!/usr/bin/env python3
# Author: Janky

import argparse
import os
import subprocess
import gzip
import multiprocessing as mp
import pandas as pd
import numpy as np
import contextlib
import shutil
from datetime import datetime

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
		+' --genomeFastaFiles '+ref_fasta \
		+' --sjdbGTFfile '+ref_gtf \
		+' --sjdbOverhang 100' \
		+' --genomeSAindexNbases '+str(args.genomeSAindexNbases)

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def STAR(args):
	STAR_index = os.path.abspath(args.STAR_index)
	cmd = 'STAR'+' --runThreadN '+str(args.nthread) \
		+' --genomeDir '+STAR_index \
		+' --outFileNamePrefix '+out_dir+'/'+args.prefix+'_' \
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

def Trinity(args):
	cmd = 'Trinity'+' --genome_guided_bam '+align_bam \
		+' --genome_guided_max_intron '+str(args.max_intron) \
		+' --CPU '+str(args.nthread) \
		+' --max_memory '+args.nRAMassem \
		+' --output '+'./trinity_all'

	if args.stranded_type:
		cmd += ' --SS_lib_type '+args.stranded_type

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')
	shutil.copyfile('./trinity_all/Trinity-GG.fasta',Trinity_fasta)
	shutil.rmtree('./trinity_all')

def create_GMAP_index(args):
	GMAP_index = os.path.abspath(args.GMAP_index)
	cmd = 'gmap_build'+' -D '+GMAP_index \
		+' -d '+args.GMAP_index_name \
		+' '+ref_fasta

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def Remap(args):
	GMAP_index = os.path.abspath(args.GMAP_index)
	cmd = 'gmap'+' -t '+str(args.nthread) \
		+' -D '+GMAP_index \
		+' -d '+args.GMAP_index_name \
		+' -f gff3_gene' \
		+' --max-intronlength-middle='+str(args.max_intron) \
		+' --no-chimeras' \
		+' --min-identity='+str(args.min_identity) \
		+' --min-trimmed-coverage='+str(args.min_coverage) \
		+' '+Trinity_fasta \
		+' > '+Trinity_gff3

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'gffread'+' -i '+str(args.max_intron) \
		+' --sort-alpha -T' \
		+' -o '+Trinity_gtf \
		+' '+Trinity_gff3

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def StringTie(args):
	cmd = 'stringtie'+' -p '+str(args.nthread) \
		+' -G '+ref_gtf \
		+' -o '+STRG_gtf \
		+' '+align_bam

	if args.stranded_type:
		cmd += ' --'+args.stranded_type.lower()

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def TEAM(args):
	cmd = 'python3 '+os.path.join(script_dir, 'team.py') \
		+' -m '+str(args.merge) \
		+' -S '+STRG_gtf+' -T '+Trinity_gtf \
		+' -r '+ref_TE_bed \
		+' -b '+align_bam \
		+' -o '+out_dir \
		+' -p '+args.prefix

	if args.stranded_type:
		cmd += ' -s'

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def TEEN(args):
	cmd = 'python3 '+os.path.join(script_dir, 'teen.py') \
		+' -r '+ref_fasta \
		+' -a '+ref_gtf \
		+' -e '+ref_TE_bed \
		+' -t '+str(args.nthread) \
		+' -fq1 '+fastq1 \
		+' -fq2 '+fastq2 \
		+' -o '+out_dir \
		+' -p '+args.prefix

	if args.kallisto:
		cmd += ' --quant kallisto'
	if args.rsem:
		cmd += ' --quant rsem'
	if args.stranded_type:
		cmd += ' -s '+args.stranded_type
	if args.TE_exon:
		cmd += ' --TE_exon '+args.TE_exon

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def detect(args):
	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] TEA identification start.', flush=True)

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

		STAR(args)

		cmd = 'rm -rf '+args.prefix+'__STARpass1 '+args.prefix+'__STARgenome '+args.prefix+'_SJ.out.tab '+args.prefix+'_Log.progress.out '+args.prefix+'_Log.out Log.out'
		subprocess.check_call(cmd, shell=True, executable='/bin/bash')

		create_BAM_index(args)

		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish.', flush=True)


		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to de novo assemble transcripts by Trinity.', flush=True)

		Trinity(args)
		Remap(args)

		cmd = 'rm -rf '+Trinity_gff3+' '+Trinity_fasta
		subprocess.check_call(cmd, shell=True, executable='/bin/bash')

		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish.', flush=True)


		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to assemble transcripts by StringTie.', flush=True)

		StringTie(args)

		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish.', flush=True)


		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to run TEAM.', flush=True)

		TEAM(args)

		cmd = 'gffread -T --sort-alpha -o '+out_gtf+' '+TEAM_gtf
		subprocess.check_call(cmd, shell=True, executable='/bin/bash')

		cmd = 'rm -rf '+STRG_gtf+' '+Trinity_gtf+' '+TEAM_gtf
		subprocess.check_call(cmd, shell=True, executable='/bin/bash')

		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish.', flush=True)

	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] TEA identification is done.', flush=True)

def quant(args):
	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] TEA quantification start.', flush=True)

	TEEN(args)

	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] TEA quantification is done.', flush=True)


parser = argparse.ArgumentParser(description='TEA: pipeline for Transposable element analysis')
subparsers = parser.add_subparsers(help='sub-command help')

parser_detect = subparsers.add_parser('detect', help='detect help')
parser_detect.add_argument('-fq1', '--fastq1', help='Read1 in FASTQ format (required)', required=True)
parser_detect.add_argument('-fq2', '--fastq2', help='Read1 in FASTQ format (required)', required=True)
parser_detect.add_argument('-e', '--TE', help='TE position in BED format (required)', required=True)
parser_detect.add_argument('-r', '--ref_genome', help='Reference genome in FASTA format (required)', required=True)
parser_detect.add_argument('-a', '--annotation', help='Genome annotation in GTF format (required)', required=True)
parser_detect.add_argument('-s', '--stranded_type', help='Strand-specific RNA-seq read orientation: RF or FR', choices=['RF', 'FR'])
parser_detect.add_argument('-o', '--output_dir', default='.', help='Output directory (default: .)')
parser_detect.add_argument('-p', '--prefix', default='TEA', help='Prefix for output file name (default: TEA)')
parser_detect.add_argument('-m', '--merge', default=1, type=int, help='Merge pattern. 1: local, 2: global (default: 1)', choices=[1, 2])
parser_detect.add_argument('-S', '--STAR_index', default='./STAR_index', help='Path to the directory where STAR index generated (default: STAR_index)')
parser_detect.add_argument('-G', '--GMAP_index', default='./GMAP_index', help='Path to the directory where GMAP index generated (default: GMAP_index)')
parser_detect.add_argument('-g', '--GMAP_index_name', default='GRCh38', help='GMAP index name (default: GRCh38)')
parser_detect.add_argument('-t', '--nthread', type=int, default=1, help='Number of threads to run TEA (default: 1)')
parser_detect.add_argument('--genomeSAindexNbases', default=14, type=int, help='length (bases) of the SA pre-indexing string for creating STAR index. Typically between 10 and 15. For small genomes, this parameter must be scaled down to min(14, log2(GenomeLength)/2-1)')
parser_detect.add_argument('--nthreadsort', type=int, help='Number of threads for BAM sorting')
parser_detect.add_argument('--nRAMsort', type=int, default=10000000000, help='Maximum available RAM (bytes) for sorting BAM (default: 10000000000).')
parser_detect.add_argument('--nRAMassem', default='10G', help='Maximum available RAM (Gb) for assembly (default: 10G)')
parser_detect.add_argument('--max_intron', default=200000, type=int, help='Maximum intron length of transcripts (default: 200000)')
parser_detect.add_argument('--min_identity', default=0.95, help='Minimum identity of assembled transcripts (default: 0.95)')
parser_detect.add_argument('--min_coverage', default=0.95, help='Minimum coverage of assembled transcripts (default: 0.95)')

parser_detect.set_defaults(func=detect)

parser_quant = subparsers.add_parser('quant', help='quant help')
parser_quant.add_argument('-fq1', '--fastq1', help='Read1 in FASTQ format (required)', required=True)
parser_quant.add_argument('-fq2', '--fastq2', help='Read1 in FASTQ format (required)', required=True)
parser_quant.add_argument('-e', '--TE', help='TE position in BED format (required)', required=True)
parser_quant.add_argument('-r', '--ref_genome', help='Reference genome in FASTA format (required)', required=True)
parser_quant.add_argument('-a', '--annotation', help='Genome annotation in GTF format (required)', required=True)
parser_quant.add_argument('-s', '--stranded_type', help='Strand-specific RNA-seq read orientation: RF or FR', choices=['RF', 'FR'])
parser_quant.add_argument('-o', '--output_dir', default='.', help='Output directory (default: .)')
parser_quant.add_argument('-p', '--prefix', default='TEA', help='Prefix for output file name (default: TEA)')
parser_quant.add_argument('-t', '--nthread', type=int, default=1, help='Number of threads to run TEA (default: 1)')
parser_quant.add_argument('--kallisto', help='Specific Quantification by kallisto', action='store_true')
parser_quant.add_argument('--rsem', help='Specific Quantification by RSEM', action='store_true')
parser_quant.add_argument('--TE_exon', help='TE exon annotation in BED format (Only for quantification analysis)')

parser_quant.set_defaults(func=quant)

args = parser.parse_args()
script_dir = os.path.abspath(os.path.dirname(__file__))
out_dir = os.path.abspath(args.output_dir)

ref_fasta = os.path.abspath(args.ref_genome)
ref_gtf = os.path.abspath(args.annotation)
ref_TE_bed = os.path.abspath(args.TE)

fastq1 = os.path.abspath(args.fastq1)
fastq2 = os.path.abspath(args.fastq2)

align_bam = out_dir+'/'+args.prefix+'_Aligned.sortedByCoord.out.bam'

Trinity_fasta = out_dir+'/'+args.prefix+'_Trinity.fasta'
Trinity_gff3 = out_dir+'/'+args.prefix+'_Trinity.gff3'
Trinity_gtf = out_dir+'/'+args.prefix+'_Trinity.gtf'
STRG_gtf = out_dir+'/'+args.prefix+'_STRG.gtf'
TEAM_gtf = out_dir+'/'+args.prefix+'_TEAM.gtf'
out_gtf = out_dir+'/'+args.prefix+'.gtf'

if not os.path.exists(out_dir):
    os.makedirs(out_dir)


print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running TEA on {0:d} threads.'.format(args.nthread), flush=True)

args.func(args)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Congratulations!!! TEA Finished.', flush=True)
