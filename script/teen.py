#!/usr/bin/env python3
# Author: Janky

import argparse
import os
import subprocess
import gzip
import multiprocessing as mp
import pandas as pd
import numpy as np
from collections import defaultdict
import contextlib
import shutil
from datetime import datetime


@contextlib.contextmanager
def cd(cd_path):
	saved_path = os.getcwd()
	os.chdir(cd_path)
	yield
	os.chdir(saved_path)


def create_kallisto_index(args):
	cmd = 'cat '+ref_gtf+' '+ref_TE_gtf+' >merged.gtf'
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'gffread'+' -g '+ref_fasta \
		+' -w '+transcript_fasta \
		+' merged.gtf'

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'kallisto index'+' -i '+Index \
		+' '+transcript_fasta

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def kallisto(args):
	cmd = 'kallisto quant'+' -i '+Index \
		+' -o '+out_dir \
		+' -t '+str(args.nthread)

	if args.stranded_type:
		cmd += ' --'+args.stranded_type.lower()+'-stranded'

	cmd += ' ' + fastq1 + ' ' + fastq2

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def create_rsem_index(args):
	cmd = 'cat '+ref_gtf+' '+ref_TE_gtf+' >merged.gtf'
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'rsem-prepare-reference'+' -p '+str(args.nthread) \
		+' --star --gtf merged.gtf' \
		+' '+ref_fasta \
		+' '+Index

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def rsem(args):
	cmd = 'rsem-calculate-expression'+' -p '+str(args.nthread) \
		+' --paired-end --star '

	if args.stranded_type=='RF':
		cmd += ' --forward-prob 0'
	elif args.stranded_type=='FR':
		cmd += ' --forward-prob 1'

	if fastq1.endswith('.gz'):
		cmd += ' --star-gzipped-read-file '+fastq1+' '+fastq2
	else:
		cmd += ' '+fastq1+' '+fastq2

	cmd += ' '+Index+' '+args.prefix

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def create_STAR_index(args):
	cmd = 'STAR'+' --runMode genomeGenerate' \
		+' --runThreadN '+str(args.nthread) \
		+' --genomeDir '+ref_dir \
		+' --genomeFastaFiles '+ref_fasta \
		+' --sjdbGTFfile '+ref_gtf \
		+' --sjdbOverhang 100'

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def STAR(args):
	cmd = 'STAR'+' --runThreadN '+str(args.nthread) \
		+' --genomeDir '+ref_dir \
		+' --outFileNamePrefix '+args.prefix+'_' \
		+' --readFilesIn '+fastq1+' '+fastq2 \
		+' --outFilterType BySJout' \
		+' --outFilterIntronMotifs RemoveNoncanonical' \
		+' --outSAMtype BAM SortedByCoordinate' \
        	+' --outSAMattributes NH HI AS nM NM' \
		+' --twopassMode Basic' \
		+' --outSAMstrandField intronMotif'

	if fastq1.endswith('.gz'):
		cmd += ' --readFilesCommand zcat'

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'samtools index'+' -@ '+str(args.nthread)+' '+args.prefix+'_Aligned.sortedByCoord.out.bam'

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')


def telescope(args):
	align_bam = args.prefix+'_Aligned.sortedByCoord.out.bam'

	cmd = 'telescope assign --reassign_mode average --outdir ./' \
		+' --attribute transcript_id '+align_bam+' '+ref_TE_gtf

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	TE_dict = defaultdict()
	fho = open(unit_quant_out,'w')
	print("%s\t%s" % ('TE_ID', 'count'), file=fho, flush=True)
	with open('telescope-telescope_report.tsv','r') as fh:
		for line in fh.readlines():
			line = line.strip().split('\t')
			if line[0][0]=='#': continue
			if line[0]=='transcript' or line[0]=='__no_feature': continue

			ID = line[0]
			count = float(line[2])

			if count > 0:
				TE_dict[ID] = count
				print("%s\t%d" % (ID, count), file=fho, flush=True)
	fho.close()

	family_dict = defaultdict()
	with open(ref_TE_gtf,'r') as fh:
		for line in fh.readlines():
			line = line.strip().split('\t')

			if line[0][0]=='#': continue

			attributes = defaultdict()
			for a in line[8].replace('"', '').split(';')[:-1]:
				kv = a.strip().split(' ')
				attributes[kv[0]] = kv[1]

			TE_id = attributes['transcript_id']
			TE_name = attributes['gene_id']
			family_id = TE_id.split('_')[0]

			if TE_id in TE_dict:
				if family_id not in family_dict:
					family_dict[family_id] = TE_dict[TE_id]
				else:
					family_dict[family_id] = family_dict[family_id] + TE_dict[TE_id]

	fho = open(family_quant_out,'w')
	print("%s\t%s" % ('family_ID', 'count'), file=fho, flush=True)
	for family_id in family_dict:
		print("%s\t%d" % (family_id, family_dict[family_id]), file=fho, flush=True)
	fho.close()

def TEEN(args):
	TE_exon = pd.read_table(TE_exon_bed, header=None)
	TE_exon.columns = ['chr', 'start', 'end', 'transcript_id', 'gene_id', 'strand', 'exon_type', 'exon_class']
	TE_exon['exon_ID'] = TE_exon['chr']+':'+TE_exon['start'].astype('str')+'-'+TE_exon['end'].astype('str')+'('+TE_exon['strand']+')'
	TE_exon['exon_length'] = TE_exon['end'] - TE_exon['start']

	quant = pd.read_table(transcript_quant_out)

	if args.quant == 'kallisto':
		quant = quant[['target_id', 'length', 'eff_length', 'est_counts', 'tpm']]
	elif args.quant == 'rsem':
		quant = quant[['transcript_id', 'length', 'effective_length', 'expected_count', 'TPM']]

	quant.columns = ['transcript_id', 'length', 'eff_length', 'count', 'TPM']
	transcript_quant = quant
	transcript_quant['count'] = round(transcript_quant['count'],2)
	transcript_quant['TPM'] = round(transcript_quant['TPM'],2)
	transcript_quant.to_csv(transcript_quant_out, sep='\t', index=0)

	exon_quant = pd.merge(TE_exon, quant, on='transcript_id')
	exon_quant['count'] = exon_quant['exon_length']*exon_quant['count']/exon_quant['eff_length']
	exon_count = exon_quant.groupby('exon_class')['count'].apply(sum).reset_index()
	exon_count['count'] = round(exon_count['count'],2)
	exon_TPM = exon_quant.groupby('exon_class')['TPM'].apply(sum).reset_index()
	exon_TPM['TPM'] = round(exon_TPM['TPM'],2)
	exon_cluster_quant = pd.merge(exon_count, exon_TPM, on='exon_class')
	exon_cluster_quant['exon_cluster'] = 'exon_cluster_'+exon_cluster_quant['exon_class'].astype('str')
	exon_cluster_quant[['exon_cluster', 'count', 'TPM']].to_csv(exon_quant_out, sep='\t', index=0)
	os.remove(TE_exon_self_overlap)

parser = argparse.ArgumentParser(description='TEEN: Transposable Element Expression')
parser.add_argument('-fq1', '--fastq1', help='Read1 in FASTQ format (required)')
parser.add_argument('-fq2', '--fastq2', help='Read1 in FASTQ format (required)')
parser.add_argument('--TE_gtf', help='TE position in GTF format (required)')
parser.add_argument('-l', '--level', default=1, help='TE level: 1 for exon-level, 2 for unit-level and family-level (default: 1)')
parser.add_argument('-s', '--stranded_type', help='Strand-specific RNA-seq read orientation: RF or FR')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory (default: .)')
parser.add_argument('-p', '--prefix', default='TEEN', help='Prefix for output file name (default: TEEN)')
parser.add_argument('-r', '--ref_genome', help='Reference genome in FASTA format (required)')
parser.add_argument('-a', '--annotation', help='Genome annotation in GTF format (required)')
parser.add_argument('--TE_exon', help='TE exon annotation in BED format')
parser.add_argument('-q', '--quant', help='Quantification pattern: kallisto or rsem (required for exon-level quantification)')
parser.add_argument('-i', '--index', default='./TEEN_index/TEEN', help='Index name (default: ./TEEN_index/TEEN)')
parser.add_argument('-t', '--nthread', type=int, default=1, help='Number of threads to run TEEN (default: 1)')

args = parser.parse_args()
out_dir = os.path.abspath(args.output_dir)
Index = os.path.abspath(args.index)
ref_dir = os.path.abspath(os.path.dirname(args.index))

ref_fasta = os.path.abspath(args.ref_genome)
ref_gtf = os.path.abspath(args.annotation)
ref_TE_gtf = os.path.abspath(args.TE_gtf)

fastq1 = os.path.abspath(args.fastq1)
fastq2 = os.path.abspath(args.fastq2)

transcript_quant_out = out_dir+'/'+args.prefix+'.transcript.quant.out'
exon_quant_out = out_dir+'/'+args.prefix+'.TE.exon.quant.out'
unit_quant_out = out_dir+'/'+args.prefix+'.TE.unit.quant.out'
family_quant_out = out_dir+'/'+args.prefix+'.TE.family.quant.out'

if not os.path.exists(out_dir):
    os.makedirs(out_dir)


print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running TEEN on {0:d} threads.'.format(args.nthread), flush=True)


with cd(out_dir):
	if int(args.level) == 2:
		if not os.path.exists(os.path.abspath(os.path.dirname(Index))):
			os.makedirs(os.path.abspath(os.path.dirname(Index)))

			print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to create index.', flush=True)

			if not args.annotation:
				print('ERROR: Lack annotation file (--annotation)')
			if not args.ref_genome:
				print('ERROR: Lack reference genome (--ref_genome)')

			create_STAR_index(args)

			print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish.', flush=True)

		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to do quantification.', flush=True)

		STAR(args)
		
		telescope(args)

	elif args.quant == 'kallisto':
		if not os.path.exists(Index):
			if not os.path.exists(os.path.abspath(os.path.dirname(Index))):
				os.makedirs(os.path.abspath(os.path.dirname(Index)))

			print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to create index.', flush=True)

			if not args.annotation:
				print('ERROR: Lack annotation file (--annotation)')
			if not args.ref_genome:
				print('ERROR: Lack reference genome (--ref_genome)')

			transcript_fasta = os.path.abspath(os.path.dirname(Index))+'/'+args.prefix+'.transcripts.fa'
			create_kallisto_index(args)

			print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish.', flush=True)

		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to do quantification.', flush=True)

		kallisto(args)
		shutil.copyfile('abundance.tsv', transcript_quant_out)
		cmd = 'rm -rf abundance.h5 abundance.tsv run_info.json'
		subprocess.check_call(cmd, shell=True, executable='/bin/bash')
		TEEN(args)

	elif args.quant == 'rsem':
		if not os.path.exists(os.path.abspath(os.path.dirname(Index))):
			os.makedirs(os.path.abspath(os.path.dirname(Index)))

			print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to create index.', flush=True)

			if not args.annotation:
				print('ERROR: Lack annotation file (--annotation)')
			if not args.ref_genome:
				print('ERROR: Lack reference genome (--ref_genome)')

			create_rsem_index(args)

			print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish.', flush=True)

		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to do quantification.', flush=True)

		rsem(args)
		shutil.copyfile(args.prefix+'.isoforms.results', transcript_quant_out)
		cmd = 'rm -rf '+args.prefix+'.isoforms.results '+args.prefix+'.genes.results '+args.prefix+'.log '+args.prefix+'.stat '+args.prefix+'.transcript.bam'
		subprocess.check_call(cmd, shell=True, executable='/bin/bash')
		TEEN(args)

	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish.', flush=True)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] TE quantification is done.', flush=True)
