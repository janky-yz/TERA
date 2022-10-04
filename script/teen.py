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
	cmd = 'gffread'+' -g '+ref_fasta \
		+' -w '+transcript_fasta \
		+' '+ref_gtf

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
	cmd = 'rsem-prepare-reference'+' -p '+str(args.nthread) \
		+' --star --gtf '+ref_gtf \
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

def extract_TE_exon(args):
	exon_dict = {}
	with open(ref_gtf,'r') as fh:
		for line in fh.readlines():
			line = line.strip().split('\t')

			if line[0][0]=='#': continue # skip header

			chrn = line[0]
			entry_type = line[2]
			exon = (int(line[3]),int(line[4]))
			strand = line[6]

			attributes = defaultdict()
			for a in line[8].replace('"', '').split(';')[:-1]:
				kv = a.strip().split(' ')
				attributes[kv[0]] = kv[1]

			if entry_type == 'exon':
				transcript_id = attributes['transcript_id']
				gene_id = attributes['gene_id']
				tid = chrn+":"+transcript_id+":"+gene_id+":"+strand
				if tid not in exon_dict:
					exon_dict[tid] = [exon]
				else:
					exon_dict[tid].append(exon)

	fho = open(transcript_exon_bed,'w')
	for tid in exon_dict:
		chrn,transcript_id,gene_id,strand = tid.split(":")
		exon_list = exon_dict[tid]
		if len(exon_list) > 1:
			for i in range(len(exon_list)):
				if (i == 0 and strand == '+') or (i == len(exon_list)-1 and strand =='-'):
					print("%s\t%d\t%d\t%s\t%s\t%s\t%s" % (chrn, int(exon_list[i][0])-1, int(exon_list[i][1]), transcript_id, gene_id, strand, 'start'), file=fho, flush=True)
				elif (i == 0 and strand == '-') or (i == len(exon_list)-1 and strand =='+'):
					print("%s\t%d\t%d\t%s\t%s\t%s\t%s" % (chrn, int(exon_list[i][0])-1, int(exon_list[i][1]), transcript_id, gene_id, strand, 'end'), file=fho, flush=True)
				else:
					print("%s\t%d\t%d\t%s\t%s\t%s\t%s" % (chrn, int(exon_list[i][0])-1, int(exon_list[i][1]), transcript_id, gene_id, strand, 'middle'), file=fho, flush=True)
		else:
			print("%s\t%d\t%d\t%s\t%s\t%s\t%s" % (chrn, int(exon_list[0][0])-1, int(exon_list[0][1]), transcript_id, gene_id, strand, 'SE'), file=fho, flush=True)
	fho.close()

	cmd = "bedtools intersect -a "+transcript_exon_bed+" -b "+ref_TE_bed+" -s -wo | awk '{if($NF>20){print $0}}' | cut -f 1-7 | uniq >"+TE_exon_bed
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def exon_compare(exon_overlap, cutoff):
	keep_middle = (exon_overlap[6]=="middle") & (abs(exon_overlap[1]-exon_overlap[8])<=cutoff) & (abs(exon_overlap[2]-exon_overlap[9])<=cutoff)
	keep_start_rev = (exon_overlap[6]=="start") & (abs(exon_overlap[1]-exon_overlap[8])<=cutoff) & (exon_overlap[5]=="-")
	keep_start_for = (exon_overlap[6]=="start") & (abs(exon_overlap[2]-exon_overlap[9])<=cutoff) & (exon_overlap[5]=="+")
	keep_end_rev = (exon_overlap[6]=="end") & (abs(exon_overlap[2]-exon_overlap[9])<=cutoff) & (exon_overlap[5]=="-")
	keep_end_for = (exon_overlap[6]=="end") & (abs(exon_overlap[1]-exon_overlap[8])<=cutoff) & (exon_overlap[5]=="+")
	keep_SE = exon_overlap[6]=="SE"
	keep = keep_middle | keep_start_rev | keep_start_for | keep_end_rev | keep_end_for | keep_SE

	return(keep)

def TEEN(args):
	cmd = "bedtools intersect -a "+TE_exon_bed+" -b "+TE_exon_bed+" -s -wo | awk '{if($NF>20 && $7==$14){print $0}}' >"+TE_exon_self_overlap
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	exon_overlap = pd.read_table(TE_exon_self_overlap, header=None)
	exon_overlap['ID1'] = exon_overlap[0]+':'+exon_overlap[1].astype('str')+'-'+exon_overlap[2].astype('str')+'('+exon_overlap[5]+')'
	exon_overlap['ID2'] = exon_overlap[7]+':'+exon_overlap[8].astype('str')+'-'+exon_overlap[9].astype('str')+'('+exon_overlap[12]+')'

	exon_filt = exon_overlap[exon_overlap['ID1'] != exon_overlap['ID2']]
	keep_exon = exon_compare(exon_filt, args.exon_diff)
	exon_filt = exon_filt[keep_exon]

	exon_filt['exon_ID1'] = exon_filt[['ID1', 'ID2']].apply(min, axis=1)
	exon_filt['exon_ID2'] = exon_filt[['ID1', 'ID2']].apply(max, axis=1)
	ID_match = exon_filt[['exon_ID1', 'exon_ID2']].drop_duplicates().reset_index()

	exon_ID = []
	exon_class = []
	index = 0
	for i in range(ID_match.shape[0]):
		if((ID_match.loc[i,'exon_ID1'] not in exon_ID) & (ID_match.loc[i,'exon_ID2'] not in exon_ID)):
			exon_ID += [ID_match.loc[i,'exon_ID1'], ID_match.loc[i,'exon_ID2']]
			index += 1
			exon_class += [index, index]
		if((ID_match.loc[i,'exon_ID1'] in exon_ID) & (ID_match.loc[i,'exon_ID2'] not in exon_ID)):
			exon_ID.append(ID_match.loc[i,'exon_ID2'])
			exon_class.append(exon_class[exon_ID.index(ID_match.loc[i,'exon_ID1'])])
		if((ID_match.loc[i,'exon_ID1'] not in exon_ID) & (ID_match.loc[i,'exon_ID2'] in exon_ID)):
			exon_ID.append(ID_match.loc[i,'exon_ID1'])
			exon_class.append(exon_class[exon_ID.index(ID_match.loc[i,'exon_ID2'])])

	cluster1 = pd.DataFrame({'exon_ID':exon_ID, 'exon_class':exon_class})
	all_ID = exon_overlap['ID1'].drop_duplicates()
	remain_ID = sorted(list(set(all_ID).difference(set(cluster1['exon_ID']))))
	cluster2 = pd.DataFrame({'exon_ID':remain_ID, 'exon_class':np.arange(index+1,index+1+len(remain_ID))})
	exon_cluster = pd.concat([cluster1, cluster2])

	TE_exon = pd.read_table(TE_exon_bed, header=None)
	TE_exon.columns = ['chr', 'start', 'end', 'transcript_id', 'gene_id', 'strand', 'exon_type']
	TE_exon['exon_ID'] = TE_exon['chr']+':'+TE_exon['start'].astype('str')+'-'+TE_exon['end'].astype('str')+'('+TE_exon['strand']+')'
	TE_exon['exon_length'] = TE_exon['end'] - TE_exon['start']
	TE_exon_cluster = pd.merge(TE_exon, exon_cluster, on='exon_ID')

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

	exon_quant = pd.merge(TE_exon_cluster, quant, on='transcript_id')
	exon_quant['count'] = exon_quant['exon_length']*exon_quant['count']/exon_quant['eff_length']
	exon_count = exon_quant.groupby('exon_class')['count'].apply(sum).reset_index()
	exon_count['count'] = round(exon_count['count'],2)
	exon_TPM = exon_quant.groupby('exon_class')['TPM'].apply(sum).reset_index()
	exon_TPM['TPM'] = round(exon_TPM['TPM'],2)
	exon_cluster_quant = pd.merge(exon_count, exon_TPM, on='exon_class')
	exon_cluster_quant['exon_cluster'] = 'exon_cluster_'+exon_cluster_quant['exon_class'].astype('str')
	exon_cluster_quant[['exon_cluster', 'count', 'TPM']].to_csv(exon_quant_out, sep='\t', index=0)
	os.remove(TE_exon_self_overlap)

	cmd = "bedtools intersect -a "+TE_exon_bed+" -b "+ref_TE_bed+" -s -wo | awk '{if($NF>20){print $0}}' >"+TE_exon_ref_overlap
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	TE_exon_anno = pd.read_table(TE_exon_ref_overlap, header=None)
	TE_exon_anno['ID'] = TE_exon_anno[0]+':'+TE_exon_anno[1].astype('str')+'-'+TE_exon_anno[2].astype('str')+'('+TE_exon_anno[5]+'):'+TE_exon_anno[3]
	TE_ID = TE_exon_anno.groupby('ID')[10].apply(lambda x: ','.join(x)).reset_index()
	TE_name = TE_exon_anno.groupby('ID')[11].apply(lambda x: ','.join(x)).reset_index()
	TE_family = TE_exon_anno.groupby('ID')[13].apply(lambda x: ','.join(x)).reset_index()
	TE_class = TE_exon_anno.groupby('ID')[14].apply(lambda x: ','.join(x)).reset_index()
	TE_anno = pd.merge(TE_ID, TE_name, on='ID')
	TE_anno = pd.merge(TE_anno, TE_family, on='ID')
	TE_anno = pd.merge(TE_anno, TE_class, on='ID')
	TE_anno.columns = ['ID', 'TE_ID', 'TE_name', 'TE_family', 'TE_class']

	TE_exon_cluster['ID'] = TE_exon_cluster['chr']+':'+TE_exon_cluster['start'].astype('str')+'-'+TE_exon_cluster['end'].astype('str')+'('+TE_exon_cluster['strand']+'):'+TE_exon_cluster['transcript_id']
	TE_exon_anno = pd.merge(TE_exon_cluster, TE_anno, on='ID')
	TE_exon_anno['exon_cluster'] = 'exon_cluster_'+TE_exon_anno['exon_class'].astype('str')
	TE_exon_anno = TE_exon_anno[['chr', 'start', 'end', 'transcript_id', 'gene_id', 'strand', 'exon_type', 'exon_cluster', 'TE_ID','TE_name','TE_family','TE_class']]
	TE_exon_anno = TE_exon_anno.sort_values(by=['chr', 'start', 'end'])
	TE_exon_anno.to_csv(TE_exon_bed, sep='\t', header=0, index=0)
	os.remove(TE_exon_ref_overlap)


parser = argparse.ArgumentParser(description='TEEN: Transposable Element Expression')
parser.add_argument('-fq1', '--fastq1', help='Read1 in FASTQ format (required)')
parser.add_argument('-fq2', '--fastq2', help='Read1 in FASTQ format (required)')
parser.add_argument('-e', '--TE', help='TE position in BED format (required)')
parser.add_argument('-s', '--stranded_type', help='Strand-specific RNA-seq read orientation: RF or FR')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory (default: .)')
parser.add_argument('-p', '--prefix', default='TEEN', help='Prefix for output file name (default: TEEN)')
parser.add_argument('-r', '--ref_genome', help='Reference genome in FASTA format (required)')
parser.add_argument('-a', '--annotation', help='Genome annotation in GTF format (required)')
parser.add_argument('--TE_exon', help='TE exon annotation in BED format')
parser.add_argument('-q', '--quant', help='Quantification pattern: kallisto or rsem (required)')
parser.add_argument('-i', '--index', default='./TEEN_index/TEEN', help='Index name (default: ./TEEN_index/TEEN)')
parser.add_argument('-t', '--nthread', type=int, default=1, help='Number of threads to run TEEN (default: 1)')
parser.add_argument('-d', '--exon_diff', type=int, default=10, help='Maximum difference (bp) of exon ends (default: 10)')

args = parser.parse_args()
out_dir = os.path.abspath(args.output_dir)
Index = os.path.abspath(args.index)
ref_rsem = os.path.abspath(os.path.dirname(args.index))

ref_fasta = os.path.abspath(args.ref_genome)
ref_gtf = os.path.abspath(args.annotation)
ref_TE_bed = os.path.abspath(args.TE)

fastq1 = os.path.abspath(args.fastq1)
fastq2 = os.path.abspath(args.fastq2)

TE_exon_self_overlap = args.prefix+'_TE_exon_self_overlap.txt'
TE_exon_ref_overlap = args.prefix+'_TE_exon_ref_overlap.txt'

TE_exon_bed = out_dir+'/'+args.prefix+'_TE_exon_anno.bed'
transcript_quant_out = out_dir+'/'+args.prefix+'_transcript_quant.out'
exon_quant_out = out_dir+'/'+args.prefix+'_TE_exon_quant.out'

if not os.path.exists(out_dir):
    os.makedirs(out_dir)


print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running TEEN on {0:d} threads.'.format(args.nthread), flush=True)


with cd(out_dir):
	if args.quant == 'kallisto':
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

	if args.TE_exon:
		cmd = 'cut -f1-7 '+os.path.abspath(args.TE_exon)+' >'+TE_exon_bed
		subprocess.check_call(cmd, shell=True, executable='/bin/bash')
	else:
		transcript_exon_bed = out_dir+'/'+args.prefix+'_transcript_exon.bed'
		extract_TE_exon(args)
		os.remove(transcript_exon_bed)

	TEEN(args)

	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish.', flush=True)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] TE quantification is done.', flush=True)
