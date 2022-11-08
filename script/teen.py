#!/usr/bin/env python3
# Author: Janky

import argparse
import os
import subprocess
import gzip
import multiprocessing as mp
import pandas as pd
import numpy as np
import csv
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


def kallisto(args, TEEN_index, fastq1, fastq2):
	cmd = 'kallisto quant'+' -i '+TEEN_index+'/TEEN.index' \
		+' -o '+out_dir \
		+' -t '+str(args.nthread)

	if args.stranded_type:
		cmd += ' --'+args.stranded_type.lower()+'-stranded'

	cmd += ' ' + fastq1 + ' ' + fastq2

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def rsem(args, TEEN_index, fastq1, fastq2):
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

	cmd += ' '+TEEN_index+'/TEEN'+' '+args.prefix

	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def extract_TE_exon(input_gtf, ref_TE_bed, transcript_exon_bed, TE_exon_bed):
	exon_dict = defaultdict()
	with open(input_gtf,'r') as fh:
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

def gtf_deal(ref_TE_bed, ref_gtf, TE_gtf, TEEN_gtf, TEEN_TE_exon_bed):
	ref_sort_gtf = 'ref_sort.gtf'
	ref_transcript_exon_bed = 'ref_transcript_exon.bed'
	ref_TE_exon_bed = 'ref_TE_exon.bed'
	ref_nonTE_gtf = 'ref_nonTE.gtf'

	cmp_gtf = 'gffcmp.annotated.gtf'
	TE_cmp_gtf = 'TE_cmp.gtf'
	merge_gtf = 'ref_TE_merge.gtf'

	cmd = 'gffread'+' -T --sort-alpha -o'+ref_sort_gtf \
		+' '+ref_gtf
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	extract_TE_exon(ref_sort_gtf, ref_TE_bed, ref_transcript_exon_bed, ref_TE_exon_bed)
	gtf = pd.read_table(ref_sort_gtf, header=None)
	gtf['transcript_id'] = gtf.apply(lambda x:x[8].split('"')[1], axis=1)
	gtf['gene_id'] = gtf.apply(lambda x:x[8].split('"')[3], axis=1)
	bed = pd.read_table(ref_TE_exon_bed, header=None)
	gtf_nonTE = gtf[~gtf.transcript_id.isin(bed[3].drop_duplicates())]
	gtf_nonTE.iloc[:,:9].to_csv(ref_nonTE_gtf, quoting=csv.QUOTE_NONE, sep='\t', header=0, index=0)
	gene_map = gtf_nonTE[['transcript_id', 'gene_id']].drop_duplicates()

	cmd = 'gffcompare'+' -T -r '+ref_nonTE_gtf \
		+' '+TE_gtf
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	gene_dict = defaultdict()

	fho = open(TE_cmp_gtf, 'w')
	with open(cmp_gtf,'r') as fh:
		for line in fh.readlines():
			line = line.strip().split('\t')

			entry_type = line[2]

			attributes = defaultdict()
			for a in line[8].replace('"', '').split(';')[:-1]:
				kv = a.strip().split(' ')
				attributes[kv[0]] = kv[1]

			if entry_type == 'transcript':
				transcript_id = attributes['transcript_id']
				class_code = attributes['class_code']
				if class_code in ['=','c','j','e','i','o']:
					cmp_ref = attributes['cmp_ref']
					gene_id = gene_map.loc[gene_map.transcript_id==cmp_ref,'gene_id'].to_list()[0]
				else:
					gene_id = attributes['gene_id']

				gene_dict[transcript_id] = gene_id
			elif entry_type == 'exon':
				transcript_id = attributes['transcript_id']
				gene_id = gene_dict[transcript_id]
			attr = 'transcript_id '+'"'+transcript_id+'";'+' gene_id '+'"'+gene_id+'";'
			print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], attr), file=fho, flush=True)

	cmd = 'cat '+ref_nonTE_gtf+' '+TE_cmp_gtf+' >'+merge_gtf
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'gffread'+' -T --sort-alpha -o '+TEEN_gtf+' '+merge_gtf
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	TE_exon = pd.read_table(TEEN_TE_exon_bed, header=None)
	TE_exon[4] = TE_exon.apply(lambda x:gene_dict[x[3]], axis=1)
	TE_exon.to_csv(TEEN_TE_exon_bed, header=0, index=0, sep='\t')

	subprocess.check_call('rm *cmp* ref*', shell=True, executable='/bin/bash')

def exon_compare(exon_overlap, cutoff):
	keep_middle = (exon_overlap[6]=="middle") & (abs(exon_overlap[1]-exon_overlap[8])<=cutoff) & (abs(exon_overlap[2]-exon_overlap[9])<=cutoff)
	keep_start_rev = (exon_overlap[6]=="start") & (abs(exon_overlap[1]-exon_overlap[8])<=cutoff) & (exon_overlap[5]=="-")
	keep_start_for = (exon_overlap[6]=="start") & (abs(exon_overlap[2]-exon_overlap[9])<=cutoff) & (exon_overlap[5]=="+")
	keep_end_rev = (exon_overlap[6]=="end") & (abs(exon_overlap[2]-exon_overlap[9])<=cutoff) & (exon_overlap[5]=="-")
	keep_end_for = (exon_overlap[6]=="end") & (abs(exon_overlap[1]-exon_overlap[8])<=cutoff) & (exon_overlap[5]=="+")
	keep_SE = exon_overlap[6]=="SE"
	keep = keep_middle | keep_start_rev | keep_start_for | keep_end_rev | keep_end_for | keep_SE

	return(keep)

def TE_cluster(ref_TE_bed, TEEN_TE_exon_bed, TEEN_TE_exon_cluster_out, TEEN_TE_exon_anno_bed, cutoff):
	TE_exon_self_overlap = 'TEEN_TE_exon_self_overlap.txt'
	TE_exon_ref_overlap = 'TEEN_TE_exon_ref_overlap.txt'
	TE_exon_cluster_out = 'TEEN_TE_exon_cluster.out'

	cmd = "bedtools intersect -a "+TEEN_TE_exon_bed+" -b "+TEEN_TE_exon_bed+" -s -wo | awk '{if($NF>20 && $7==$14){print $0}}' >"+TE_exon_self_overlap
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	exon_overlap = pd.read_table(TE_exon_self_overlap, header=None)
	exon_overlap['ID1'] = exon_overlap[0]+':'+exon_overlap[1].astype('str')+'-'+exon_overlap[2].astype('str')+'('+exon_overlap[5]+')'
	exon_overlap['ID2'] = exon_overlap[7]+':'+exon_overlap[8].astype('str')+'-'+exon_overlap[9].astype('str')+'('+exon_overlap[12]+')'

	exon_filt = exon_overlap[exon_overlap['ID1'] != exon_overlap['ID2']]
	keep_exon = exon_compare(exon_filt, cutoff)
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

	TE_exon = pd.read_table(TEEN_TE_exon_bed, header=None)
	TE_exon.columns = ['chr', 'start', 'end', 'transcript_id', 'gene_id', 'strand', 'exon_type']
	TE_exon['exon_ID'] = TE_exon['chr']+':'+TE_exon['start'].astype('str')+'-'+TE_exon['end'].astype('str')+'('+TE_exon['strand']+')'
	TE_exon['exon_length'] = TE_exon['end'] - TE_exon['start']
	TE_exon_cluster = pd.merge(TE_exon, exon_cluster, on='exon_ID')
	TE_exon_cluster.to_csv(TEEN_TE_exon_cluster_out, sep='\t', index=0)

	cmd = "bedtools intersect -a "+TEEN_TE_exon_bed+" -b "+ref_TE_bed+" -s -wo | awk '{if($NF>20){print $0}}' >"+TE_exon_ref_overlap
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
	TE_exon_anno.to_csv(TEEN_TE_exon_anno_bed, sep='\t', header=0, index=0)
	subprocess.check_call('rm *overlap.txt', shell=True, executable='/bin/bash')

def index(args):
	ref_fasta = os.path.abspath(args.ref_genome)
	ref_TE_bed = os.path.abspath(args.ref_TE)
	ref_gtf = os.path.abspath(args.ref_anno)
	TE_gtf = os.path.abspath(args.TE_anno)

	TEEN_gtf = out_dir+'/TEEN.gtf'
	TEEN_fasta = out_dir+'/TEEN.fasta'
	TEEN_TE_exon_bed = out_dir+'/TEEN_TE_exon.bed'
	TEEN_TE_exon_cluster_out = out_dir+'/TEEN_TE_exon_cluster.out'
	TEEN_TE_exon_anno_bed = out_dir+'/TEEN_TE_exon_anno.bed'

	if args.TE_exon:
		cmd = 'cut -f1-7 '+os.path.abspath(args.TE_exon)+' >'+TEEN_TE_exon_bed
		subprocess.check_call(cmd, shell=True, executable='/bin/bash')
	else:
		transcript_exon_bed = out_dir+'/TEEN_transcript_exon.bed'
		extract_TE_exon(TE_gtf, ref_TE_bed, transcript_exon_bed, TEEN_TE_exon_bed)
		os.remove(transcript_exon_bed)

	with cd(out_dir):
		gtf_deal(ref_TE_bed, ref_gtf, TE_gtf, TEEN_gtf, TEEN_TE_exon_bed)
		TE_cluster(ref_TE_bed, TEEN_TE_exon_bed, TEEN_TE_exon_cluster_out, TEEN_TE_exon_anno_bed, args.exon_diff)

		cmd = 'gffread'+' -g '+ref_fasta \
			+' -w '+TEEN_fasta \
			+' '+TEEN_gtf

		subprocess.check_call(cmd, shell=True, executable='/bin/bash')

		if args.rsem:
			cmd = 'rsem-prepare-reference'+' -p '+str(args.nthread) \
				+' --star --gtf '+TEEN_gtf \
				+' '+ref_fasta \
				+' TEEN'
		else:
			cmd = 'kallisto index'+' -i TEEN.index' \
				+' '+TEEN_fasta

		subprocess.check_call(cmd, shell=True, executable='/bin/bash')

def quant(args):
	TEEN_index = os.path.abspath(args.index)
	fastq1 = os.path.abspath(args.fastq1)
	fastq2 = os.path.abspath(args.fastq2)

	TEEN_TE_exon_cluster_out = TEEN_index+'/TEEN_TE_exon_cluster.out'
	transcript_quant_out = out_dir+'/'+args.prefix+'_transcript_quant.out'
	exon_quant_out = out_dir+'/'+args.prefix+'_TE_exon_quant.out'

	with cd(out_dir):
		if args.rsem:
			rsem(args, TEEN_index, fastq1, fastq2)
			shutil.copyfile(args.prefix+'.isoforms.results', transcript_quant_out)
			cmd = 'rm -rf '+args.prefix+'.isoforms.results '+args.prefix+'.genes.results '+args.prefix+'.log '+args.prefix+'.stat '+args.prefix+'.transcript.bam'
			subprocess.check_call(cmd, shell=True, executable='/bin/bash')

			quant = pd.read_table(transcript_quant_out)
			quant = quant[['transcript_id', 'length', 'effective_length', 'expected_count', 'TPM']]
		else:
			kallisto(args, TEEN_index, fastq1, fastq2)
			shutil.copyfile('abundance.tsv', transcript_quant_out)
			cmd = 'rm -rf abundance.h5 abundance.tsv run_info.json'
			subprocess.check_call(cmd, shell=True, executable='/bin/bash')

			quant = pd.read_table(transcript_quant_out)
			quant = quant[['target_id', 'length', 'eff_length', 'est_counts', 'tpm']]

		quant.columns = ['transcript_id', 'length', 'eff_length', 'count', 'TPM']
		transcript_quant = quant
		transcript_quant['count'] = round(transcript_quant['count'],2)
		transcript_quant['TPM'] = round(transcript_quant['TPM'],2)
		transcript_quant.to_csv(transcript_quant_out, sep='\t', index=0)

		TE_exon_cluster = pd.read_table(TEEN_TE_exon_cluster_out)
		exon_quant = pd.merge(TE_exon_cluster, quant, on='transcript_id')
		exon_quant['count'] = exon_quant['exon_length']*exon_quant['count']/exon_quant['eff_length']
		exon_count = exon_quant.groupby('exon_class')['count'].apply(sum).reset_index()
		exon_count['count'] = round(exon_count['count'],2)
		exon_TPM = exon_quant.groupby('exon_class')['TPM'].apply(sum).reset_index()
		exon_TPM['TPM'] = round(exon_TPM['TPM'],2)
		exon_cluster_quant = pd.merge(exon_count, exon_TPM, on='exon_class')
		exon_cluster_quant['exon_cluster'] = 'exon_cluster_'+exon_cluster_quant['exon_class'].astype('str')
		exon_cluster_quant[['exon_cluster', 'count', 'TPM']].to_csv(exon_quant_out, sep='\t', index=0)


parser = argparse.ArgumentParser(description='TEEN: Transposable Element Expression')
subparsers = parser.add_subparsers(help='sub-command help')

parser_index = subparsers.add_parser('index', help='index help')
parser_index.add_argument('-r', '--ref_genome', help='Reference genome in FASTA format (required)')
parser_index.add_argument('-e', '--ref_TE', help='Reference TE position in BED format (required)')
parser_index.add_argument('--ref_anno', help='Reference genome annotation in GTF format (required)')
parser_index.add_argument('--TE_anno', help='TE annotation in GTF format (required)')
parser_index.add_argument('--TE_exon', help='TE exon annotation in BED format')
parser_index.add_argument('--kallisto', help='Specific Quantification by kallisto', action='store_true')
parser_index.add_argument('--rsem', help='Specific Quantification by RSEM', action='store_true')
parser_index.add_argument('-o', '--output_dir', default='./TEEN_index/', help='Output directory (default: ./TEEN_index/)')
parser_index.add_argument('-t', '--nthread', type=int, default=1, help='Number of threads (default: 1)')
parser_index.add_argument('-d', '--exon_diff', type=int, default=10, help='Maximum difference (bp) of exon ends (default: 10)')

parser_index.set_defaults(func=index)

parser_quant = subparsers.add_parser('quant', help='quant help')
parser_quant.add_argument('-fq1', '--fastq1', help='Read1 in FASTQ format (required)')
parser_quant.add_argument('-fq2', '--fastq2', help='Read1 in FASTQ format (required)')
parser_quant.add_argument('-s', '--stranded_type', help='Strand-specific RNA-seq read orientation: RF or FR')
parser_quant.add_argument('-o', '--output_dir', default='./TEEN_quant', help='Output directory (default: ./TEEN_quant)')
parser_quant.add_argument('-p', '--prefix', default='TEEN', help='Prefix for output file name (default: TEEN)')
parser_quant.add_argument('--kallisto', help='Specific Quantification by kallisto', action='store_true')
parser_quant.add_argument('--rsem', help='Specific Quantification by RSEM', action='store_true')
parser_quant.add_argument('-i', '--index', default='./TEEN_index/', help='Index directory (default: ./TEEN_index/)')
parser_quant.add_argument('-t', '--nthread', type=int, default=1, help='Number of threads to run TEEN (default: 1)')

parser_quant.set_defaults(func=quant)

args = parser.parse_args()
out_dir = os.path.abspath(args.output_dir)

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running TEEN on {0:d} threads.'.format(args.nthread), flush=True)

args.func(args)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Done.', flush=True)
