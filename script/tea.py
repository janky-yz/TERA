#!/usr/bin/env python3
# Author: Janky

import argparse
import os
import csv
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

def interval_union(intervals):
    intervals.sort(key=lambda x: x[0])
    union = [intervals[0]]
    for i in intervals[1:]:
        if i[0] <= union[-1][1]:
            if i[1] > union[-1][1]:
                union[-1][1] = i[1]
        else:
            union.append(i)
    return union

def deal_ref_anno(args):
	exon_dict = {}
	with open(gene_gtf,'r') as fh:
		for line in fh.readlines():
			line = line.strip().split('\t')

			if line[0][0]=='#': continue

			chrn = line[0]
			entry_type = line[2]
			exon = [int(line[3]),int(line[4])]
			strand = line[6]

			attributes = defaultdict()
			for a in line[8].replace('"', '').split(';')[:-1]:
				kv = a.strip().split(' ')
				attributes[kv[0]] = kv[1]

			if entry_type == 'exon':
				transcript_id = attributes['transcript_id']
				gene_id = attributes['gene_id']
				gid = chrn+":"+gene_id+":"+strand
				if gid not in exon_dict:
					exon_dict[gid] = [exon]
				else:
					exon_dict[gid].append(exon)

	fho1 = open(gene_merged_exon_bed,'w')
	fho2 = open(gene_bed, 'w')
	for gid in exon_dict:
		chrn,gene_id,strand = gid.split(":")
		exon_list = interval_union(exon_dict[gid])
		print("%s\t%d\t%d\t%s\t%s\t%s" % (chrn, int(exon_list[0][0])-1, int(exon_list[len(exon_list)-1][1]), gene_id, gene_id, strand), file=fho2, flush=True)
		for i in range(len(exon_list)):
			print("%s\t%d\t%d\t%s\t%s\t%s" % (chrn, int(exon_list[i][0])-1, int(exon_list[i][1]), gene_id, gene_id, strand), file=fho1, flush=True)
	fho1.close()
	fho2.close()

def extract_TE_exon(args):
	exon_dict = {}
	with open(input_gtf,'r') as fh:
		fho = open(TE_transcript_bed, 'w')

		for line in fh.readlines():
			line = line.strip().split('\t')

			if line[0][0]=='#': continue

			chrn = line[0]
			entry_type = line[2]
			exon = (int(line[3]),int(line[4]))
			strand = line[6]

			attributes = defaultdict()
			for a in line[8].replace('"', '').split(';')[:-1]:
				kv = a.strip().split(' ')
				attributes[kv[0]] = kv[1]

			if entry_type == 'transcript':
				transcript_id = attributes['transcript_id']
				gene_id = attributes['gene_id']
				print("%s\t%d\t%d\t%s\t%s\t%s" % (chrn, int(line[3])-1, int(line[4]), transcript_id, gene_id, strand), file=fho, flush=True)

			if entry_type == 'exon':
				transcript_id = attributes['transcript_id']
				gene_id = attributes['gene_id']
				tid = chrn+":"+transcript_id+":"+gene_id+":"+strand
				if tid not in exon_dict:
					exon_dict[tid] = [exon]
				else:
					exon_dict[tid].append(exon)
		fho.close()

	if args.TE_exon:
		cmd = "cut -f 1-7 "+os.path.abspath(args.TE_exon)+" >"+TE_exon_bed
		subprocess.check_call(cmd, shell=True, executable='/bin/bash')
	else:
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
		os.remove(transcript_exon_bed)

def exon_compare(exon_overlap, cutoff):
	keep_middle = (exon_overlap[6]=="middle") & (abs(exon_overlap[1]-exon_overlap[8])<=cutoff) & (abs(exon_overlap[2]-exon_overlap[9])<=cutoff)
	keep_start_rev = (exon_overlap[6]=="start") & (abs(exon_overlap[1]-exon_overlap[8])<=cutoff) & (exon_overlap[5]=="-")
	keep_start_for = (exon_overlap[6]=="start") & (abs(exon_overlap[2]-exon_overlap[9])<=cutoff) & (exon_overlap[5]=="+")
	keep_end_rev = (exon_overlap[6]=="end") & (abs(exon_overlap[2]-exon_overlap[9])<=cutoff) & (exon_overlap[5]=="-")
	keep_end_for = (exon_overlap[6]=="end") & (abs(exon_overlap[1]-exon_overlap[8])<=cutoff) & (exon_overlap[5]=="+")
	keep_SE = exon_overlap[6]=="SE"
	keep = keep_middle | keep_start_rev | keep_start_for | keep_end_rev | keep_end_for | keep_SE

	return(keep)

def TEA(args):
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
	TE_exon_cluster = pd.merge(TE_exon, exon_cluster, on='exon_ID')
	TE_exon_cluster = TE_exon_cluster[['chr', 'start', 'end', 'transcript_id', 'gene_id', 'strand', 'exon_type', 'exon_class']]
	TE_exon_cluster.to_csv(TE_exon_bed, sep='\t', quoting=csv.QUOTE_NONE, header=0, index=0)

	cmd = 'bedtools intersect -a '+TE_transcript_bed+' -b '+gene_bed+' -v >'+intergenic_bed
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'bedtools intersect -a '+TE_transcript_bed+' -b '+gene_bed+' -s -v >'+intergenic_antisense_bed
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'cat '+intergenic_antisense_bed+' '+intergenic_bed+' | sort | uniq -u >'+antisense_bed
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'bedtools intersect -a '+TE_transcript_bed+' -b '+gene_merged_exon_bed+' -s -v >'+intergenic_antisense_intronic_bed
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'cat '+intergenic_antisense_bed+' '+intergenic_antisense_intronic_bed+' | sort | uniq -u >'+intronic_bed
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = "bedtools intersect -a "+ref_TE_bed+" -b "+TE_exon_bed+" -s -wo  | awk '{if($NF>20){print $0}}' >"+TE_overlap
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'Rscript '+os.path.join(script_dir, 'overlap.R')+' '+args.prefix
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'bedtools intersect -a '+overlap_pos_bed+' -b '+gene_merged_exon_bed+' -s -f 1 -v >'+chimeric_bed
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'Rscript '+os.path.join(script_dir, 'anno.R')+' '+args.prefix
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'rm -rf '+intergenic_bed+' '+intergenic_antisense_bed+' '+antisense_bed+' '+intergenic_antisense_intronic_bed+' '+intronic_bed+' '+TE_overlap+' '+chimeric_bed+' '+overlap_pos_bed
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'rm -rf '+gene_bed+' '+gene_merged_exon_bed+' '+TE_exon_self_overlap+' '+TE_transcript_bed
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

parser = argparse.ArgumentParser(description='TEA: Transposable Element Annotation')
parser.add_argument('-i', '--input', help='Input file of TE transcripts in GTF format (required)')
parser.add_argument('--TE_bed', help='TE position in BED format (required)')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory (default: .)')
parser.add_argument('-p', '--prefix', default='TEA', help='Prefix for output file name (default: TEA)')
parser.add_argument('-a', '--annotation', help='Genome annotation in GTF format (required)')
parser.add_argument('--TE_exon', help='TE exon annotation in BED format')
parser.add_argument('-d', '--exon_diff', type=int, default=5, help='Maximum difference (bp) of exon ends (default: 5)')

args = parser.parse_args()
script_dir = os.path.abspath(os.path.dirname(__file__))
out_dir = os.path.abspath(args.output_dir)

gene_gtf = os.path.abspath(args.annotation)
gene_merged_exon_bed = out_dir+'/gene_merged_exon.bed'
gene_bed = out_dir+'/gene.bed'

ref_TE_bed = os.path.abspath(args.TE_bed)
input_gtf = os.path.abspath(args.input)

TE_exon_self_overlap = args.prefix+'_TE_exon_self_overlap.txt'
TE_exon_ref_overlap = args.prefix+'_TE_exon_ref_overlap.txt'

TE_transcript_bed = out_dir+'/'+args.prefix+'.TE.tran.bed'
transcript_exon_bed = out_dir+'/'+args.prefix+'_transcript_exon.bed'
TE_exon_bed = out_dir+'/'+args.prefix+'_TE_exon_anno.bed'
intergenic_bed = out_dir+'/'+args.prefix+'.TE.intergenic.bed'
intergenic_antisense_bed = out_dir+'/'+args.prefix+'.TE.intergenic.antisense.bed'
antisense_bed = out_dir+'/'+args.prefix+'.TE.antisense.bed'
intergenic_antisense_intronic_bed = out_dir+'/'+args.prefix+'.TE.intergenic.antisense.intronic.bed'
intronic_bed = out_dir+'/'+args.prefix+'.TE.intronic.bed'
chimeric_bed = out_dir+'/'+args.prefix+'.TE.chimeric.bed'
TE_overlap = out_dir+'/'+args.prefix+'.TE.overlap.txt'
overlap_pos_bed = out_dir+'/'+args.prefix+'.TE.overlap.pos.bed'

if not os.path.exists(out_dir):
    os.makedirs(out_dir)


print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running TEA.', flush=True)


with cd(out_dir):
#	deal_ref_anno(args)

	extract_TE_exon(args)

	TEA(args)

	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish.', flush=True)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] TE annotation is done.', flush=True)
