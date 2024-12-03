#!/usr/bin/env python3
# Author: Janky

import argparse
import re
import random
import os
import csv
import shutil
import subprocess
import pandas as pd
import numpy as np
from datetime import datetime
from collections import defaultdict
import contextlib

@contextlib.contextmanager
def cd(cd_path):
	saved_path = os.getcwd()
	os.chdir(cd_path)
	yield
	os.chdir(saved_path)

def tmpdir(n=16):
	code = ''
	for i in range(n):
		rand_num =random.randint(0,9)
		rand_lower = chr(random.randint(97,122))
		rand_upper = chr(random.randint(65,90))
		s = str(random.choice([rand_num, rand_lower, rand_upper]))
		code += s
	dir_name = 'tmp.'+code
	return dir_name

def extract_exon(gtf_file, exon_dict):
	with open(gtf_file,'r') as fh:
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
	return exon_dict

def extract_ME(exon_dict, ME_file):
	fho = open(ME_file,'w')
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
	fho.close()

def extract_SE(exon_dict, SE_file):
	fho = open(SE_file,'w')
	for tid in exon_dict:
		chrn,transcript_id,gene_id,strand = tid.split(":")
		exon_list = exon_dict[tid]
		if len(exon_list) == 1:
			print("%s\t%d\t%d\t%s\t%s\t%s" % (chrn, int(exon_list[0][0])-1, int(exon_list[0][1]), transcript_id, gene_id, strand), file=fho, flush=True)
	fho.close()

def overlap(args):
	sort_cmd = 'bedtools sort -i '+SE_file+' >'+SE_sort_file
	subprocess.check_call(sort_cmd, shell=True, executable='/bin/bash')

	if args.stranded:
		merge_cmd = 'bedtools merge -i '+SE_sort_file+' -s -c 6 -o distinct >'+SE_merge_file
	else:
		merge_cmd = 'bedtools merge -i '+SE_sort_file+' -c 6 -o distinct | cut -d "," -f1 >'+SE_merge_file
	subprocess.check_call(merge_cmd, shell=True, executable='/bin/bash')

	fho = open(SE_file, 'w')
	with open(SE_merge_file,'r') as fh:
		for line in fh.readlines():
			line = line.strip().split('\t')
			chrn = line[0]
			start = int(line[1])
			end = int(line[2])
			strand = line[3]
			transcript_id = chrn+':'+str(start)+'-'+str(end)+'('+strand+')'
			print("%s\t%d\t%d\t%s\t%s\t%s\t%s" % (chrn, start, end, transcript_id, transcript_id, strand, 'SE'), file=fho, flush=True)
	fho.close()

	if args.stranded:
		overlap_cmd = 'bedtools intersect -a '+SE_file+' -b '+ME_file+' -s -wa -wb >'+overlap_file
	else:
		overlap_cmd = 'bedtools intersect -a '+SE_file+' -b '+ME_file+' -wa -wb >'+overlap_file
	subprocess.check_call(overlap_cmd, shell=True, executable='/bin/bash')

def exon_merge(overlap):
	start = overlap.groupby('ID')['Mstart'].apply(min).reset_index()
	start.columns = ['ID','min']
	end = overlap.groupby('ID')['Mend'].apply(max).reset_index()
	end.columns = ['ID','max']
	start_end = pd.merge(start,end,on='ID')
	data = pd.merge(overlap,start_end,on='ID')
	data['start'] = data['Sstart']
	data['start_type'] = 'SE'
	data['end'] = data['Send']
	data['end_type'] = 'SE'
	start = data[data['Mstart']==data['min']]
	start = start.copy()
	start.loc[start['start']>=start['min'], 'start'] = start.loc[start['start']>=start['min'], 'min']
	start.loc[start['start']>=start['min'], 'start_type'] = start.loc[start['start']>=start['min'], 'Mtype']
	start_stat = start[['ID','Schr','tid','start','start_type']].drop_duplicates()
	end = data[data['Mend']==data['max']]
	end = end.copy()
	end.loc[end['end']<=end['max'], 'end'] = end.loc[end['end']<=end['max'], 'max']
	end.loc[end['end']<=end['max'], 'end_type'] = end.loc[end['end']<=end['max'], 'Mtype']
	end_stat = end[['ID','tid','end','end_type','strand']].drop_duplicates()
	start_stat = start_stat.rename(columns={'Schr':'chr', 'tid':'start_tid'})
	end_stat = end_stat.rename(columns={'tid':'end_tid'})

	stat = pd.merge(start_stat, end_stat, on='ID')
	stat['type'] = 'SE'
	keep_mid = (stat['start_type']=="middle") & (stat['end_type']=="middle")
	stat.loc[keep_mid,'type'] = 'middle'
	keep_start_for = (stat['start_type']!="middle") & (stat['start_type']!="end") & (stat['strand']=="+") & (stat['end_type']!="SE") & (stat['end_type']!="end")
	keep_start_rev = (stat['start_type']!="SE") & (stat['start_type']!="end") & (stat['strand']=="-") & (stat['end_type']!="middle") & (stat['end_type']!="end")
	keep_start = keep_start_for | keep_start_rev
	stat.loc[keep_start,'type'] = 'start'
	stat.loc[keep_start_for,'start_tid'] = pd.NA
	stat.loc[keep_start_rev,'end_tid'] = pd.NA

	keep_end_for = (stat['start_type']!="SE") & (stat['start_type']!="start") & (stat['strand']=="+") & (stat['end_type']!="middle") & (stat['end_type']!="start")
	keep_end_rev = (stat['start_type']!="middle") & (stat['start_type']!="start") & (stat['strand']=="-") & (stat['end_type']!="SE") & (stat['end_type']!="start")
	keep_end = keep_end_for | keep_end_rev
	stat.loc[keep_end,'type'] = 'end'
	stat.loc[keep_end_for,'end_tid'] = pd.NA
	stat.loc[keep_end_rev,'start_tid'] = pd.NA
	stat.loc[stat['type']=='SE','start_tid'] = pd.NA
	stat.loc[stat['type']=='SE','end_tid'] = pd.NA

	return(stat)

def pos_deal(start, end, start_pos_list, end_pos_list):
	if pd.isnull(start_pos_list):
		exon_pos = start
	else:
		start_pos_list = start_pos_list.split(',')
		start_pos = start_pos_list[:start_pos_list.index(start)]
		exon_pos = ','.join(start_pos)+','+start

	if pd.isnull(end_pos_list):
		exon_pos = exon_pos+','+end
	else:
		end_pos_list = end_pos_list.split(',')
		end_pos = end_pos_list[end_pos_list.index(end)+1:]
		exon_pos = exon_pos+','+end+','+','.join(end_pos)

	return(exon_pos)

def TEAM(args):
	ME = pd.read_table(ME_file, header=None)
	ME.columns = ['chr', 'start', 'end', 'tid', 'gid', 'strand', 'type']
	ME['ID'] = ME['chr']+':'+ME['start'].astype('str')+'-'+ME['end'].astype('str')+'('+ME['strand']+')'
	ME['pos'] = ME['start'].astype('str')+','+ME['end'].astype('str')
	ME_pos = ME.groupby('tid')['pos'].apply(lambda x: ','.join(x)).reset_index()
	ME_pos = ME_pos.drop_duplicates('pos')
	ME = ME[ME.tid.isin(ME_pos.tid)]

	overlap = pd.read_table(overlap_file, header=None)
	overlap = overlap[[0,1,2,12,6,7,8,9,10,13]].drop_duplicates()
	overlap.columns = ['Schr', 'Sstart', 'Send', 'strand', 'Stype', 'Mchr', 'Mstart', 'Mend', 'tid', 'Mtype']
	overlap['ID'] = overlap['Schr']+':'+overlap['Sstart'].astype('str')+'-'+overlap['Send'].astype('str')
	overlap = overlap[overlap.tid.isin(ME_pos.tid)]

	SE = pd.read_table(SE_file, header=None)
	SE.columns = ['chr', 'start', 'end', 'tid', 'gid', 'strand', 'type']
	SE['ID'] = SE['chr']+':'+SE['start'].astype('str')+'-'+SE['end'].astype('str')
	SE = SE[~SE.ID.isin(overlap.ID)]

	keep_rm = (overlap['Sstart']-overlap['Mstart'] >= -10) & (overlap['Send']-overlap['Mend'] <= 10)
	rm_ID = overlap[keep_rm]['ID'].drop_duplicates()
	filt = overlap[~overlap.ID.isin(rm_ID)]

	stat_forward = exon_merge(filt[filt['strand']=="+"])
	stat_reverse = exon_merge(filt[filt['strand']=="-"])
	stat = pd.concat([stat_forward, stat_reverse])

	stat_uniq = stat[["ID", "chr", "start", "end", "start_tid", "end_tid", "strand", "type"]].drop_duplicates()
	novel_SE = stat_uniq[stat_uniq['type']=='SE']
	del novel_SE['ID']
	del SE['ID']
	novel_SE.columns = ["chr", "start", "end", "tid", "gid", "strand", "type"]
	SE_out = pd.concat([novel_SE, SE])
	SE_out['tid'] = SE_out['gid'] = SE_out['chr']+':'+SE_out['start'].astype('str')+'-'+SE_out['end'].astype('str')
	SE_out = SE_out.drop_duplicates('tid')

	novel_ME = stat_uniq[stat_uniq['type']!='SE'][["chr", "start", "end", "start_tid", "end_tid", "strand"]]
	novel_ME = pd.merge(novel_ME, ME_pos, left_on='start_tid', right_on='tid', how='left')
	novel_ME = pd.merge(novel_ME, ME_pos, left_on='end_tid', right_on='tid', how='left')

	novel_ME['exon_pos'] = novel_ME.apply(lambda x:pos_deal(str(x['start']), str(x['end']), x['pos_x'], x['pos_y']), axis=1)
	novel_ME['index'] = np.arange(1,novel_ME.shape[0]+1)
	novel_ME['gid'] = "TEAM_novel_MG"+novel_ME['index'].astype('str')
	novel_ME['tid'] = "TEAM_novel_MT"+novel_ME['index'].astype('str')
	novel_ME['exon_num'] = novel_ME['exon_pos'].map(lambda x:int(len(x.split(','))/2))
	novel_ME['chr_list'] = novel_ME.apply(lambda x:','.join(np.repeat(x['chr'],x['exon_num'])), axis=1)
	novel_ME['strand_list'] = novel_ME.apply(lambda x:','.join(np.repeat(x['strand'],x['exon_num'])), axis=1)
	novel_ME['tid_list'] = novel_ME.apply(lambda x:','.join(np.repeat(x['tid'],x['exon_num'])), axis=1)
	novel_ME['gid_list'] = novel_ME.apply(lambda x:','.join(np.repeat(x['gid'],x['exon_num'])), axis=1)
	novel_ME['type_list'] = novel_ME['exon_num'].map(lambda x:'start'+','+','.join(np.repeat('middle',x-2))+','+'end')
	novel_ME.loc[novel_ME['strand']=='-','type_list'] = novel_ME.loc[novel_ME['strand']=='-','exon_num'].map(lambda x:'end'+','+','.join(np.repeat('middle',x-2))+','+'start')
	novel_ME.loc[(novel_ME['strand']=='-') & (novel_ME['exon_num']==2),'type_list'] = 'end,start'
	novel_ME.loc[(novel_ME['strand']=='+') & (novel_ME['exon_num']==2),'type_list'] = 'start,end'
	novel_ME=novel_ME.drop_duplicates('exon_pos')

	pos_list = ','.join(novel_ME['exon_pos']).split(',')
	ME_data = pd.DataFrame(np.reshape(pos_list, (int(len(pos_list)/2),2)), columns=['start', 'end'])
	ME_data['start'] = pd.to_numeric(ME_data['start'])
	ME_data['end'] = pd.to_numeric(ME_data['end'])
	ME_data['chr'] = ','.join(novel_ME['chr_list']).split(',')
	ME_data['strand'] = ','.join(novel_ME['strand_list']).split(',')
	ME_data['tid'] = ','.join(novel_ME['tid_list']).split(',')
	ME_data['gid'] = ','.join(novel_ME['gid_list']).split(',')
	ME_data['type'] = ','.join(novel_ME['type_list']).split(',')
	ME_data = ME_data[["chr", "start", "end", "tid", "gid", "strand", "type"]]
	ME = ME[["chr", "start", "end", "tid", "gid", "strand", "type"]]
	ME_out = pd.concat([ME_data,ME])
	exon_out = pd.concat([ME_out,SE_out])
	exon_out.loc[exon_out['strand']=='.','strand'] = '+'
	exon_out = exon_out.sort_values(by=['chr', 'start', 'end'])
	exon_out.to_csv(exon_file, sep='\t', header=0, index=0)

	cmd = "bedtools intersect -a "+exon_file+" -b "+ref_TE_bed+" -s -wo | awk '{if($NF>20){print $0}}' | cut -f 1-7 | uniq >"+TE_exon_bed
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')
	cmd = 'bedtools coverage -a '+TE_exon_bed+' -b '+bam_file+' -split >'+cov_file
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')
	TE_exon = pd.read_table(cov_file, header=None)
	TE_tid = TE_exon.loc[TE_exon[10]>=0.80,3].drop_duplicates()
	tid_rm = TE_exon.loc[~TE_exon[3].isin(TE_tid),3].drop_duplicates()
	exon_out = exon_out.loc[~exon_out.tid.isin(tid_rm),]
	TE_exon = TE_exon.loc[TE_exon[10]>=0.80,:6].reset_index(drop=True)
	TE_exon.columns = ['chr', 'start', 'end', 'tid', 'gid', 'strand', 'type']

	transcript_start = exon_out.groupby('tid')['start'].apply(min).reset_index()
	transcript_start.columns = ['tid','start']
	transcript_end = exon_out.groupby('tid')['end'].apply(max).reset_index()
	transcript_end.columns = ['tid','end']
	transcript_start_end = pd.merge(transcript_start, transcript_end, on='tid')
	transcript = pd.merge(exon_out[["chr", "tid", "gid", "strand"]].drop_duplicates(), transcript_start_end, on='tid')
	transcript = transcript[["chr", "start", "end", "tid", "gid", "strand"]]
	transcript = transcript.sort_values(by=['chr', 'start', 'end'])
	transcript_forward = transcript[transcript['strand'] == "+"]
	transcript_reverse = transcript[transcript['strand'] == "-"]
	transcript = pd.concat([transcript_forward, transcript_reverse]).reset_index(drop=True)
	transcript['index'] = transcript.index

	GID = TID = transcript_min = transcript_max = 0
	chrom = 'chr1'
	strand = '+'
	for i in range(transcript.shape[0]):
		if (transcript.loc[i,'chr'] != chrom) or (transcript.loc[i,'start'] > transcript_max) or (transcript.loc[i,'strand'] != strand):
			chrom = transcript.loc[i,'chr']
			transcript_min = transcript.loc[i,'start']
			transcript_max = transcript.loc[i,'end']
			strand = transcript.loc[i,'strand']
			GID += 1
			TID = 1
		else:
			TID += 1
			if transcript.loc[i,'end'] > transcript_max:
				transcript_max = transcript.loc[i,'end']
		transcript.loc[i,'gene_id'] = 'TEAM_G'+str(GID)
		if  transcript.loc[i,'tid'] in TE_tid:
			transcript.loc[i,'transcript_id'] = 'TEAM_G'+str(GID)+'_TET'+str(TID)
		else:
			transcript.loc[i,'transcript_id'] = 'TEAM_G'+str(GID)+'_T'+str(TID)

	TE_exon_out = pd.merge(TE_exon, transcript[['tid', 'transcript_id', 'gene_id', 'index']], on='tid')
	TE_exon_out = TE_exon_out[["chr", "start", "end", "transcript_id", "gene_id", "strand", "type"]]
	TE_exon_out.to_csv(TE_exon_bed, sep='\t', header=0, index=0)

	cmd = "bedtools intersect -a "+TE_exon_bed+" -b "+ref_TE_bed+" -s -wo | awk '{if($NF>20){print $0}}' >"+TE_exon_ref_overlap
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	TE_exon_anno = pd.read_table(TE_exon_ref_overlap, header=None)
	TE_exon_anno['ID'] = TE_exon_anno[0]+':'+TE_exon_anno[1].astype('str')+'-'+TE_exon_anno[2].astype('str')+'('+TE_exon_anno[3]+')'
	TE_ID = TE_exon_anno.groupby('ID')[10].apply(lambda x: ','.join(x)).reset_index()
	TE_name = TE_exon_anno.groupby('ID')[11].apply(lambda x: ','.join(x)).reset_index()
	TE_family = TE_exon_anno.groupby('ID')[13].apply(lambda x: ','.join(x)).reset_index()
	TE_class = TE_exon_anno.groupby('ID')[14].apply(lambda x: ','.join(x)).reset_index()
	TE_anno = pd.merge(TE_ID, TE_name, on='ID')
	TE_anno = pd.merge(TE_anno, TE_family, on='ID')
	TE_anno = pd.merge(TE_anno, TE_class, on='ID')
	TE_anno.columns = ['ID', 'TE_ID', 'TE_name', 'TE_family', 'TE_class']
	TE_exon_anno = pd.merge(TE_exon_anno, TE_anno, on='ID')
	TE_exon_anno = TE_exon_anno[[0,1,2,3,4,5,6,'TE_ID','TE_name','TE_family','TE_class']].drop_duplicates()
	TE_exon_anno = TE_exon_anno.sort_values(by=[0,1,2])
	TE_exon_anno.to_csv(out_bed, sep='\t', header=0, index=0)

	exon_out = pd.merge(exon_out, transcript[['tid', 'transcript_id', 'gene_id', 'index']], on='tid')
	exon_out['entry'] = 'exon'
	transcript['entry'] = 'transcript'
	transcript_exon = pd.concat([transcript[['chr', 'entry', 'start', 'end', 'transcript_id', 'gene_id', 'strand', 'index']], exon_out[['chr', 'entry', 'start', 'end', 'transcript_id', 'gene_id', 'strand', 'index']]])
	transcript_exon = transcript_exon.sort_values(by=['index', 'entry', 'start'], ascending=[True, False, True])
	transcript_exon['attr'] = 'transcript_id "'+transcript_exon['transcript_id'].astype('str')+'"; '+' gene_id "'+transcript_exon['gene_id'].astype('str')+'";'
	transcript_exon['start'] = transcript_exon['start']+1
	transcript_exon['anno'] = 'TEAM'
	transcript_exon['score'] = '.'
	transcript_exon = transcript_exon.loc[transcript_exon.transcript_id.isin(TE_exon_out.transcript_id),]
	transcript_exon = transcript_exon[['chr', 'anno', 'entry', 'start', 'end', 'score', 'strand', 'score', 'attr']]
	transcript_exon.to_csv(out_gtf, sep='\t', quoting=csv.QUOTE_NONE, header=0, index=0)


parser = argparse.ArgumentParser(description='TEAM: Transposable Element Associated Merge')
parser.add_argument('--S1', help="Specify a gtf file from StringTie")
parser.add_argument('--S2', help="Specify a gtf file from SERVE")
parser.add_argument('-s', '--stranded', help="Strand-specific RNA-seq", action="store_true")
parser.add_argument('-r', '--ref_TE', help="TE reference file in BED format")
parser.add_argument('-b', '--bam', help="Specify a bam file from STAR aligner")
parser.add_argument('-m', '--merge', default=1, type=int, help='Merge pattern. 1: local, 2: global (default: 1)', choices=[1,2])
parser.add_argument('-o', '--outdir', default='.', help='Output directory (default: .)')
parser.add_argument('-p', '--prefix', default='TEAM', help="Prefix for output file (default: TEAM)")

args = parser.parse_args()
script_dir = os.path.abspath(os.path.dirname(__file__))
out_dir = os.path.abspath(args.outdir)
STRG_gtf = os.path.abspath(args.S1)
SERVE_gtf = os.path.abspath(args.S2)
ref_TE_bed = os.path.abspath(args.ref_TE)
bam_file = os.path.abspath(args.bam)

os.chdir(out_dir)
tmp_dir = tmpdir(16)
while os.path.exists(tmp_dir):
	tmp_dir = tmpdir(16)
os.makedirs(tmp_dir)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to do TE-exon merge.', flush=True)

with cd(tmp_dir):
	SE_file = args.prefix+'.SE.bed'
	SE_sort_file = args.prefix+'.sorted.SE.bed'
	SE_merge_file = args.prefix+'.merge.SE.bed'
	ME_file = args.prefix+'.ME.bed'
	overlap_file = args.prefix+'.overlap.txt'
	exon_file = args.prefix+'.exon.bed'
	cov_file = args.prefix+'.cov.txt'
	TE_exon_bed = args.prefix+'.TE.exon.bed'
	out_bed = args.prefix+'_TE_exon.bed'
	out_gtf = args.prefix+'.gtf'
	TE_exon_ref_overlap = args.prefix+'_TE_exon_ref_overlap.txt'

	exon_dict = defaultdict()

	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] StringTie GTF Input.', flush=True)
	exon_dict = extract_exon(STRG_gtf, exon_dict)

	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] SERVE GTF Input.', flush=True)
	exon_dict = extract_exon(SERVE_gtf, exon_dict)

	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] TEAM merge.', flush=True)
	extract_ME(exon_dict, ME_file)
	extract_SE(exon_dict, SE_file)

	overlap(args)

	if args.merge == 2:
		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Global merge pattern open.', flush=True)
		cmd = 'cut -f1-6 '+overlap_file+' | sort | uniq >'+SE_file
		subprocess.check_call(cmd, shell=True, executable='/bin/bash')
		cmd = 'cut -f8-13 '+overlap_file+' | sort | uniq >>'+SE_file
		subprocess.check_call(cmd, shell=True, executable='/bin/bash')
		overlap(args)

	TEAM(args)

	shutil.copy(out_bed, out_dir)
	shutil.copy(out_gtf, out_dir)

shutil.rmtree(tmp_dir)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish.', flush=True)
