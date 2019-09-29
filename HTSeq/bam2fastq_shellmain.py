#!/usr/bin/env python
import os
import sys
import pandas as pd
import numpy as np
index_file=sys.argv[1]
#index_file='/home2/s421955/cleaned_data/samples_Sato_v2.csv'
#index_file='/home2/s421955/cleaned_data/samples_TCGA.csv'
data_path=sys.argv[2]
#data_path='/home2/s421955/projects/retrovirus/data'
#output data path
index=pd.read_table(index_file,header=0,sep=',')
j=0
for i in range(0,len(index)):
	line=index.iloc[i,:]
	if line['Data_type'].find('RNA')<0:
		continue
	if pd.isnull(line['Seq_ID']):
		print('Null Seq_ID:',line['Root_tree'])
		continue
	prefix='_'.join([line['Root_tree'],line['Seq_ID']])
	newfilename='.'.join([line['Root_tree'],'sh'])
	out_folder=os.path.join(data_path,prefix)
	input_bam=os.path.join(line['Path'],line['Seq_ID'])
	if not os.path.exists(input_bam):
		print('Bam not exist:',line['Seq_ID'])
		continue
		#That path folder doesn't contain compatitable bamfile
	if not os.path.exists(out_folder):
		os.mkdir(out_folder)
	if len(prefix)>0:
		cmd=' '.join(['perl','bam2fastq.pl',input_bam,out_folder,'40'])
		newR1=os.path.join(out_folder,'.'.join([(line['Seq_ID']+'_R1'),'fastq']))
		newR2=os.path.join(out_folder,'.'.join([(line['Seq_ID']+'_R2'),'fastq']))
		cmd_mv1=' '.join(['mv',os.path.join(out_folder,'fastq1.fastq'),newR1])
		cmd_mv2=' '.join(['mv',os.path.join(out_folder,'fastq2.fastq'),newR2])
		cmd_gzip1=' '.join(['gzip',newR1])
		cmd_gzip2=' '.join(['gzip',newR2])
		cmd_Restar='python Realignment_main.py /home2/s421955/data/genome/hg38mm10/STAR ' \
            '/home2/s421955/data/genome/hg38mm10/STAR_HERV_092717 /home2/s421955/data/genome/hg38mm10/hg38mm10.gtf /home2/s421955/projects/retrovirus/data '+(newR1+'.gz')+' '+(newR2+'.gz')+' '+out_folder+' '+prefix+' test'
	else:
		print('Null prefix:',line['Seq_ID'])
	if j==0:
		newfile=open(newfilename,'w')
		newfile.write('#!/bin/bash \n #SBATCH --job-name=serial \n #SBATCH --partition=super \n #SBATCH --nodes=1 \n #SBATCH --time=10-00:00:00 \n #SBATCH --output=./sbatch_output_%j\n')
	if j<3:
		newfile.write(cmd)
		newfile.write('\n')
		newfile.write(cmd_mv1)
		newfile.write('\n')
		newfile.write(cmd_mv2)
		newfile.write('\n')
		newfile.write(cmd_gzip1)
		newfile.write('\n')
		newfile.write(cmd_gzip2)
		newfile.write('\n')
		newfile.write(cmd_Restar)
		newfile.write('\n')
		j=j+1
	if j==3 or i==len(index)-1:
		newfile.close()
		j=0
