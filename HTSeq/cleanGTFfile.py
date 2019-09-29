#! usr/bin/env python
# This script is to delete all suspective ERVs in hg38mm10.gtf file
# Separate hg38mm10.gtf into two parts for viewing the deleted lines.
# Use dictionary to match with ERV index
import numpy as np
import os
import re
#ERVind=pd.read_table('/Users/zhangze/Google Drive/codes and scripts to share/refID_ERVinhg38.txt',sep=' ',header=None,names=None,index_col=None)
ERVind=open('/Users/zhangze/Google Drive/codes and scripts to share/refID_ERVinhg38.txt','r')
dict_ERV={}
# re.split for multiple split markers
for line in ERVind.readlines():
    tmp=re.split('[.|\n]',line)
    print tmp
    dict_ERV[tmp[0]]=None
#oldGTF=open('/Users/zhangze/Google Drive/hg38mm10.gtf')
oldGTF=open('/Users/zhangze/Desktop/hg38mm10_deleted.gtf')
#double check since hg38mm10.gtf+hg38mm10_deleted=original hg38mm10.gtf
filename_filtered='/Users/zhangze/Desktop/hg38mm10_filt2.gtf'
filename_deleted='/Users/zhangze/Desktop/hg38mm10_removed.gtf'
newfile_filt=open(filename_filtered,'w')
newfile_delet=open(filename_deleted,'w')

for line in oldGTF:
    tmp_GTF=line.split()
    tmp_GeneID=re.split('["|;|,]',tmp_GTF[9])
    tmp_transID = re.split('["|;|,]', tmp_GTF[11])
    #tmp_GeneID[1]/tmp_transID[1] are the strs expected
    if dict_ERV.has_key(tmp_GeneID[1]) or dict_ERV.has_key(tmp_transID[1]):
        newfile_delet.write(line)
    else:
        newfile_filt.write(line)

newfile_filt.close()
newfile_delet.close()
