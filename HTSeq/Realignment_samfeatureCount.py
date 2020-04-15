#!/usr/bin/env python
#This file is to count sam features manually, NH:n = n unique alignments
import sys
import re
import pandas as pd
import os
from collections import Counter
samfile=sys.argv[1]
normfile=sys.argv[2]
mode=sys.argv[3]

sam=open(samfile,'r')
ID=[]
feature=[]
NH=[]
strand=[]
status=[]
tmp='TBD'
for lines in sam:
    if lines.find('@')==0:
        continue
    else:
        line=re.split('[\t|\n| ]',lines)
        if line[8]=='-'+tmp:
            status.append(-1)
        else:
            status.append(1)
        ID.append(line[0])
        feature.append(line[2])
        NH_num=re.split('[:]',line[11])
        NH.append(NH_num[2])
        strand.append(line[8])
        tmp=line[8]
info={'ID':ID,'feature':feature,'NH':NH,'strand':strand,'status':status}
info=pd.DataFrame(info)
if mode=='paired':
    if len(info[info.status==1])!=len(info[info.status==-1]):
    	print('Warning: unmatched paired-end sequences.')
    info_pos=info[info.status==1].sort_values(by='feature')
else:
    info_pos=info.sort_values(by='feature')
#-M
feature_freqs=Counter(info_pos['feature'])
d=pd.DataFrame.from_dict(feature_freqs,orient='index')
d.rename(columns={0:'frequencies'})
if mode=='single':
    d['frequencies']=d['frequencies']/2
#-M + fraction
#features=info_pos['feature'].unique()
#info_pos['NH']=[1/int(c) for c in info_pos['NH'].tolist()]
#feature_freqs=[sum(info_pos[info_pos.feature==c].NH) for c in features]
#d=pd.DataFrame(data={'frequencies':feature_freqs},index=features)
d.to_csv(samfile+'.featurefreqs.csv')
if not os.path.isfile(samfile+'.ctrlfreqs.csv') or os.stat(samfile+'.ctrlfreqs.csv').st_size<100:
	print('Norm rewriting')
	norm=pd.read_table(normfile,index_col=0,header=2,sep='\t',comment='#',usecols=[0,6])
	norm.columns=['freq']
	norm=norm[norm.freq!=0]
	feature_norm=norm.index.unique()
	fnorm_freqs=[sum(norm[norm.index==c].freq) for c in feature_norm]
	df=pd.DataFrame(data={'frequencies':fnorm_freqs},index=feature_norm)
	df.index.name=None
	df.to_csv(samfile+'.ctrlfreqs.csv')
	sam.close()
