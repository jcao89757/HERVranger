#! usr/bin/env python

# 061617 ZeZhang
# This file is to (1) filter unassigned reads in featureCount results. (2) extract unmapped reads from unmatted reads
# (3) Combine the two part into new fasta files (4) Align them to STAR_061617 (5) FeatureCounts again?
# Try not to use linux cmds. Try portal pycharm (same as matlab)separate sam files into R1 and R2 (check examples from
# original R1 and R2). SAM files may contain same seqID and same seqs, but aligned to different locations.
# Just keep the pair. Scripts start from sys.argv[1]
# Move every possible arguments to sys.argv[n]
# Packed script is considered to treat with only on pair.
# Maybe sam files contain fragments, ie, split R1 and mate to different features.
# Seems like all multiple records in test.fcindex come from Unassigned_MultiMapping.
import re
import os
import sys
import shutil
from Realignment_funct import seqExtract
import time
global pattern
#AlignRef='/home2/s421955/data/genome/hg38mm10/STAR_HERV_061517'
#Data_path='/home2/s421955/projects/retrovirus/data'
#Fasta_path='/archive/BICF/shared/Kidney/rna/RAW'
#name='9267090_1_R1.fastq.gz',for test
start_time=time.time()
AlignRef= sys.argv[1]
Data_path=sys.argv[2]
name_R1=sys.argv[3]
name_R2=sys.argv[4]
Fasta_path=sys.argv[5]
SAMfile=sys.argv[6]
feature_details=sys.argv[7]
output_prefix=sys.argv[8]
#To be added: defalt path is the current path
#To be added: names with paths should get treated

## Built in paths
#SAMfile = os.path.join(Data_path, name_R1, (name_R1+ 'Aligned.out.sam'))
#feature_details = SAMfile + '.featureCounts'
unmated1 = os.path.join(Data_path, output_prefix, (output_prefix + 'Unmapped.out.mate1'))
unmated2 = os.path.join(Data_path, output_prefix, (output_prefix + 'Unmapped.out.mate2'))
## Above paths are basicly unchanged in different users' cases

seqExtract(Fasta_path, Data_path, name_R1,name_R2, feature_details, unmated1, unmated2, output_prefix)

print('Running time: ')
print(time.time()-start_time, 'seconds')
