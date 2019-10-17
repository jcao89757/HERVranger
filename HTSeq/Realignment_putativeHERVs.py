#! usr/bin/env python

# 061617 ZeZhang
# This file is to (1) filter unassigned reads in featureCount results. (2) extract unmapped reads from unmatted reads
# (3) Combine the two part into new fasta files (4) Align them to STAR_061617 

import re
import os
import sys
import shutil
from Realignment_funct import seqExtract
import time
global pattern

start_time=time.time()
AlignRef= sys.argv[1]
Data_path=sys.argv[2]
name_R1=sys.argv[3]
name_R2=sys.argv[4]
Fasta_path=sys.argv[5]
SAMfile=sys.argv[6]
feature_details=sys.argv[7]
output_prefix=sys.argv[8]
mode=sys.argv[9]

## Built in paths
unmated1 = os.path.join(Data_path, output_prefix, (output_prefix + 'Unmapped.out.mate1'))
unmated2 = os.path.join(Data_path, output_prefix, (output_prefix + 'Unmapped.out.mate2'))

seqExtract(Fasta_path, Data_path, name_R1,name_R2, feature_details, unmated1, unmated2, output_prefix, mode)

print('Running time: ')
print(time.time()-start_time, 'seconds')
