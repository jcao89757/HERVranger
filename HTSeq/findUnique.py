
#Use together with awk '{print $1}' SAM19944317_37049_R1.fastq.gzAligned.out.sam.featureCounts 
#> SAM19944317_37049_R1.fastq.gzAligned.out.sam.featureCounts.index
# If print 'find one' then not all same IDs form consequtively.

#! usr/bin/env python 
f=open('/home2/s421955/projects/retrovirus/data/19424793_1_R1.fastq.gz/19424793_1_R1.fastq.gzAligned.out.sam.featureCounts.index','r')
dict_unique={}
initpattern=f.next()
dict_unique[initpattern]=0
i=0
for line in f:
    i=i+1
    if i%10000==0:
        print 'still working'
    tmp=line
    if not tmp.find(initpattern):# tmp is a new feature
        initpattern=tmp
        dict_unique[initpattern]=0
        if dict_unique.has_key(tmp):
            print 'find one'
            break

