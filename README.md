# HERVranger
An alignment tool targeting the footprints of HERVs in human genomes
## Introduction
HERVranger is a bioinformatics pipeline for the analysis of the footprints of human endogenous retroviruses (HERVs) from human genomes.
## Getting started
## Dependencies
Linux (x86_64-redhat-linux-gn) shell (4.2.46(2))

Python (>=2.7.16)

STAR (2.5.2b preferred)

featureCounts (subread-1.5.0 only)

Only run the jobs on 256GB nodes.

**Python packages**

os, sys, re, shutil, time, pandas, collections, itertools

## Guided tutorial
### Input data

library1: a STAR reference library build from HERV-removed Hg38 genome reference downloaded from NCBI RefSeq. Please use [/project/SCCC/Wang_lab/shared/HERV_Ref/STAR].

library2: a STAR reference library build from HERV RNA-Seqs. Please use [/project/SCCC/Wang_lab/shared/HERV_Ref/STAR_HERV_092717].

ref.gtf: a featureCounts reference GTF file with HERV removed. Please use [/project/SCCC/Wang_lab/shared/HERV_Ref/hg38mm10.gtf]

output: a path to build a result folder (named by Prefix) to write result files and intermediate files, total usage could be over 10 GB for one pair of RNA-Seqs.

R1.fastq.gz, R2.fastq.gz: a pair of pair-end RNA-Seqs.

To apply HERVranger to single-strand alignments, please input R1.fastq.gz NA, i.e. replace the R2.fastq.gz with a string 'NA' and keep the space between them.

Path_to_saved_RNA_Seqs: where the pair of RNA-Seqs stores.

Prefix: a user decided sample name.

### Example shell cmds
```{r}
python Realignment_maint.py /path/to/library1 /path/to/library2 /path/to/ref.gtf /path/output R1.fastq.gz R2.fastq.gz /path/to/saved/RNA_Seqs Prefix
```
If an index file like the example [samples_example.xlsx] is provided, the [jobupload.py] script may be applied to generate large batches of job scripts. The index file is expected to contain one 'Root_tree' column, one 'Seq_ID' column and one 'Path' column.
## Version update
1.0.0: First release. (09-28-2019)
1.1.0: Update. Single-strand alignment mode was provided. (10-17-2019)
