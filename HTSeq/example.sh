#!/bin/bash
#SBATCH --job-name=serial                                            # job name
#SBATCH --partition=256GB                                            # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                                    # number of nodes requested by user
#SBATCH --time=20-00:00:00                                            # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=./sbatch_output_%j                                 # redirect both standard output and erro output to the same file
#SBATCH --error=./sbatch_error_%j
#source ~/.bash_profile

##########  JOBSTART  ###################
#Pair-end
python /home2/s421955/projects/retrovirus/code/HERVranger/Realignment_main.py /project/SCCC/Wang_lab/shared/HERV_Ref/STAR /project/SCCC/Wang_lab/shared/HERV_Ref/STAR_HERV_092717 /project/SCCC/Wang_lab/shared/HERV_Ref/hg38mm10.gtf /home2/s421955/projects/retrovirus/data SRR6468300_1.fastq.gz SRR6468300_2.fastq.gz /project/bioinformatics/Xiao_lab/shared/genomics/dbgap/dbGaP-17356/ RCC_106
#Single-end
python /home2/s421955/projects/retrovirus/code/HERVranger/Realignment_main.py /project/SCCC/Wang_lab/shared/HERV_Ref/STAR /project/SCCC/Wang_lab/shared/HERV_Ref/STAR_HERV_092717 /project/SCCC/Wang_lab/shared/HERV_Ref/hg38mm10.gtf /home2/s421955/projects/retrovirus/data SRR6468300_1.fastq.gz NA /project/bioinformatics/Xiao_lab/shared/genomics/dbgap/dbGaP-17356/ RCC_106
