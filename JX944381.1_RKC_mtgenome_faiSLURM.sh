#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=fai_JX944381.1_RKC_mtgenome
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=laura.timm@noaa.gov
#SBATCH --output=/home/ltimm/rkc_mito/job_outfiles/fai_JX944381.1_RKC_mtgenome.out

module unload bio/samtools/1.11
module load bio/samtools/1.11

samtools faidx /home/ltimm/ref_genomes/JX944381.1_RKC_mtgenome.fasta
