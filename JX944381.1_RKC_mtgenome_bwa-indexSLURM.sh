#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --job-name=bwa_index_JX944381.1_RKC_mtgenome
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=laura.timm@noaa.gov
#SBATCH --output=/home/ltimm/rkc_mito/job_outfiles/bwa-index_JX944381.1_RKC_mtgenome.out

module unload aligners/bwa/0.7.17
module load aligners/bwa/0.7.17

bwa index -p /home/ltimm/rkc_mito/bwa/JX944381.1_RKC_mtgenome /home/ltimm/ref_genomes/JX944381.1_RKC_mtgenome.fasta
