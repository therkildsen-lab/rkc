#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --job-name=multiQC
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=laura.timm@noaa.gov
#SBATCH --output=/home/ltimm/rkc_mito/job_outfiles/PCAM-mito-raw_multiQC.out

source /home/ltimm/bin/hydraQC/bin/activate
multiqc /home/ltimm/rkc_mito/fastqc/raw/
