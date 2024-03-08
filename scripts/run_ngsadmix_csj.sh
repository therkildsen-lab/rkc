#!/bin/bash
## This script is used to run ngsAdmix. See http://www.popgen.dk/software/index.php/NgsAdmix for details. 

INPUT_PATH=$1 # Path to the directory where the beagle formatted genotype likelihood file is stored. An example for the Greenland cod data is: /workdir/cod/greenland-cod/angsd/
BEAGLE=$2 # Name the beagle formatted genotype likelihood file. An example for the Greenland cod data is: bam_list_realigned_mincov_contamination_filtered_mindp368_maxdp928_minind167_minq20_downsampled_unlinked.beagle.gz
MINMAF=$3 # Minimum allele frequency filter
MININD=$4
MINK=$5 # Minimum value of K
MAXK=$6 # Maximum value of K
THREADS=${7-8} # Number of threads to use, default value is 8, but the program can use a lot more if they are made available
NGSADMIX=${8:-'/programs/NGSadmix/NGSadmix'} # Path to NGSAdmix, default value is '/programs/NGSadmix/NGSadmix'

PREFIX=${BEAGLE%%.*}

for ((K = $MINK; K <= $MAXK; K++)); do
	#run ngsAdmix
	echo $K
	$NGSADMIX \
	-likes $INPUT_PATH$BEAGLE \
	-K $K \
	-P $THREADS \
	-o $INPUT_PATH'ngsadmix_'$PREFIX'_k'$K \
	-minMaf $MINMAF \
	-minInd $MININD
done
