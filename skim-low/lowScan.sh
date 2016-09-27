#!/bin/csh -f
#$ -cwd
#$ -j y
#$ -o /global/homes/w/wisecg/low-skim/output
#$ -P majorana

# Example usage options
# 1. qsub RunScan.csh Your_File.txt
# 2. qsub -l debug=1 RunScan.csh Your_File.txt
# 3. qsub -l h_vmem=2G -q mndl_prod.q RunScan.csh

source /global/homes/w/wisecg/env/EnvBatch.sh
cd /global/homes/w/wisecg/low-skim/
echo job start
date
./skim_mjd_data 1 $1 ./output/
echo job end
date
