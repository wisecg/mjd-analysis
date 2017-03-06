#!/bin/bash
#$ -cwd
#$ -j y
#$ -o /global/homes/w/wisecg/thresholds/logs/
#$ -P majorana
source /global/homes/w/wisecg/env/EnvBatch.sh
cd /global/homes/w/wisecg/thresholds

echo "Job Start:"
date
echo "Node:  "$HOSTNAME
echo "Job ID:  "$JOB_ID

echo "./thresholds -d $1 $2"
./thresholds -d $1 $2

echo "Job Complete:"
date
