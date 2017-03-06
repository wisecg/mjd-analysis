#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o /global/u1/w/wisecg/channelSel/reportQuality/output
#$ -P majorana

# This script is called from "submitJobs.sh"

# source /etc/profile.d/modules.sh
source /global/homes/w/wisecg/EnvBatch.sh

# USAGE:
# 0. Make sure you've recently ran a make && make clean on PDSF
# NOTE: don't use a path or an extension for the file
#       reportQuality looks in /runs for .txt files
# 1. qsub RunScan.sh Your_File
# 2. qsub -l debug=1 RunScan.sh Your_File
# 3. qsub RunScan.sh Your_File

echo Running reportQuality with input file:
echo $1

cd /global/u1/w/wisecg/channelSel/reportQuality
echo job start
date
MkCookie mjdb.phy.ornl.gov 443 majorana \$t@ckingPb
./reportQuality $1 -C
echo job end
date
echo "Cletus codes good."