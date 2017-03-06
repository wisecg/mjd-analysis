#!/bin/bash
# C. Wiseman, USC.

function runThresh()
{
  # DS map (from DataSetInfo.hh): {0,76},{1,51},{3,24},{4,22},{5,80}
  dsNum=0
  dsMax=76
  echo "running MkCookie"
  # you have to run MkCookie with a password here.  Ask clint for the line to put here.
  echo "make -s"
  make -s
  for ((i=0;i<=${dsMax};i++)); do
    echo $i
    # ./thresh-job.sh ${dsNum} $i
    qsub thresh-job.sh ${dsNum} $i
  done

  # qsub thresh-job.sh 0 54
}

function mergeThresh()
{
  # DS map (from DataSetInfo.hh): {0,76},{1,51},{3,24},{4,22},{5,80}
  dsNum=0
  dsMax=76
  outFile="thresholdsDS${dsNum}.root"
  args=()
  for ((i=0; i<=${dsMax}; i++)); do
      args+=("thresholdsDS${dsNum}_$i.root")
  done
  echo "${args[@]}"
  thisDir=`pwd`
  cd data
  hadd $outFile ${args[@]}
  cd $thisDir
}

function findExposure()
{
  arr=(0 1 3 4 5)
  for i in "${arr[@]}"; do
    echo $i
    make -s && ./thresholds -e final/thresholdsDS$i.root $i
  done
}

# =================================
# runThresh
# mergeThresh
findExposure