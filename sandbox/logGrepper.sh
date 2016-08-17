#!/bin/bash

#cd /global/project/projectdirs/majorana/data/production/mjdProcessingLogs
#result=$(grep -rl "Warning! Found no LED events" .)
#result2=$(grep -oP 'run \s*\K\d+' $result)
#echo "$result"

echo "Start!"
while read p; do
  #echo $p
  grep -oP 'run \s*\K\d+' $p
done <LEDErrorLogs.txt
