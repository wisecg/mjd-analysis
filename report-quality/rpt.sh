# !/bin/bash
# A cute little wrapper for "reportQuality"
# allowing for easy parallelizing, etc.
# Clint Wiseman, USC/Majorana
# 7/24/2016

# Set environment
UNAME=$(uname)
if [ "${UNAME}" == "Linux" ]; then
	source ${HOME}/EnvBatch.sh
	# printf "PATH:\n$PATH\nLD_LIBRARY_PATH:\n$LD_LIBRARY_PATH\n"
elif [ "${UNAME}" == "Darwin" ]; then
	source $HOME/dev/EnvSetup.sh
	# printf "PATH:\n$PATH\nDYLD_LIBRARY_PATH:\n$DYLD_LIBRARY_PATH\n"
fi

numArgs=$#
if [ $numArgs = 0 ]; then
	printf "    Usage : ./rpt.sh [option] [arguments]\n\n"
	printf "    -split  [OutputFile] [N.lines] [Prefix]: Split a run list into sublists\n"
	printf "    -range  [OutputFile] [Prefix] [Lo] [Hi] : Save start/stop info for many run lists\n"
	printf "    -scpdsf [OutputFile] [Prefix] [Lo] [Hi] : Download slow controls info using run info on PDSF\n"
	printf "    -scloc  [OutputFile] [Prefix] [Lo] [Hi] : Download slow controls info using a range info file\n"
	printf "    -qsub   [Prefix] [Lo] [Hi] : Submit batch jobs for the channel quality report\n"
	printf "    -rpt    [Prefix] [Lo] [Hi] : Do the channel quality report on this node\n"
	printf "    -merge  [OutputFile] [Prefix] [Lo] [Hi] : Merge output files\n\n"
	printf "    ReportQuality requires MkCookie.\n    Have you run it during this session?\n\n"
fi

# (-split) Split up a run list file
if [ $numArgs = 4 ] && [ "$1" = "-split" ]; then
	splitFile=$2
	lines=$3
	prefix=$4
	cd ./runs
	if [ "${UNAME}" == "Linux" ]; then
		echo "This doesn't work on PDSF, stupid."
	elif [ "${UNAME}" == "Darwin" ]; then
		gsplit -dl ${lines} --additional-suffix=.txt ${splitFile} ${prefix}
	fi
	cd ..
	echo "Split, son!"
fi

# (-range / -L) Save Range Info. Run this on PDSF
if [ $numArgs = 5 ] && [ "$1" = "-range" ]; then
	make clean && make
	rangeFile=$2
	prefix=$3
	lo=$4
	hi=$5
	for ((i=lo; i<=hi; i++)); do
		if (( i < 10 )); then
			tmp[i+1]="$prefix""0$i"
		elif (( i >= 10 )); then
			tmp[i+1]="$prefix""$i"
		fi
	done
	fileList=${tmp[@]}
	echo "./reportQuality ${rangeFile} -L ${fileList}"
fi

# (-scloc / -D) Get SC data, using the range info file.
if [ $numArgs = 5 ] && [ "$1" = "-scloc" ]; then
	make clean && make
	slowFile=$2
	prefix=$3
	lo=$4
	hi=$5
	for ((i=lo; i<=hi; i++)); do
		if (( i < 10 )); then
			eval "./reportQuality ${prefix}0$i -D $slowFile"
		elif (( i >= 10 )); then
			eval "./reportQuality ${prefix}$i -D $slowFile"
		fi
	done
fi

# (-qsub / -S) Submit QSUB jobs
if [ $numArgs = 4 ] && [ "$1" = "-qsub" ]; then
	make clean && make
	prefix=$2
	lo=$3
	hi=$4
	for ((i=lo; i<=hi; i++)); do
		if (( i < 10 )); then
			eval "qsub -l h_vmem=2G -q mndl_prod.q runScan.sh $prefix""0$i"
		elif (( i >= 10 )); then
			eval "qsub -l h_vmem=2G -q mndl_prod.q runScan.sh $prefix""$i"
		fi
	done
fi

# (-rpt / -C) Do channel quality report on this local node
if [ $numArgs = 4 ] && [ "$1" = "-rpt" ]; then
	make clean && make
	prefix=$2
	lo=$3
	hi=$4
	for ((i=lo; i<=hi; i++)); do
		if (( i < 10 )); then
			eval "./reportQuality ${prefix}0$i -C"
		elif (( i >= 10 )); then
			eval "./reportQuality ${prefix}$i -C"
		fi
	done
fi

# (-merge / -M) Merge Output Files
if [ $numArgs = 5 ] && [ "$1" = "-merge" ]; then
	make clean && make
	mergeFile=$2
	prefix=$3
	lo=$4
	hi=$5
	for ((i=lo; i<=hi; i++)); do
		if (( i < 10 )); then
			tmp[i+1]="$prefix""0$i"
		elif (( i >= 10 )); then
			tmp[i+1]="$prefix""$i"
		fi
	done
	fileList=${tmp[@]}
	eval "./reportQuality $mergeFile -M ${fileList}"
fi