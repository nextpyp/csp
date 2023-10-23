#!/bin/bash
#
# Main function for running frealign. This function calls the function
# for each iteration in the frealign alignament (search) and
# reconstruction.
#
# $Id: frealign_dispatcher.sh 163 2010-01-07 17:08:16Z fefo $

# Test whether command-line argument is present (non-empty).
if ! [ -n "$2" ]; then
    echo -e "NAME\n\t`basename $0`\n"
    echo -e "SYNOPSIS\n\t`basename $0` [-m] max_iter [-c] cutoff_percentage [-i ini_iter] [-p pause_time]"
	echo -e "\nDESCRIPTION\n\tIterate until the max_iter iteration is reach.\n"
	echo -e "MANDATORY PARAMETERS"
	echo -e "\t-m \033[4mmax_iter\033[0m\n\t\tMaximum number of iterations."
	echo -e "\t-c \033[4mcutoff_percentage\033[0m\n\t\tPercentage of particle to preserve."
	echo -e "\nOPTIONAL PARAMETERS"
	echo -e "\t-i \033[4mini_iter\033[0m\n\t\tStart from iteration ini_iter. Default: 2"
	echo -e "\t-p \033[4mXt\033[0m\n\t\tPause time, X is a number and t={s|m|h} (seconds, minutes, hours). Default: 5m"
    exit 1
else
	echo "[`date '+%Y%m%d-%H%M%S'`] `basename $0` $@"
fi

# Initial values for optional variables
pause_time='5m'
ini_iter=2

# Process command line arguments
index=1
while [ $# -gt 0 ]; do
    case "$1" in
        -i) ini_iter=$2; shift;;
        -m) max_iter=$2; shift;;
		-c) cutoff_percentage=$2; shift;;
        -p) pause_time=$2; shift;;
        *)  if [ "$index" == "1" ]; then max_iter=$1; fi
            if [ "$index" == "2" ]; then cutoff_percentage=$1; fi
            if [ "$index" == "3" ]; then pause_time=$1; fi
            index=$((index+1));;
    esac
    shift
done

# Create log directory if needed and name the logfile
mkdir -p log
dts=`date '+%Y%m%d_%H%M%S'`
logfile=log/${dts}_frealign_dispatcher.log

# Set variables values
last_iteration_completed=$((ini_iter-1))
if [ "$max_iter" -lt "$ini_iter" ]; then
	max_iter=$ini_iter
	echo "[WARN] max_iter was set to $ini_iter" > $logfile  
fi

# Check if all the parameter files are present. Note: the parameters
# files are suppose to be name as frealign_parameters_XX
fparamsprefix=frealign_parameters
OK_TO_CONTINUE=1
for (( iter=$ini_iter; iter <= $max_iter; iter++ )); do
	fparams=${fparamsprefix}_`printf %02d $iter`
	if ! [ -e $fparams ]; then
		echo "[ERROR] `basename $0`: can not find $fparams"
		OK_TO_CONTINUE=0
	fi
done
if ! [ $OK_TO_CONTINUE ]; then exit 1; fi

# Set the number of the last completed iteration to 1
echo $last_iteration_completed > .last_iteration_completed

# Run from iter=2 to max_iter.
# Waits until previous iteration is completed to run the next
# one. This is controlled by a number in the file
# .last_iteration_completed corresponding to the last iteration
# completed. The update of .last_iteration_completed is done as the
# last step in the script frealign_iter_implement.sh
for (( iter=$ini_iter; iter <= $max_iter; iter++ )); do
	# Run
	comando="${SPA_DIR}/frealign/frealign_iter.sh $iter $cutoff_percentage"
	echo "[`date '+%Y%m%d_%H%M%S'`] $comando" > $logfile 
	$comando >& $logfile 
	# Wait
	while [ "$last_iteration_completed" -ne "$iter" ]; do
		sleep $pause_time
		last_iteration_completed=`cat .last_iteration_completed`
	done
done