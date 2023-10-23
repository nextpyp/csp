#!/bin/bash
#
# Generate one parameter file per iteration. The parameters are read
# from the file <fparams_all> with this format per line:
# <parameter_name> <value_for_iter_02> <value_for_iter_03> ... <value_for_last_iter>
# For example:
# res_refinement   10.0 10.0 9.0 8.0 7.0 6.0 5.0 4.0
#
# Warning: the fiel <fparams_all> must end with an empty line.
#
# Input: 
# - fparams_all: file with global parameter for all the
# iterations. The name of the output files will be based in the name
# of this file without the _all
#
# $Id: generate_parameters_files.sh 164 2010-01-07 17:14:44Z fefo $

if ! [ -n "$1" ]; then
    echo -e "Usage: `basename $0` file_with_all_parameters"
    exit 1
fi  
fparams_all=$1
fparams=`echo $fparams_all | sed -e 's/_all//'`
fparams_base=${fparams}_base

# Read the configuration for all the iterations
ind=0
exec 0<"$fparams_all"
while read -r line; do
    # The following allows comments, everything after a '#' is ignored
    line=`echo $line | sed -e 's/#.*$//'`
    tmp_array=( $line )
    # The following if allows empty lines
    if [ -n "$tmp_array" ]; then
	p_names[$ind]=${tmp_array[0]}
	if [ "$ind" -ne 0 ]; then
	    # Check the number of iterations in each variable
	    tmp_max_iter=`echo ${#tmp_array[@]} -1 | bc -l`
	    if [ "$max_iter" -ne "$tmp_max_iter" ]; then
		echo "[ERROR] `basename $0`: the number of iteration for ${p_names[$ind]} is incoherent with previous ones."
		exit 1
	    fi
	else
	    max_iter=`echo ${#tmp_array[@]} -1 | bc -l`
	fi
	for (( iter=0; iter<max_iter; iter++ )); do
	    p_ind=$(( max_iter * ind + iter ))
	    p_values[$p_ind]=${tmp_array[$(( iter + 1 ))]}
	done
	ind=$((ind+1))
    fi
done

# Display some information
echo -n 'Variables to process: '
echo `echo ${p_names[@]} | sed -e 's/ /, /g'`

# Copy the base configuration file and substitute the values of the
# variable for the correspondent iteration
tmpfilename=temporalfile
for (( iter=0; iter<max_iter; iter++ )); do
    fparams_iter=${fparams}_`printf %02d $(( iter + 2 ))`
    rm -f fparams_iter # Delete the previous one if exists or not.
    cp ${fparams_base} ${fparams_iter}
    echo -n "Generating parameters file for iteration $((iter+2)): $fparams_iter"
    for (( i=0; i<${#p_names[@]}; i++ )); do
        #echo "sed -e 's/${p_names[$i]}.*/${p_names[$i]}\t${p_values[(( i + iter ))]}/"
	#echo $i, $iter: $(( max_iter * i + iter )) ${p_values[$(( max_iter * i + iter ))]}
	sed -e "s/${p_names[$i]}.*/${p_names[$i]}\t\t${p_values[$(( max_iter * i + iter ))]}/" ${fparams_iter} > $tmpfilename
	mv $tmpfilename ${fparams_iter}
    done
    echo " [Done]"
done



