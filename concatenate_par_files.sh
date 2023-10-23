#!/bin/bash

function usage {
    echo "Usage: `basename $0` file1 file2 output_file"
    echo -e "\t file1 and file2 (in any order)." 
    echo -e "\t Overwrites the output_file" 
}

if ! [ -n "$3" ]; then usage; exit; fi

file1=${3}_tmp1
cat $1 | grep -v C > ${file1}
file2=${3}_tmp2
cat $2 | grep -v C > ${file2}
output=$3
rm -rf $output

# keep header of FREALIGN's run .par
cat $2 | grep C > ${output}

l1=`cat ${file1} | grep -v C | head -1 | awk '{ print length($0) }'`
l2=`cat ${file2} | grep -v C | head -1 | awk '{ print length($0) }'`

if [ "$l1" -gt "$l2" ]; then
    fileA=${file2}
    fileB=${file1}
    nA=$l2
    nB=$l1
else
    fileA=${file1}
    fileB=${file2}
    nA=$l1
    nB=$l2
fi

nlinesA=`cat $fileA | wc -l`
nlinesB=`cat $fileB | wc -l`

if [ "$nlinesA" -ne "$nlinesB" ]; then
    echo "ERROR: Files have different number of lines."
    exit
fi

# Concatenate the last part of the 
for (( k = 1; k<=$nlinesA; k++ ))
do
    lineA=`awk '{if($1=='${k}') print $0}' ${fileA}`
    lineB=`awk '{if($1=='${k}') print $0}' ${fileB}`
    echo "${lineA:0:$nA}${lineB:$nA}" >> ${output}
done

rm -f ${file1} ${file2}
