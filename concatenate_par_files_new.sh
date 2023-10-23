#!/bin/bash

function usage {
    echo "Usage: `basename $0` file1 file2 output_file"
    echo -e "\t file1 and file2 (in any order)." 
    echo -e "\t Overwrites the output_file" 
}

if ! [ -n "$3" ]; then usage; exit; fi

output=$3
rm -rf $output

l1=`cat $1 | grep -v C | head -1 | awk '{ print length($0) }'`
l2=`cat $2 | grep -v C | head -1 | awk '{ print length($0) }'`

if [ "$l1" -eq "$l2" ]; then
    echo "ERROR: Files must have different number of columns."
    exit
else
    if [ "$l1" -gt "$l2" ]; then
        standard=$2
        extended=$1
    else
        standard=$1
        extended=$2
    fi
fi

nlinesA=`cat $standard | grep -v C | wc -l`
nlinesB=`cat $extended | grep -v C | wc -l`

if [ "$nlinesA" -ne "$nlinesB" ]; then
    echo "ERROR: Files have different number of lines."
    exit
fi

cat $extended | grep C > $output
cat $extended | grep -v C | cut -b104- > .extended
cat $standard | grep -v C > .standard
paste .standard .extended >> $output
rm -f .standard .extended