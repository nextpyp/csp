#!/bin/bash
#
#   Script for ploting phase residuals resulting from FREALIGN refinement.
#
# $Id: frealign_plt.sh 142 2009-12-16 03:38:21Z fefo $

if ! [ -n "$2" ]; then
    echo "Usage: `basename $0` iteration_number dataset nodisplay"
    exit 1
else
	echo "[`date '+%Y%m%d-%H%M%S'`] `basename $0` $@"
fi

iter=$1
dataset=$2
iterf=`printf %02d $iter`

prs=maps/${dataset}_CSP_${iterf}_PR

# extract phase residuals and sort them
cat scratch/${dataset}_CSP_${iterf}.parx | grep -v "C" | sort -k 12 | awk '{print $12}' > ${prs}.txt

# extract phase residuals and sort them
cat scratch/${dataset}_CSP_${iterf}.parx | grep -v "C" | awk '{print $12}' > ${prs}_unsorted.txt

gnuplot << EOF
set terminal png
set output "${prs}.png"
set xlabel "Particle Number"
set ylabel "Phase Residual"
set title " ${dataset}, Iteration ${iterf}"
plot "${prs}_unsorted.txt", "${prs}.txt"
EOF

#mv ${prs}_unsorted.txt ${prs}.txt
rm -f ${prs}_unsorted.txt ${prs}.txt

if ! [ -n "$3" ]; then
    display ${prs}.png
fi