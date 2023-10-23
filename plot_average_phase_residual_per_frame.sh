#!/bin/bash
#
# Version: $Id: create_new_stacks.sh 511 2010-09-16 23:14:02Z fefo $
#

# Test whether command-line argument is present (non-empty).
if ! [ -n "$1" ]; then 
    echo "Usage: `basename $0` parfile"
    exit 1
fi  

if ! [ -f ${parfile} ]; then
    echo "[ERROR]: File ${parfile} does not exist."
    exit 1
fi

parfile=$1

plotfile=${parfile}_plt
rm -f $plotfile

micro=$2

gnuplot="plot "
frames=`cat $parfile | grep -v C | awk -v micro=$micro '{ if ($8==micro) print $17}' | sort -n | tail -1`
let frames=frames+1
echo Found $frames frames
for (( i=0; i<$frames; i++ ))
do
    cat $parfile | awk -v film=$i -v micro=$micro '{ if ($8==micro && $17==film) print $12 }' | awk '{sum+=$1} END { print sum/NR}' >> $plotfile
    cat $parfile | awk -v film=$i -v micro=$micro '{ if ($8==micro && $17==film) print $12 }' | sort -n > ${plotfile}_${i}
    if [ ${i} -eq ${frames} ]
    then
        gnuplot="$gnuplot \"${plotfile}_${i}\" w l"
    else
        gnuplot="$gnuplot \"${plotfile}_${i}\" w l,"
    fi
done

file=`echo ${parfile%.*}`

# Plot FSC graph
gnuplot << EOF
set terminal png font "Candara Bold" 14
set output "${file}_PR_frames_${micro}.png"
set xlabel "Frame Number"
set ylabel "Phase Residual"
set title "Average Phase Residual / Frame"
set boxwidth 0.5
set style fill solid
plot "${plotfile}" with boxes title '${parfile}'
EOF

# cat $plotfile
rm -f ${plotfile}*