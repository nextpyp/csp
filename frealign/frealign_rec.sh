#!/bin/bash
#
# Control script to submit multiple jobs on a cluster using the Sun
# Grid Engine. Each job processes N particles. N is specified in the
# 'mparameters' file as 'increment'.
#
# $Id: frealign_rec.sh 142 2009-12-16 03:38:21Z fefo $

if ! [ -n "$3" ]; then
    echo "Usage: `basename $0` iteration_number cutoff dataset"
    exit 1
else
	echo "[`date '+%Y%m%d-%H%M%S'`] `basename $0` $@"
fi

iter=$1
cutoff=$2
dataset=$3

iterformat=`printf %02d ${iter}`

boff=`echo ${cutoff}/2 | bc -l`

cd scratch

image_stack=../${dataset}_stack.mrc
last=`header ${image_stack} | grep sections | awk '{print $9}'`
let cutoff=last*cutoff/100
thresh=`cat ${dataset}_CSP_${iterformat}.parx | grep -v C | sort -k 12 | head -$cutoff | tail -1 | awk '{print $12}'`
echo "Will keep $cutoff particles from a total of $last using PR = $thresh"

${CSPDIR}/frealign/mreconstruct.sh ${iter} ${dataset} ${thresh} ${boff}

cd - > /dev/null
