#!/bin/bash
#
# This function is called from the C++ program to run the Frealign
# reconstruction.
#
# Version: $Id: run_frealign_reconstruction.sh 493 2010-08-21 12:53:09Z fefo $
#

if ! [ -n "$1" ]; then
    echo "Usage: `basename $0` parx_file"
    exit 1
else
	echo "[`date '+%Y%m%d-%H%M%S'`] `basename $0` $@"
fi  

inDotPar=$1

# store copy of extended .par file containing internal state
cp frealign/scratch/$inDotPar frealign/scratch/${inDotPar}t

# Get the basename for this experiment
dataset=`echo $inDotPar | sed -e "s/\(.*\)_CSP.*parx.*/\1/"`

iter=`cat .csp_current_iteration`
iterf=`printf %02d $iter`
cutoff=100

# Change folder
cd frealign

# refine shifts once more
${CSPDIR}/frealign/frealign_ali.sh ${iter} ${dataset}

# Run frealign reconstruction
${CSPDIR}/frealign/frealign_rec.sh ${iter} ${cutoff} ${dataset}

# Generate phase residual plots
${CSPDIR}/frealign/frealign_plt.sh ${iter} ${dataset} 1

# cp scratch/$inDotPar maps/${inDotPar}f
${CSPDIR}/concatenate_par_files.sh scratch/$inDotPar scratch/${inDotPar}t ${dataset}_01.parx maps/$inDotPar
rm -f scratch/${inDotPar}t

# Get back
cd - > /dev/null

# merge all log files into single file
cd frealign/log
rm -f ${dataset}_CSP_${iterf}_msearch_n.tgz
listoffiles=`ls ${dataset}_CSP_T??_?????_????????_${iterf}_msearch_n.log`
tar cvfz ${dataset}_CSP_${iterf}_msearch_n.tgz ${listoffiles} > /dev/null
rm -rf ${listoffiles}
cd - > /dev/null

# get rid of unnecesary .parx files
listoffiles=`ls frealign/scratch/${dataset}_CSP_T??_?????_????????.parx`
# rm -f ${listoffiles}

# Save copy of most current parameters to maps directory
# mv frealign/scratch/${dataset}_CSP_${iterf}.parx frealign/maps
