#!/bin/bash
#
# Generate empty particle parameter list for FREALIGN refinement
# using the defocus values specified for each micrograph.
# Script takes 2 parameters: dataset + txt file with defocus values
# e.g. frealign_ctf.sh gp140_17b parnum_defocus.txt
#
# $Id: frealign_generate_par_file.sh 124 2009-12-06 02:16:40Z fefo $


if ! [ -n "$6" ]; then
    echo -e "Usage: `basename $0` total_particles lower_index upper_index defocus micrograph_index dataset"
    exit 1
fi  

total_particles=$1
lower_index=$2
upper_index=$3
defocus=$4
arg=`head -$5 ../micrographs | tail -1`
dataset=$6

magnification=`grep REALMAG ../microscope_parameters | awk '{print $2}'`

cd ../scratch/$arg

${SPA_DIR}/frealign_v8.08/bin/frealign_v8.exe << eot
M,-3,F,F,F,F,0,F,F,F,0,F,1	!CFORM,IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,FBEAUT,FCREF,FMATCH,IFSC,FSTAT,IBLOW
125,0.0,3.11,0.07,0,100,0.0,300.0,1,1	!RO,RI,PSIZE,WGH,XSTD,PBC,BOFF,DANG,ITMAX,IPMAX
0,0,0,0,0								!MASK
${lower_index},${upper_index}								!IFIRST,ILAST 
0
1.0, 14.0, 25, 62, 2.0, 120, 0., 0.		!RELMAG,DSTEP,TARGET,THRESH,CS,AKV,TX,TY
15,100,15,150,0				!RREC,RMAX1,RMAX2,DFSIG,RBFACT
../../${dataset}_stack.mrc
dummy_match.mrc
${total_particles},${magnification},1,${defocus},${defocus},45.,0					!NIN,ABSMAGPIN,IFILMIN,DFMID1IN,DFMID2IN,ANGASTIN,MORE
dummy_${arg}.par
dummy.shft
-100., 0., 0., 0., 0., 0., 0., 0.					!terminator with RELMAG=-100.0 to skip 3D reconstruction
../model.mrc
dummy_weights
dummy_map1
dummy_map2
dummy_phasediffs
dummy_pointspread
eot

mv dummy_${arg}.par ..