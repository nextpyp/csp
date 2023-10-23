#!/bin/bash
#
# $Id: mreconstruct.sh 142 2009-12-16 03:38:21Z fefo $

if ! [ -n "$4" ]; then
    echo "Usage: `basename $0` iteration_number dataset PR_cutoff BOFF"
    exit 1
else
	echo "[`date '+%Y%m%d-%H%M%S'`] `basename $0` $@"
fi  

iter=$1
dataset=$2
thresh=$3
boff=$4

let prev=$iter-1

iterf=`printf %02d $iter`
prevf=`printf %02d $prev`

# Parameters files
fparams=../frealign_parameters_${iterf}
if ! [ -f $fparams ]; then
    echo "[ERROR] `basename $0`: Cannot find frealign parameter file $fparams"
    exit 1
fi

mparams=../../microscope_parameters
if ! [ -f $mparams ]; then
    echo "[ERROR] `basename $0`: Cannot find frealign parameter file $mparams"
    exit 1
fi

image_stack=../${dataset}_stack.mrc
radius=`grep radius $fparams | awk '{print $2}'`
target=`grep thresh_refine $fparams | awk '{print $2}'`
pbc=`grep PBC $fparams | awk '{print $2}'`
dang=`grep DANG $fparams | awk '{print $2}'`
itmax=`grep ITMAX $fparams | awk '{print $2}'`
mode=`grep MODE $fparams | awk '{print $2}'`
fpart=`grep FPART $fparams | awk '{print $2}'`
flip=`grep FLIP $fparams | awk '{print $2}'`
rrec=`grep res_reconstruction $fparams | awk '{print $2}'`
rref=`grep res_refinement $fparams | awk '{print $2}'`
rlowref=`grep res_low_refinement $fparams | awk '{print $2}'`
xstd=`grep XSTD $fparams | awk '{print $2}'`
fmag=`grep FMAG $fparams | awk '{print $2}'`
fdef=`grep FDEF $fparams | awk '{print $2}'`
fastig=`grep FASTIG $fparams | awk '{print $2}'`
rbfact=`grep RBFACT $fparams | awk '{print $2}'`
iewald=`grep IEWALD $fparams | awk '{print $2}'`
sym=`grep SYMM $fparams | awk '{print $2}'`
dfsig=`grep dfsig $fparams | awk '{print $2}'`
iblow=`grep IBLOW $fparams | awk '{print $2}'`
fmatch=`grep FMATCH $fparams | awk '{print $2}'`
w_gh=`grep WGH $mparams | awk '{print $2}'`
kV1=`grep kV1 $mparams | awk '{print $2}'`
cs=`grep CS $mparams | awk '{print $2}'`
realmag=`grep REALMAG $mparams | awk '{print $2}'`
pixel_size=`grep pixel_size $mparams | awk '{print $2}'`
dstep=`echo "$pixel_size*$realmag/10000.0" | bc -l`

cp -p ../maps/${dataset}_CSP_${prevf}.mrc ${dataset}_CSP_${iterf}.mrc

last=`header ${image_stack} | grep sections | awk '{print $9}'`

${CSPDIR}/frealign/frealign_v8.exe << eot > ../log/${dataset}_CSP_${iterf}_mreconstruct.log
M,0,F,F,F,F,${iewald},T,F,${fmatch},0,F,${iblow}		!CFORM,IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,FBEAUT,FCREF,FMATCH,IFSC,FSTAT,IBLOW
${radius},0.,${pixel_size},${w_gh},${xstd},${pbc},${boff},${dang},${itmax},10	!RO,RI,PSIZE,WGH,XSTD,PBC,BOFF,DANG,ITMAX,IPMAX
1,1,1,1,1								!MASK
1,${last}								!IFIRST,ILAST
${sym}									!ASYM symmetry card
1.0, ${dstep}, ${target}, ${thresh}, ${cs}, ${kV1}, 0., 0.		!RELMAG,DSTEP,TARGET,THRESH,CS,AKV,TX,TY
${rrec},${rlowref},${rref},${dfsig},${rbfact}				!RREC,RMAX1,RMAX2,DFSIG,RBFACT
${image_stack}
/dev/null
${dataset}_CSP_${iterf}.parx
${dataset}_CSP_${iterf}_dummy.par
${dataset}_CSP_${iterf}_dummy.shft
0., 0., 0., 0., 0., 0., 0., 0.						!terminator with RELMAG=0.0
${dataset}_CSP_${iterf}.mrc
${dataset}_CSP_${iterf}_weights_${iterf}
${dataset}_CSP_${iterf}_map1.mrc
${dataset}_CSP_${iterf}_map2.mrc
${dataset}_CSP_${iterf}_phasediffs
${dataset}_CSP_${iterf}_pointspread
eot

# Clean up this mess
rm -f ${dataset}_CSP_${iterf}_weights_${iterf}
rm -f ${dataset}_CSP_${iterf}_dummy.par
rm -f ${dataset}_CSP_${iterf}_dummy.shft
rm -f ${dataset}_CSP_${iterf}_phasediffs
rm -f ${dataset}_CSP_${iterf}_pointspread

echo 'mreconstruct.com finished' >> ../log/${dataset}_CSP_${iterf}_mreconstruct.log
date

# Compute FSC
proc3d ${dataset}_CSP_${iterf}_map1.mrc ${dataset}_CSP_${iterf}_map2.mrc apix=${pixel_size} fsc=${dataset}_CSP_${iterf}_fsc.txt

# Plot graph
gnuplot << EOF
set terminal png
set output "${dataset}_CSP_${iterf}_fsc.png"
set xlabel "Resolution"
set title "$FSC curve for ${dataset}_CSP, Iteration ${iterf}"
set xr [0:.15]
plot "${dataset}_CSP_${iterf}_fsc.txt" using 1:2 with linespoints
EOF

mv ${dataset}_CSP_${iterf}.mrc ${dataset}_CSP_${iterf}_fsc.png ${dataset}_CSP_${iterf}_fsc.txt ../maps
rm ${dataset}_CSP_${iterf}_map1.mrc ${dataset}_CSP_${iterf}_map2.mrc
