#!/bin/bash
#
# Script for running frealing refinement called from frealign_ali.sh
#
# $Id: msearch_n.sh 172 2010-06-16 12:37:32Z abarte $

if ! [ -n "$5" ]; then
    echo "Usage: `basename $0` dataset image_stack iteration reference dot_par_in"
    exit 1
#else
#	echo "[`date '+%Y%m%d-%H%M%S'`] `basename $0` $@"
fi  

source ${SPA_DIR}/spa.bashrc

dataset=$1
image_stack=$2
iter=$3
reference=$4
dot_par_in=$5

fparams=/scratch/frealign_parameters_`printf %02d $iter`
if ! [ -f $fparams ]; then
	echo ERROR - Frealign parameter file $fparams not found.
	exit 1
fi

mode=`grep MODE $fparams | awk '{print $2}'`
radius=`grep radius $fparams | awk '{print $2}'`
target=`grep thresh_refine $fparams | awk '{print $2}'`
thresh=`grep thresh_reconst $fparams | awk '{print $2}'`
pbc=`grep PBC $fparams | awk '{print $2}'`
boff=`grep BOFF $fparams | awk '{print $2}'`
dang=`grep DANG $fparams | awk '{print $2}'`
itmax=`grep ITMAX $fparams | awk '{print $2}'`
fpart=`grep FPART $fparams | awk '{print $2}'`
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
# mask=`grep MASK $fparams | awk '{print $2}'`

# over-ride parameters for CSP refinement
mask=0,0,0,1,1 										# force refinement of shifts only
mode=1
sym=0

mparams=/scratch/microscope_parameters
if ! [ -f $mparams ]; then
	echo ERROR - Frealign parameter file $mparams not found.
	exit 1
fi

pixel_size=`grep pixel_size $mparams | awk '{print $2}'`

# EXPERIMENTAL

# retrieve number of current inner iteration
inner_iter=`cat /scratch/.csp_current_inner_iteration`

# retrieve resolution limits for refinement
rlowref=`grep ^LowResolutionForRefinement /scratch/parameters.config | awk '{print $'${inner_iter}'}'`
rref=`grep ^HighResolutionForRefinement /scratch/parameters.config | awk '{print $'${inner_iter}'}'`

# retrieve binning for refinement
binning=`grep ^BinningForRefinement /scratch/parameters.config | awk '{print $'${inner_iter}'}'`
if [ -z $binning ]
then
	binning=1
fi

if [ "$binning" -gt "1" ]
then
	pixel_size=`echo "${pixel_size}*${binning}" | bc -l`
	bin=`printf %02d ${binning}`
	bin_image_stack=`echo ${image_stack} | sed -e "s|_stack|_stack_bin${bin}|"`
	bin_reference=`echo ${reference} | sed -e "s|.mrc|_bin${bin}.mrc|"`
	image_stack=${bin_image_stack}
	reference=${bin_reference}
fi
mw=`grep particle_mw /scratch/.parameters.config | awk '{print $2}'`
w_gh=`grep WGH $mparams | awk '{print $2}'`
kV1=`grep kV1 $mparams | awk '{print $2}'`
cs=`grep CS $mparams | awk '{print $2}'`
realmag=`grep REALMAG $mparams | awk '{print $2}'`
dstep=`echo "$pixel_size*$realmag/10000.0" | bc -l`

particles=`${IMOD_DIR}/bin/header ${image_stack} | grep sections | awk '{print $9}'`

# run FREALIGN

##${SPA_DIR}/frealign_v8.09/bin/frealign_v8_mp_do_not_write_masks.exe << eot >& ../log/${dataset}_`printf %02d ${iter}`_msearch_n.log1
##M,${mode},${fmag},${fdef},${fastig},${fpart},${iewald},T,F,${fmatch},0,F,${iblow}	!CFORM,IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,FBEAUT,FCREF,FMATCH,IFSC,FSTAT,IBLOW
#${SPA_DIR}/frealign_v8.10_intel/bin/frealign_v8_mp_do_not_write_masks.exe << eot >& /dev/null # ${dataset}_`printf %02d ${iter}`_msearch_n.log1
#M,${mode},${fmag},${fdef},${fastig},${fpart},${iewald},T,F,F,${fmatch},0,F,${iblow}	!CFORM,IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,FBEAUT,FFILT,FBFACT,FMATCH,IFSC,FSTAT,IBLOW
#${radius},0.,${pixel_size},${w_gh},${xstd},${pbc},${boff},${dang},${itmax},10	!RO,RI,PSIZE,WGH,XSTD,PBC,BOFF,DANG,ITMAX,IPMAX
#$mask									!MASK
#1,${particles}								!IFIRST,ILAST 
#${sym}
#1.0, ${dstep}, ${target}, ${thresh}, ${cs}, ${kV1}, 0., 0.		!RELMAG,DSTEP,TARGET,THRESH,CS,AKV,TX,TY
#${rrec},${rlowref},${rref},${dfsig},${rbfact}				!RREC,RMAX1,RMAX2,DFSIG,RBFACT
#${image_stack}
#${dataset}_`printf %02d ${iter}`_match.mrc
#${dot_par_in}
#${dot_par_in}
#/dev/null
#-100., 0., 0., 0., 0., 0., 0., 0.					!terminator with RELMAG=-100.0 to skip 3D reconstruction
#${reference}
#${dataset}_`printf %02d ${iter}`_weights_`printf %02d ${iter}`
#${dataset}_`printf %02d ${iter}`_map1.mrc
#${dataset}_`printf %02d ${iter}`_map2.mrc
#${dataset}_`printf %02d ${iter}`_phasediffs
#${dataset}_`printf %02d ${iter}`_pointspread
#eot

${SPY_DIR}/frealign_v9.08/bin/frealign_v9.exe << eot >& /dev/null # ${dataset}_`printf %02d ${iter}`_msearch_n.log1
M,${mode},${fmag},${fdef},${fastig},${fpart},${iewald},T,F,F,${fmatch},0,F,${iblow},1	!CFORM,IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,FBEAUT,FFILT,FBFACT,FMATCH,IFSC,FDUMP,IMEM,INTERP
${radius},0.,${pixel_size},${mw},${w_gh},${xstd},${pbc},${boff},${dang},${itmax},10	!RO,RI,PSIZE,MW,WGH,XSTD,PBC,BOFF,DANG,ITMAX,IPMAX
$mask									!MASK
1,${particles}								!IFIRST,ILAST 
${sym}
1.0, ${dstep}, ${target}, ${thresh}, ${cs}, ${kV1}, 0., 0.		!RELMAG,DSTEP,TARGET,THRESH,CS,AKV,TX,TY
${rrec},${rlowref},${rref},${rref},${dfsig},${rbfact}			!RREC,RMIN,RMAX,RCLAS,DFSTD,RBFACT
${image_stack}
${dataset}_`printf %02d ${iter}`_match.mrc
${dot_par_in}
${dot_par_in}
/dev/null
-100., 0., 0., 0., 0., 0., 0., 0.					!terminator with RELMAG=-100.0 to skip 3D reconstruction
${reference}
${dataset}_`printf %02d ${iter}`_weights_`printf %02d ${iter}`
${dataset}_`printf %02d ${iter}`_map1.mrc
${dataset}_`printf %02d ${iter}`_map2.mrc
${dataset}_`printf %02d ${iter}`_phasediffs
${dataset}_`printf %02d ${iter}`_pointspread
eot

#${SPA_DIR}/frealign_v8.08/bin/frealign_v8.exe << eot >& ../log/${dataset}_`printf %02d ${iter}`_msearch_n.log2
# ${SPA_DIR}/frealign_v8.08/bin/frealign_v8.exe << eot >& /dev/null
# M,${mode},${fmag},${fdef},${fastig},${fpart},${iewald},T,F,${fmatch},0,F,${iblow}	!CFORM,IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,FBEAUT,FCREF,FMATCH,IFSC,FSTAT,IBLOW
# ${radius},0.,${pixel_size},${w_gh},${xstd},${pbc},${boff},${dang},${itmax},10	!RO,RI,PSIZE,WGH,XSTD,PBC,BOFF,DANG,ITMAX,IPMAX
# $mask									!MASK
# 1,${particles}								!IFIRST,ILAST 
# ${sym}
# 1.0, ${dstep}, ${target}, ${thresh}, ${cs}, ${kV1}, 0., 0.		!RELMAG,DSTEP,TARGET,THRESH,CS,AKV,TX,TY
# ${rrec},${rlowref},${rref},${dfsig},${rbfact}				!RREC,RMAX1,RMAX2,DFSIG,RBFACT
# ${image_stack}
# ${dataset}_`printf %02d ${iter}`_match.mrc
# ${dot_par_in}
# ${dot_par_in}
# /dev/null
# -100., 0., 0., 0., 0., 0., 0., 0.					!terminator with RELMAG=-100.0 to skip 3D reconstruction
# ${reference}
# ${dataset}_`printf %02d ${iter}`_weights_`printf %02d ${iter}`
# ${dataset}_`printf %02d ${iter}`_map1.mrc
# ${dataset}_`printf %02d ${iter}`_map2.mrc
# ${dataset}_`printf %02d ${iter}`_phasediffs
# ${dataset}_`printf %02d ${iter}`_pointspread
# eot
#
#${SPA_DIR}/frealign_v8.08/bin/frealign_v8.exe << eot >& ../log/${dataset}_`printf %02d ${iter}`_msearch_n.log3
# ${SPA_DIR}/frealign_v8.08/bin/frealign_v8.exe << eot >& /dev/null
# M,${mode},${fmag},${fdef},${fastig},${fpart},${iewald},T,F,${fmatch},0,F,${iblow}	!CFORM,IFLAG,FMAG,FDEF,FASTIG,FPART,IEWALD,FBEAUT,FCREF,FMATCH,IFSC,FSTAT,IBLOW
# ${radius},0.,${pixel_size},${w_gh},${xstd},${pbc},${boff},${dang},${itmax},10	!RO,RI,PSIZE,WGH,XSTD,PBC,BOFF,DANG,ITMAX,IPMAX
# $mask									!MASK
# 1,${particles}								!IFIRST,ILAST 
# ${sym}
# 1.0, ${dstep}, ${target}, ${thresh}, ${cs}, ${kV1}, 0., 0.		!RELMAG,DSTEP,TARGET,THRESH,CS,AKV,TX,TY
# ${rrec},${rlowref},${rref},${dfsig},${rbfact}				!RREC,RMAX1,RMAX2,DFSIG,RBFACT
# ${image_stack}
# ${dataset}_`printf %02d ${iter}`_match.mrc
# ${dot_par_in}
# ${dot_par_in}
# /dev/null
# -100., 0., 0., 0., 0., 0., 0., 0.					!terminator with RELMAG=-100.0 to skip 3D reconstruction
# ${reference}
# ${dataset}_`printf %02d ${iter}`_weights_`printf %02d ${iter}`
# ${dataset}_`printf %02d ${iter}`_map1.mrc
# ${dataset}_`printf %02d ${iter}`_map2.mrc
# ${dataset}_`printf %02d ${iter}`_phasediffs
# ${dataset}_`printf %02d ${iter}`_pointspread
# eot
#