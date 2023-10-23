#!/bin/bash
#
# This function is called from the C++ program to run the Frealign
# refinement.
#
# Version: $Id: run_frealign_refine.sh 493 2010-08-21 12:53:09Z fefo $
#

if ! [ -n "$1" ]; then
    echo "Usage: `basename $0` parx_file"
    exit 1
# else
#	echo `basename $0` $1
fi  

inDotPar=$1

# Get the basename for this experiment
datasetName=`echo $inDotPar | sed -e "s/\(.*\)_CSP.*par.*/\1/"`
fileBaseName=`echo $inDotPar | sed -e "s/\(.*\)_I.*par.*/\1/" | sed 's|'_CSP'|''|g'`
fileBaseNameFull=`echo $inDotPar | sed -e "s/\(.*\).parx/\1/"`
iter=`cat /scratch/.csp_current_iteration`
let prev=iter-1

# echo $HOSTNAME
# ls -la /scratch

# sed -i -e "s/MASK.*/MASK                    0,0,0,1,1/" ${frealign_parameter_file}

workdir=`pwd`
cd /scratch

raw_images=/scratch/${fileBaseName}_stack.mrc
if ! [ -f ${raw_images} ]
then
    cp ${workdir}/frealign/data/${fileBaseName}_stack.mrc ${raw_images}
fi
reference=/scratch/${datasetName}_CSP_`printf %02d ${prev}`.mrc
#reference=../maps/${datasetName}_CSP_01.mrc

# save current par file
# cp $inDotPar ${fileBaseNameFull}_`printf %02d ${iter}`.parx

# append resolution table to parameter file
resfile=${workdir}/frealign/maps/${datasetName}_`printf %02d ${prev}`.rese
if [ -f $resfile ]
then
    cat ${resfile} >> ${inDotPar}
fi

# cp ${inDotPar} ${inDotPar}_before
# Run FREALIGN
# com="${PYP_DIR}/CSP/frealign/msearch_n.sh ${fileBaseNameFull} ${raw_images} $iter ${reference} ${inDotPar}"
com="${PYP_DIR}/CSP/frealign/msearch_n.py ${fileBaseNameFull} ${raw_images} $iter ${reference} ${inDotPar}"
# echo $com; 
$com
# cp ${inDotPar} ${inDotPar}_after
# cp ${inDotPar} ${inDotPar}_1
# $com
# cp ${inDotPar} ${inDotPar}_2
# $com
# cp ${inDotPar} ${inDotPar}_3
# $com
# cp ${inDotPar} ${inDotPar}_4
# $com
# cp ${inDotPar} ${inDotPar}_5

cd - > /dev/null

# exit 1

# ############################################################################

# # Save the matches file for each iteration?
# saveMatchesMRC=`grep saveMatchesMRC $paramsfile | awk '{print $2}'`

# # Build the name of the output dot par file to be readed
# outDotPar=`echo $inDotPar | sed -e 's/_${prev}_\(.*\)parx/_${iter}_\1parx/'`

# # Grep the msearch log in order to recover the cross-correlation
# # coefficient for each particle projection; save a line for PP in a
# # log file to be readed by the CSP framework.
# fileSuffix=000001_`printf %06d ${nParticles}`
# outCC=`echo ${outDotPar} | sed -e 's/.parx/_CC.log/'`
# grep --text 'CC for particle' ${fileBaseName}_2_msearch_n.log_${fileSuffix} | awk '{print $4"\t"$6"\t"$10}' > ../log/${outCC}

# # Save the msearch_n.log file for this iteration.
# outLog=`echo ${outDotPar} | sed -e 's/.parx/_msearch_n.log/'`
# cp ${fileBaseName}_2_msearch_n.log_${fileSuffix} ../log/${outLog}

# # Rename the output dot par file. Save it for the records.
# if ! [ -f "scratch/${fileBaseName}_2.par_${fileSuffix}" ]; then
    # echo "[ERROR ??] scratch/${fileBaseName}_2.par_${fileSuffix} does not exists."
    # exit 1
# else
    # cp scratch/${fileBaseName}_2.par_${fileSuffix} ${experimentName}/${outDotPar}
# fi

# # Save for the records the 
# if [[ ${saveMatchesMRC} -eq 1 ]]; then
    # outMatchDotMRC=`echo $outDotPar | sed -e 's/.parx/_match.mrc/'`
    # cp scratch/${fileBaseName}_2_match.mrc_${fileSuffix} ${experimentName}/${outMatchDotMRC}
# fi

# # Return to the folder
# cd - > /dev/null
