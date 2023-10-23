#!/bin/bash
#
# Version: $Id: compare_all_experiments.sh 566 2011-02-27 22:25:19Z fefo $
#

# function multiecho {
#     str=$1
#     # Calculates number of available columns on the terminal
#     COLUMNS=$( stty -a | head -n 1 | awk '{print $7}' | rev | cut -b 2- | rev )
#     echo $1 columns times $2
#     for i in $(seq 1 $1)
#     do
#         echo -n -e "$2"
#     done
# }

##########################################################################
function Usage {

    echo -e "\nNAME"
    echo -e "\t`basename $0`\n"

    echo -e "SYNOPSIS"
    echo -e "\t`basename $0` EXPERIMENT_BASENAME GT_MRC_FILE\n"

    echo "DESCRIPTION"
    echo -e "\tProcess all the experiments named {EXPERIMENT_BASENAME}_ini_*_noise_*"
    # echo -e "\t\${EXPERIMENT_BASENAME}_ini_\${nIni}_noise_\${nNoise}. Compare the results with"
    # echo -e "\tthe groundtruth values in groundtruth par file located in the same folder."

    # echo -e "\nOPTIONS" 

    # echo -e "\n\t-iter ITER " 
    # echo -e "\t\tProcess only the ITER iteration. (Currently not functional)" 

    # echo -e "\n\t-csp" 
    # echo -e "\t\tProcess only the CSP experiment. (Currently not functional)" 

    # echo -e "\n\t-usp" 
    # echo -e "\t\tProcess only the USP experiment. (Currently not functional)" 

    # echo -e "\n\t-nc\n\t--no-cuts" 
    # echo -e "\t\tDo not write the model cuts." 

    # echo -e "\n\t-yerrorbars" 
    # echo -e "\t\tPlot the PR with error bars." 

    # echo -e "\n\t-nmr\n\t--do-not-use-max-radius"
    # echo -e "\t\tDo not scale the polar graph to the max radius, instead use 1." 

    # echo -e "\n\t-map FILE\n\t--gt-map FILE"
    # echo -e "\t\tGroundtruth density map MRC file." 

}

#if [ "${CLUSTER_ID}" == "NIH" ]; then export DISPLAY=biowulf.nih.gov:51.0; fi

##########################################################################

# With no arguments show the usage help.
if [ "$#" == "0" ]; then Usage; exit; fi

# Default values for some variables.
outputtype=png
write_model_cuts=1
use_max_radius=1
yerrorbars=0
gt_map_mode=0
itermode=0
uspmode=-1
cspmode=-1
nIni=-1
nNoise=-1

# # Process options.
# while [ "$1" != "" ]
# do
#     case $1 in
#         -i | --ini )
#             shift
#             nIni=$1
#             ;;
#         -n | --noise )
#             shift
#             nNoise=$1
#             ;;            
#         -iter ) 
#             shift
#             nIter=$1
#             itermode=1
#             ;;
#         -usp )
#             uspmode=1
#             ;;
#         -csp )
#             cspmode=1
#             ;;
#         -ps )
#             outputtype=ps
#             ;;
#         -svg )
#             outputtype=svg
#             ;;
#         -nc | --no-cuts )
#             write_model_cuts=0
#             ;;
#         -yerrorbars )
#             yerrorbars=1
#             ;;
#         -nmr | --do-not-use-max-radius )
#             use_max_radius=0
#             ;;
#         -h )
#             Usage
#             exit
#             ;;
#         -map | --gt-map )
#             shift
#             GTMapFile=`pwd`/$1
#             gt_map_mode=1
#             ;;
#         * )
#             expBaseName=`basename $1`
#             ;;
#     esac
#     shift
# done

expBaseName=`basename $1`
GTMapFile=`pwd`/$2

# Create results folder.
mkdir -p results

# Column number for indexes.
particleColumn=\$14
micrographColumn=\$17
tiltSeriesColumn=\$8

# Change folder
CSPDIR=$CSPDIR
# cd ${expName}/frealign/maps

# OUTDIR=`mktemp -d -p .` #| sed -e 's:/tmp::'`
OUTDIR=`pwd`/results/compare_all
mkdir -p $OUTDIR
EXPERIMENTS=${OUTDIR}/${expBaseName}_experiments
rm -f ${EXPERIMENTS}
for k in `ls -d ${expBaseName}_ini_??_noise_??`
do
    expDir=${k}/frealign/maps
    # ls ${expDir}/${k}_CSP*.parx
    echo ${k} >> ${EXPERIMENTS}
done
nroExp=`wc -l ${EXPERIMENTS} | awk '{print $1}'`
index=0
#tput civis # Do not show the cursor.
scilab_file=${OUTDIR}/${expBaseName}_SCILAB
rm -f $scilab_file

# Maximum values for ini and noise
maxIni=$(awk 'FS="_" {print $4}' ${EXPERIMENTS} | sort -n -r | head -1)
maxNoise=$(awk 'FS="_" {print $6}' ${EXPERIMENTS} | sort -n -r | head -1)

# Replace the maxNoise in the scilab script that uses it
# (CompareAllExperiments.sci)
sed -e "s|REPLACE_MAX_NOISE|${maxNoise}|g" ${CSPDIR}/scilab/CompareAllExperiments.sci.base > ${CSPDIR}/scilab/CompareAllExperiments.sci

#cat ${EXPERIMENTS}
#exit

# Remove previous intermediate files.
rm -f ${OUTDIR}/${expBaseName}_SCILAB_{CSP,USP}_*

exec 0<"${EXPERIMENTS}"
while read -r line
do
    # echo $line
    expName=`echo $line | awk '{print $1}'`
    expDir=${expName}/frealign/maps
 
    # Groundtruth filename.
    GTParFile=${expBaseName}_01.parx
    # Check existence of the input files 
    if ! [ -f ${GTParFile} ]; then
        echo "[ERROR] `basename $0`: File ${GTParFile} does not exist."
        rm -rf ${OUTDIR}
        exit 1
    fi
    nLines=`wc -l ${GTParFile} | awk '{print $1}'`

    # Get the number of particle, micrograph and tilt series.
    numberOfParticles=`awk '{print '${particleColumn}'}' ${GTParFile} | sort -nr | head -1`
    numberOfMicrographs=`awk '{print '${micrographColumn}'}' ${GTParFile} | sort -nr | head -1`
    numberOfTiltSeries=`awk '{print '${tiltSeriesColumn}'}' ${GTParFile} | sort -nr | head -1`
 
    # Which ini and noise?
    kIni=`echo $expName | sed -e 's:.*ini_\(..\).*:\1:g'`
    kNoise=`echo $expName | sed -e 's:.*noise_\(..\).*:\1:g'`
    
    # Create an array with all the iterations numbers.
    maxCSPIter=`ls ${expDir}/${expName}_CSP*.parx | sed -e 's:'${expDir}/${expName}'_CSP_\(.*\).parx:\1:g' | bc -l | sort -n `
    maxUSPIter=`ls ${expDir}/${expName}_USP*.par | sed -e 's:'${expDir}/${expName}'_USP_\(.*\).par:\1:g' | bc -l | sort -n `
    arrayCSPIter=( $maxCSPIter )
    arrayUSPIter=( $maxUSPIter )

    # Use just the first and last iteration.
    size=${#arrayCSPIter[@]}
    arrayCSPIter=( ${arrayCSPIter[$size-1]} )
    size=${#arrayUSPIter[@]}
    arrayUSPIter=( ${arrayUSPIter[$size-1]} )

    # CSP and USP experiment basenames.
    CSPExpName=${expDir}/${expName}_CSP
    USPExpName=${expDir}/${expName}_USP

    doIterations="1"
    for n in ${doIterations}
    do   
        
        # Iteration numbers, par{x} and fsc files.
        CSPIter=`printf %02d ${arrayCSPIter[${n}-1]}`
        USPIter=`printf %02d ${arrayUSPIter[${n}-1]}`
        CSPParFile=${CSPExpName}_${CSPIter}.parx
        USPParFile=${USPExpName}_${USPIter}.par

        CSPFSCFile=${CSPExpName}_${CSPIter}_fsc.txt
        USPFSCFile=${USPExpName}_${USPIter}_fsc.txt
        
        CSPGTMRCFile=${CSPExpName}_${CSPIter}.mrc
        USPGTMRCFile=${USPExpName}_${USPIter}.mrc
        CSPGTFSCFile=${CSPExpName}_${CSPIter}_gt_fsc.txt
        USPGTFSCFile=${USPExpName}_${USPIter}_gt_fsc.txt

        # Log file where to write the information.
        log_file=${OUTDIR}/${expName}_DATA_`printf %02d $n`
        rm -f $log_file

        # Show progress.
        index=$(( $index + 1 ))
        perc=`echo "$index/$nroExp*100" | bc -l`
        echo -ne "[`printf %3.0f%% $perc`] Processing ${expName}, CSP( ${CSPIter} ) & USP( ${USPIter} ): Extracting data.            \r"

        # Process line by line.
        for k in $(seq 1 1 ${nLines})
        do
            # Tilt series ID was added.
            # awk '{if($1 == '$k') printf "%05d %02d %02d %8.3f %8.3f %8.3f %8.3f %8.3f ", $1, '$micrographColumn', '$particleColumn', $2, $3, $4, $5, $6 }' ${GTParFile} >> ${log_file}
            awk '{if($1 == '$k') printf "%05d %02d %02d %02d %8.3f %8.3f %8.3f %8.3f %8.3f ", $1, '$tiltSeriesColumn', '$micrographColumn', '$particleColumn', $2, $3, $4, $5, $6 }' ${GTParFile} >> ${log_file}
            awk '{if($1 == '$k') printf "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f ", $2, $3, $4, $5, $6, $12, $13 }' ${CSPParFile} >> ${log_file}
            awk '{if($1 == '$k') printf "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f ", $2, $3, $4, $5, $6, $12, $13 }' ${USPParFile} >> ${log_file}
            echo "" >> ${log_file}
        done

        # # Compute the FSC between maps and the GT map.
        # if [ -f "${GTMapFile}" ]; then
            # if [ -f "${CSPGTMRCFile}" ]; then
                # proc3d ${GTMapFile} ${CSPGTMRCFile} apix=4.1 fsc=${CSPGTFSCFile} > /dev/null
            # else
                # echo ${CSPGTMRCFile} does not exists.
            # fi
            # if [ -f "${USPGTMRCFile}" ]; then
                # proc3d ${GTMapFile} ${USPGTMRCFile} apix=4.1 fsc=${USPGTFSCFile} > /dev/null
            # else
                # echo ${USPGTMRCFile} does not exists.
            # fi
        # else
            # echo ${GTMapFile} does not exists.
        # fi

        # Copy the fsc and gt_fsc files
        cp -f $CSPFSCFile ${log_file}_CSP_fsc.txt
        cp -f $USPFSCFile ${log_file}_USP_fsc.txt
        cp -f $CSPGTFSCFile ${log_file}_CSP_gt_fsc.txt
        cp -f $USPGTFSCFile ${log_file}_USP_gt_fsc.txt      

        # Run SciLab for computing metrics and results to process.
        echo -ne "[`printf %3.0f%% $perc`] Processing ${expName}, CSP( ${CSPIter} ) & USP( ${USPIter} ): Running SciLab.             \r"
        sed -e "s|REPLACE_IN_FILE_NAME|${log_file}|" ${CSPDIR}/scilab/CompareAllExperiments.sce.base > ${CSPDIR}/scilab/CompareAllExperiments.sce
        sed -i -e "s|REPLACE_INI|${kIni}|" ${CSPDIR}/scilab/CompareAllExperiments.sce
        sed -i -e "s|REPLACE_NOISE|${kNoise}|" ${CSPDIR}/scilab/CompareAllExperiments.sce
        sed -i -e "s|REPLACE_OUT_FILE_NAME|${scilab_file}|" ${CSPDIR}/scilab/CompareAllExperiments.sce
        scilab -nwni -nb -f "${CSPDIR}/scilab/CompareAllExperiments.sce" 

    done
done
echo -ne "                                                                                                        \r"
#tput cnorm # Show the cursor.

#
# Do the plots.
#
${CSPDIR}/plot_all.sh ${expBaseName} ${OUTDIR}

#
# Remove the temporal folder.
#
# echo "Remove temporal folder (${OUTDIR}) with intermediate files?"
# echo -n "(Plots can be regenerated again with these files and plot_all.sh) [y/n]: "
# read BORRAR
# if [ "$BORRAR" == "y" ]; then
#     rm -rf ${OUTDIR}
# fi
echo "Temporal folder (${OUTDIR}) with intermediate files remains in disk in case you want to regenerate the plots (use plot_all.sh)."
