#!/bin/bash
#
# Version: $Id: compare_experiments.sh 566 2011-02-27 22:25:19Z fefo $
#

##########################################################################
function Usage {

    echo -e "\nNAME"
    echo -e "\t`basename $0`\n"

    echo -e "SYNOPSIS"

    echo -e "\t`basename $0` [OPTIONS] EXPERIMENT_BASENAME -i nIni -n nNoise\n"

    echo "DESCRIPTION"
    echo -e "\tProcess the experiment(s) with log files in the folder"
    echo -e "\t\${EXPERIMENT_BASENAME}_ini_\${nIni}_noise_\${nNoise}. Compare the results with"
    echo -e "\tthe groundtruth values in groundtruth par file located in the same folder."

    echo -e "\nOPTIONS" 

    echo -e "\n\t-iter ITER " 
    echo -e "\t\tProcess only the ITER iteration. (Currently not functional)" 

    echo -e "\n\t-csp" 
    echo -e "\t\tProcess only the CSP experiment. (Currently not functional)" 

    echo -e "\n\t-usp" 
    echo -e "\t\tProcess only the USP experiment. (Currently not functional)" 

    echo -e "\n\t-nc\n\t--no-cuts" 
    echo -e "\t\tDo not write the model cuts." 

    echo -e "\n\t-yerrorbars" 
    echo -e "\t\tPlot the PR with error bars." 

    echo -e "\n\t-nmr\n\t--do-not-use-max-radius"
    echo -e "\t\tDo not scale the polar graph to the max radius, instead use 1." 

    echo -e "\n\t-map FILE\n\t--gt-map FILE"
    echo -e "\t\tGroundtruth density map MRC file." 

}

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

mIni=-1
pIni=-1

# Process options.
while [ "$1" != "" ]
do
    case $1 in
        -im | --ini )
            shift
            mIni=$1
            ;;
        -ip | --ini )
            shift
            pIni=$1
            ;;
        -i | --ini )
            shift
            nIni=$1
            ;;
        -n | --noise )
            shift
            nNoise=$1
            ;;            
        -iter ) 
            shift
            nIter=$1
            itermode=1
            ;;
        -usp )
            uspmode=1
            ;;
        -csp )
            cspmode=1
            ;;
        -ps )
            outputtype=ps
            ;;
        -svg )
            outputtype=svg
            ;;
        -nc | --no-cuts )
            write_model_cuts=0
            ;;
        -yerrorbars )
            yerrorbars=1
            ;;
        -nmr | --do-not-use-max-radius )
            use_max_radius=0
            ;;
        -h )
            Usage
            exit
            ;;
        -map | --gt-map )
            shift
            GTMapFile=`pwd`/$1
            gt_map_mode=1
            ;;
        * )
            expBaseName=`basename $1`
            ;;
    esac
    shift
done

# Process the running mode
if [ "$cspmode" == "-1" -a "$uspmode" == "-1" ]; then
    uspmode=1
    cspmode=1
fi

# Build the experiment name.
if [ "$nIni" == "all" ]
then
    echo Todos los ini.
fi
if [ "$nNoise" == "all" ]
then
    echo Todos los ruidos.
fi
if [ "$nIter" == "all" ]
then
    echo Todas las iteraciones.
fi

if (( $nIni >= 0 ))
then
expName=${expBaseName}_ini_`printf %02d ${nIni}`_noise_`printf %02d ${nNoise}`
else
	expName=${expBaseName}_ini_M_`printf %02d ${mIni}`_P_`printf %02d ${pIni}`_noise_`printf %02d ${nNoise}`
fi

# Create results folder.
mkdir -p results

# Grountruth filename.
GTParFile=`pwd`/${expBaseName}_01.parx

# Check existence of the input files 
if ! [ -f ${GTParFile} ]; then
    echo "[ERROR] `basename $0`: File ${GTParFile} does not exist."
    exit 1
fi
if ! [ -d ${expName} ]; then
    echo "[ERROR] `basename $0`: Folder ${expName} does not exist."
    exit 1
fi
if ! [ -f ${GTMapFile} ]; then
    echo "[ERROR] `basename $0`: File ${GTMapFile} does not exist."
    exit 1
fi

# Column number for indexes.
particleColumn=\$14
micrographColumn=\$17
tiltSeriesColumn=\$8

# Change folder
CSP_FOLDER=$CSPDIR
cd ${expName}/frealign/maps

# Get the number of particle, micrograph and tilt series.
numberOfParticles=`awk '{print '${particleColumn}'}' ${GTParFile} | sort -nr | head -1`
numberOfMicrographs=`awk '{print '${micrographColumn}'}' ${GTParFile} | sort -nr | head -1`
numberOfTiltSeries=`awk '{print '${tiltSeriesColumn}'}' ${GTParFile} | sort -nr | head -1`

# Remove old data previously procesed
rm -f ${expName}_DATA_*
rm -f .max_radius
echo -1 > .max_radius

# CSP and USP should have the same number of iterations.
maxIter=`ls ${expName}_CSP*.parx | sed -e 's/'${expName}'_CSP_\(.*\).parx/\1/g' | bc -l | sort -n | tail -1`
numIter=$(( maxIter - 1 ))
index=0
nLines=`cat ${GTParFile}  | grep -v C | wc -l | awk '{print $1}'`
tput civis # Do not show the cursor.
doIterations=$(seq 2 ${maxIter})
for n in ${doIterations}
do   
 
    # Read the .par{x} file and write the log file.
    iter=`printf %02d $n`
    log_file=${expName}_DATA_$iter
    rm -f ${log_file}
    if ! [ -f "${log_file}" ]
    then
        CSPParFile=${expName}_CSP_$iter.parx
        USPParFile=${expName}_USP_$iter.par
		if ! [ -f ${USPParFile} ]
		then
			USPParFile=$CSPParFile
		fi
        for k in $(seq 1 1 ${nLines})
        do
            index=$(( $index + 1 ))
            perc=`echo "$index/$nLines/$(( maxIter - 1 ))*100" | bc -l`
            echo -ne "[`printf %3.0f%% $perc`] Processing iteration $iter: `printf %5d $k` of ${nLines}     \r"
            
            awk '{if($1 == '$k') printf "%05d %02d %02d %02d %8.3f %8.3f %8.3f %8.3f %8.3f ", $1, '$tiltSeriesColumn', '$micrographColumn', '$particleColumn', $2, $3, $4, $5, $6 }' ${GTParFile} >> ${log_file}
            awk '{if($1 == '$k') printf "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f ", $2, $3, $4, $5, $6, $12, $13 }' ${CSPParFile} >> ${log_file}
            awk '{if($1 == '$k') printf "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f ", $2, $3, $4, $5, $6, $12, $13 }' ${USPParFile} >> ${log_file}
            echo "" >> ${log_file}
        done
    fi
    index=$(( n * nLines - nLines ))
    
    # Run SciLab for computing metrics and results to process.
    rm -f ${log_file}_*
    perc=`echo "$index/$nLines/9*100" | bc -l`
    echo -ne "[`printf %3.0f%% $perc`] Processing iteration $iter: Running SciLab                          \r"
    sed -e "s|REPLACE_EXP_NAME|${expName}|" ${CSP_FOLDER}/scilab/CompareExperiments.sce.base > CompareExperiments.sce
    sed -i -e "s|REPLACE_ITER_NUMBER|${n}|" CompareExperiments.sce
    # if [ "${CLUSTER_ID}" == "NIH" ]; then export DISPLAY=biowulf.nih.gov:51.0; fi
	scilab -nwni -nb -f "CompareExperiments.sce" 
	rm -f CompareExperiments.sce
    
done
tput cnorm # Show the cursor.

pixel_size=`header ${GTMapFile} | grep spacing | awk '{print $4}'`

# Compute the FSC between maps and the GT map.
if [ "$gt_map_mode" == "1" ]
then
    for XSP in CSP USP
    do
        doIterations=$(seq 2 1 ${maxIter})
        doIterations="2 ${maxIter}"
        for k in ${doIterations}
        do
            xsp_map=${expName}_${XSP}_`printf %02d ${k}`.mrc
            if [ -f "${xsp_map}" ]; then
				# map is not neccesarily aligned to GT map, we need to align them before computing the FSC
				
				aligned=${expName}_${XSP}_`printf %02d ${k}`_aligned.mrc
				# align GT to current map
				com="${HOME}/code/ETTK/AlignRigid3D ${xsp_map} ${GTMapFile} 10 180 ${aligned}"
				echo $com; echo; $com
                fsc_file=${expName}_${XSP}_`printf %02d ${k}`_gt_fsc.txt
                proc3d ${aligned} ${xsp_map} apix=${pixel_size} fsc=${fsc_file} > /dev/null
                #proc3d ${GTMapFile} ${xsp_map} apix=${pixel_size} fsc=${fsc_file} > /dev/null
				# rm -f aligned.mrc
            fi
        done
    done
fi

#
# Do the plots.
#
com="${CSP_FOLDER}/plot.sh ${expName} ${maxIter} ${outputtype} ${write_model_cuts} ${use_max_radius} ${yerrorbars} ${gt_map_mode}"
echo $com; $com

# Get back...
cd - > /dev/null
