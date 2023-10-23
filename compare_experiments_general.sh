#!/bin/bash
#
# Version: $Id: compare_experiments.sh 531 2010-10-20 23:10:53Z fefo $
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

expName=${expBaseName}

if ! [ -d maps ]; then
    echo "[ERROR] `basename $0`: Folder maps does not exist."
    exit 1
fi

if ! [ -f ${GTMapFile} ]; then
    echo "[ERROR] `basename $0`: File ${GTMapFile} does not exist."
    exit 1
fi

# Process the running mode
if [ "$cspmode" == "-1" -a "$uspmode" == "-1" ]; then
    uspmode=1
    cspmode=1
fi

# Column number for indexes.
particleColumn=\$14
micrographColumn=\$17
tiltSeriesColumn=\$8

# Change folder
CSP_FOLDER=$CSPDIR
cd maps

GTParFile=../${expName}_01.parx

echo Data statistics:

# Get the number of particle, micrograph and tilt series.
numberOfTiltSeries=`awk '{print '${tiltSeriesColumn}'}' ${GTParFile} | sort -nr | head -1`

# numberOfParticles=`awk '{print '${particleColumn}'}' ${GTParFile} | sort -nr | head -1`
numberOfParticles=0
for (( i = 0; i < $numberOfTiltSeries; i++ ))
do
	particles=`awk '{if('${tiltSeriesColumn}' == '$i') print '${particleColumn}'}' ${GTParFile} | sort -nr | head -1`
	let numberOfParticles=numberOfParticles+particles
done

# numberOfMicrographs=`awk '{print '${micrographColumn}'}' ${GTParFile} | sort -nr | head -1`
numberOfMicrographs=0
for (( i = 0; i < $numberOfTiltSeries; i++ ))
do
	micrographs=`awk '{if('${tiltSeriesColumn}' == '$i') print '${micrographColumn}'}' ${GTParFile} | sort -nr | head -1`
	let numberOfMicrographs=numberOfMicrographs+micrographs
done

echo -e "\t"Number of particles = $numberOfParticles
echo -e "\t"Number of micrographs = $numberOfMicrographs
echo -e "\t"Number of tilt series = $numberOfTiltSeries

# Remove old data previously procesed
rm -f ${expName}_DATA_*
rm -f .max_radius
echo -1 >  .max_radius

# CSP and USP should have the same number of iterations.
maxIter=`ls ${expName}_CSP_??.parx | sed -e 's/'${expName}'_CSP_\(.*\).parx/\1/g' | bc -l | sort -n | tail -1`
maxIter=25
echo maxIter = $maxIter
numIter=$(( maxIter - 1 ))
index=0
nLines=`cat ${GTParFile}  | grep -v C | wc -l | awk '{print $1}'`
echo nLines = $nLines
tput civis # Do not show the cursor.
doIterations=$(seq 2 ${maxIter})
for n in ${doIterations}
do   
 
    # Read the .par{x} file and write the log file.
    iter=`printf %02d $n`
    log_file=${expName}_DATA_$iter
    # rm -f ${log_file}
    if [ -f "${log_file}" ]
    then
		echo Information for iteration $iter already exists, skipping.
	else
        CSPParFile=${expName}_CSP_$iter.parx
        USPParFile=${expName}_USP_$iter.par
		if ! [ -f ${USPParFile} ]
		then
			USPParFile=$CSPParFile
		fi
		# for k in $(seq 1 1 ${nLines})
        # do
            # index=$(( $index + 1 ))
            # perc=`echo "$index/$nLines/$(( maxIter - 1 ))*100" | bc -l`
            # echo -ne "[`printf %3.0f%% $perc`] Processing iteration $iter: `printf %5d $k` of ${nLines}     \r"
            
            # awk '{if($1 == '$k') printf "%06d %02d %02d %02d %8.3f %8.3f %8.3f %8.3f %8.3f ", $1, '$tiltSeriesColumn', '$micrographColumn', '$particleColumn', $2, $3, $4, $5, $6 }' ${GTParFile} >> ${log_file}
			# awk '{if($1 == '$k') printf "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f ", $2, $3, $4, $5, $6, $12, $13 }' ${CSPParFile} >> ${log_file}
            # awk '{if($1 == '$k') printf "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f ", $2, $3, $4, $5, $6, $12, $13 }' ${USPParFile} >> ${log_file}
            # echo "" >> ${log_file}
        # done
		# faster version
		cat ${GTParFile} | grep -v C | awk '{ printf "%06d %02d %02d %02d %8.3f %8.3f %8.3f %8.3f %8.3f\n", $1, '$tiltSeriesColumn', '$micrographColumn', '$particleColumn', $2, $3, $4, $5, $6 }' > ${log_file}_1
		cat ${CSPParFile} | grep -v C | awk '{ printf "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", $2, $3, $4, $5, $6, $12, $13 }' > ${log_file}_2
		cat ${USPParFile} | grep -v C | awk '{ printf "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", $2, $3, $4, $5, $6, $12, $13 }' > ${log_file}_3
		paste ${log_file}_1 ${log_file}_2 ${log_file}_3 -d " " > ${log_file}
		rm -f ${log_file}_1 ${log_file}_2 ${log_file}_3
    fi
    index=$(( n * nLines - nLines ))
    
    # Run SciLab for computing metrics and results to process.
    rm -f ${log_file}_*
    perc=`echo "$index/$nLines/9*100" | bc -l`
    echo -ne "[`printf %3.0f%% $perc`] Processing iteration $iter: Running SciLab                          \r"
    sed -e "s|REPLACE_EXP_NAME|${expName}|" ${CSP_FOLDER}/scilab/CompareExperiments.sce.base > CompareExperiments.sce
    sed -i -e "s|REPLACE_ITER_NUMBER|${n}|" CompareExperiments.sce
	# export DISPLAY=biowulf.nih.gov:19.0
	scilab -nwni -nb -f "CompareExperiments.sce" 
	rm -f CompareExperiments.sce

	# compute incremental fsc's
	pixel_size=`header ${expName}_CSP_${iter}.mrc | grep spacing | awk '{print $4}'`
	if [ $n -lt $maxIter ]
	then
		let nextiter=n+1
		com="proc3d ${expName}_CSP_$iter.mrc ${expName}_CSP_`printf %02d $nextiter`.mrc apix=${pixel_size} fsc=${expName}_fsc_${iter}_`printf %02d $nextiter`.txt"
		echo $com; $com
	fi
	
done
tput cnorm # Show the cursor.

output=${expName}_DATA.txt
rm -f $output

# Compute the FSC between maps and the GT map.
if [ "$gt_map_mode" == "1" ]
then
	pixel_size=`header ${GTMapFile} | grep spacing | awk '{print $4}'`
    for XSP in CSP
    do
        doIterations=$(seq 2 1 ${maxIter})
        # doIterations="2 10"
        for k in ${doIterations}
        do
            xsp_map=${expName}_${XSP}_`printf %02d ${k}`.mrc
            if [ -f "${xsp_map}" ]; then
				# map is not neccesarily aligned to GT map, we need to align them before computing the FSC
				
				fsc_file=${expName}_${XSP}_`printf %02d ${k}`_gt_fsc.txt

				if ! [ -f ${fsc_file} ]
				then
					# align GT to current map
					~/code/ETTK/AlignRigid3D ${xsp_map} ${GTMapFile} 10 180 aligned.mrc
				
					proc3d aligned.mrc ${xsp_map} apix=${pixel_size} fsc=${fsc_file} > /dev/null
					#proc3d ${GTMapFile} ${xsp_map} apix=${pixel_size} fsc=${fsc_file} > /dev/null
					rm -f aligned.mrc
				fi
            fi

			# get all resolution measures 
			map=${expName}_${XSP}_`printf %02d ${k}`
			fsc_file=${map}_fsc.txt
			fsc_gt_file=${map}_gt_fsc.txt
			res=`${SPA_DIR}/general/fsc_cutoff.sh ${fsc_file} | awk '{print $5}'`
			resgt=`${SPA_DIR}/general/fsc_cutoff.sh ${fsc_gt_file} | awk '{print $5}'`
			rmeasure=`${SPA_DIR}/utils/rmeasure.sh ${map}.mrc | grep "Resolution at FSC = 0.5" | awk '{print $6}'`
			echo -e "$k\t$res\t$resgt\t$rmeasure" >> $output

        done
    done
fi

#
# Do the plots.
#
com="${CSP_FOLDER}/plot_general.sh ${expName} ${maxIter} ${outputtype} ${write_model_cuts} ${use_max_radius} ${yerrorbars} ${gt_map_mode}"
echo $com; $com

# Get back...
cd - > /dev/null
