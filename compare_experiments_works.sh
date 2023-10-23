#!/bin/bash
#
# Version: $Id: compare_experiments.sh 522 2010-10-14 21:48:32Z fefo $
#

##########################################################################
function Usage {

    echo -e "\nNAME"
    echo -e "\t`basename $0`\n"

    echo -e "SYNOPSIS"

    echo -e "\t`basename $0` [OPTIONS] EXPERIMENT_FOLDER\n"

    echo "DESCRIPTION"
    echo -e "\tProcess the experiment(s) with log files in the folder EXPERIMENT_FOLDER. Compare"
    echo -e "\tthe results with the groundtruth values in groundtruth par file located in the"
    echo -e "\tsame folder."

    echo -e "\nOPTIONS" 

    echo -e "\n\t-n ITER " 
    echo -e "\t\tProcess only the ITER iteration. (Currently not functional)" 

    echo -e "\n\t-csp" 
    echo -e "\t\tProcess only the CSP experiment. (Currently not functional)" 

    echo -e "\n\t-usp" 
    echo -e "\t\tProcess only the USP experiment. (Currently not functional)" 

    echo -e "\n\t-no_cuts" 
    echo -e "\t\tDo not write the model cuts." 

    echo -e "\n\t-yerrorbars" 
    echo -e "\t\tPlot the PR with error bars." 

}

# With no arguments show the usage help
if [ "$#" == "0" ]; then Usage; exit; fi

outputtype=png
experimentName=${!#}
write_model_cuts=1
yerrorbars=0
itermode=0
uspmode=-1
cspmode=-1
while [ "$1" != "" ]
do
    case $1 in
        -n ) 
            shift
            iter=$1
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
        -no_cuts )
            write_model_cuts=0
            ;;
        -yerrorbars )
            yerrorbars=1
            ;;
        -h )
            Usage
            exit
            ;;
        * )
            expName=`basename $1`
            ;;
    esac
    shift
done

# Process the running mode
if [ "$cspmode" == "-1" -a "$uspmode" == "-1" ]; then
    uspmode=1
    cspmode=1
fi

GTParFile=${expName}_01.parx

# Check existence of the input files 
if ! [ -f ${expName}/frealign/maps/${GTParFile} ]; then
    echo "[ERROR] `basename $0`: File ${expName}/frealign/maps/${GTParFile} does not exist."
    exit 1
fi
if ! [ -d ${expName} ]; then
    echo "[ERROR] `basename $0`: Folder ${expName} does not exist."
    exit 1
fi

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
# echo Number of Tilt Series: $numberOfTiltSeries
# echo Number of Micrographs: $numberOfMicrographs 
# echo Number of particles: $numberOfParticles 

# Remove old data previously procesed
rm -f data_iter_XX_*
rm -f .max_radius
echo -1 >  .max_radius

# CSP and USP should have the same number of iterations.
maxIter=`ls ${expName}_CSP*.parx | sed -e 's/'${expName}'_CSP_\(.*\).parx/\1/g' | bc -l | sort -n | tail -1`
numIter=$(( maxIter - 1 ))
index=0
nLines=`wc -l ${GTParFile} | awk '{print $1}'`
tput civis
doIterations=$(seq 2 1 ${maxIter})
for n in ${doIterations}
do    
    iter=`printf %02d $n`
    log_file=data_iter_$iter
    # rm -f ${log_file}
    if ! [ -f "${log_file}" ]
    then
        CSPParFile=${expName}_CSP_$iter.parx
        USPParFile=${expName}_USP_$iter.par
        for k in $(seq 1 1 ${nLines})
        do
            index=$(( $index + 1 ))
            perc=`echo "$index/$nLines/$(( maxIter - 1 ))*100" | bc -l`
            echo -ne "[`printf %3.0f%% $perc`] Processing iteration $iter: `printf %5d $k` of ${nLines}     \r"
            
            awk '{if($1 == '$k') printf "%05d %02d %02d %8.3f %8.3f %8.3f %8.3f %8.3f ", $1, '$micrographColumn', '$particleColumn', $2, $3, $4, $5, $6 }' ${GTParFile} >> ${log_file}
            awk '{if($1 == '$k') printf "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f ", $2, $3, $4, $5, $6, $12, $13 }' ${CSPParFile} >> ${log_file}
            awk '{if($1 == '$k') printf "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f ", $2, $3, $4, $5, $6, $12, $13 }' ${USPParFile} >> ${log_file}
            echo "" >> ${log_file}
        done
    fi
    index=$(( n * nLines - nLines ))
    
    # Run SciLab for computing metrics and results to process.
    # echo "Run SciLab"
    rm -f ${log_file}_*
    perc=`echo "$index/$nLines/9*100" | bc -l`
    echo -ne "[`printf %3.0f%% $perc`] Processing iteration $iter: Running SciLab                          \r"
    sed -e "s|REPLACE_ITER_NUMBER|${n}|" ${CSP_FOLDER}/scilab/CompareExperiments.sce.base > ${CSP_FOLDER}/scilab/CompareExperiments.sce
    scilab -nw -nb -f "${CSP_FOLDER}/scilab/CompareExperiments.sce" 
    
done
tput cnorm

# Concatenate all frobenius
toutA=`mktemp`
toutA2=`mktemp`
cp data_iter_02_frobenius_transformation ${toutA}
doIterations=$(seq 3 1 ${maxIter})
for k in ${doIterations}
do
    join ${toutA} data_iter_`printf %02d $k`_frobenius_transformation > ${toutA2}
    cp ${toutA2} ${toutA}
done
cp ${toutA} data_iter_XX_frobenius_transformation

# Compute the histogram values for the Frobenius norm of the
# difference in the transformation matrices.
nBins=10
sed -e "s|REPLACE_NCOLS|$(( 2*maxIter - 1 ))|" ${CSP_FOLDER}/scilab/FrobeniusHistogram.sce.base > ${CSP_FOLDER}/scilab/FrobeniusHistogram.sce
sed -i -e "s|REPLACE_NBINS|${nBins}|" ${CSP_FOLDER}/scilab/FrobeniusHistogram.sce
scilab -nw -nb -f "${CSP_FOLDER}/scilab/FrobeniusHistogram.sce" 

#
# Plots
#

# Set title and output type 
title=${expName}
if [ "${outputtype}" == "ps" ]
then
    outputFile=../../../results/${expName}_results.ps
    gnustr="set term postscript enhanced color eps font 'Helvetiva,10'"
    gnustr="${gnustr}; set output '${outputFile}'"
else  
    outputFile=../../../results/${expName}_results.png
    gnustr="set terminal png font '/usr/share/fonts/truetype/msttcorefonts/arial.ttf,10' size 1200, 800"
    gnustr="${gnustr}; set output '${outputFile}'"
    # gnustr="${gnustr}; set output '| display png:-'"
fi
gnustr="${gnustr}; set multiplot layout 2, 3 title '${title} results'"

# FSC curves
title="FSC curves"
xlabel="Resolution"
ylabel=""
str="set title '${title}'"
str="${str}; set xlabel '${xlabel}'"
str="${str}; set ylabel '${ylabel}'"
str="${str}; plot "
doIterations=$(seq 2 1 ${maxIter})
doIterations="2 ${maxIter}"
for k in ${doIterations}
do 
    int=$(( 255-255*${k}/${maxIter} ))
    hexx=`printf %02X $int`
    rgb="#${hexx}FF${hexx}"
    lw=$(( 3*$k/10 ))
    str="${str} '${expName}_USP_`printf %02d $k`_fsc.txt' using 1:2 with lines lw $lw lc rgbcolor '$rgb' "
    if [ "$k" == "$maxIter" ]; then 
        str="${str} ti \"USP\", "; 
    else
        str="${str} ti \"\", "; 
    fi
done
for k in ${doIterations}
do 
    int=$(( 255-255*${k}/${maxIter} ))
    hexx=`printf %02X $int`
    rgb="#FF${hexx}${hexx}"
    lw=$(( 3*$k/10 ))
    str="${str} '${expName}_CSP_`printf %02d $k`_fsc.txt' using 1:2 with lines lw $lw lc rgbcolor '$rgb' "
    if [ "$k" == "$maxIter" ]; then 
        str="${str} ti \"CSP\" "; 
    else
        str="${str} ti \"\", "; 
    fi
    # if ! [ "$k" == "$maxIter" ]; then str="${str},"; fi
done
gnustr="${gnustr}; ${str}"

# Global PR
sed -i 1d data_iter_XX_global_pr # remove fisrt line
title="Mean Phase Residual per Iteration"
xlabel=Iterations
ylabel="Mean Phase Residual"
file=data_iter_XX_global_pr
str="set title '${title}'"
str="${str}; set xlabel '${xlabel}'"
str="${str}; set ylabel '${ylabel}'"
str="${str}; set style data linespoints"
if [ "$yerrorbars" == "0" ]; then
    str="${str}; plot '${file}' using 1:2 pt 5 ti \"CSP\", '${file}' using 1:5 pt 5 ti \"USP\""
else
    str="${str}; plot '${file}' using 1:2:3:4 with yerrorbars pt 5 ti \"CSP\", '${file}' using 1:5:6:7 with yerrorbars pt 5 ti \"USP\""
fi
gnustr="${gnustr}; ${str}"

# Global DPR
sed -i 1d data_iter_XX_global_dpr # remove fisrt line
title="Mean Differential Phase Residual"
xlabel=Iterations
ylabel="Mean Differential Phase Residual per Iteration"
file=data_iter_XX_global_dpr
str="set title '${title}'"
str="${str}; set xlabel '${xlabel}'"
str="${str}; set ylabel '${ylabel}'"
str="${str}; set style data linespoints"
str="${str}; plot '${file}' using 1:2 pt 5 ti \"CSP\", "
str="${str} '${file}' using 1:3 pt 5 ti \"USP\""
gnustr="${gnustr}; ${str}"


# Frobenius of Transformation Matrices
title="Frobenius of Difference in Transformation Matrices"
xlabel="Particle Projection Index"
ylabel=""
file=data_iter_XX_frobenius_transformation
str="set title '${title}'"
str="${str}; set logscale y"
str="${str}; set xlabel '${xlabel}'"
str="${str}; set ylabel '${ylabel}'"
str="${str}; set style data points"
str="${str}; set key title 'CSP iter'"
str="${str}; plot "
# doIterations=$(seq 2 1 ${maxIter})
doIterations="2 ${maxIter}"
for k in ${doIterations}
do 
    # str="${str} '${file}' using 1:$(( 2*k - 2 )) lt ${k} pt ${k} ti \"Iter ${k}\" "
    str="${str} '${file}' using 1:$(( 2*k - 2 )) pt 5 ps .2 ti \"${k}\" "
    # str="${str} '${file}' using 1:$(( 2*k - 2 )) lt 1 pt 1 ti \"\", "
    # str="${str} '${file}' using 1:$(( 2*k - 1 )) lt 2 pt 2 ti \"\""
    if [ "${k}" -lt "${maxIter}" ]; then str="${str}, "; fi
done
gnustr="${gnustr}; ${str}"

# Histogram of Frobenius of Difference Transformation Matrices
title="Histogram of Frobenius of Difference in Transformation Matrices"
xlabel=""
ylabel=""
str="set title '${title}'"
str="${str}; set xlabel '${xlabel}'"
str="${str}; set ylabel '${ylabel}'"
str="${str}; unset logscale y"
str="${str}; set key title 'CSP iter'"
str="${str}; set boxwidth 0.9 absolute"
str="${str}; set style fill solid 1.00 noborder"
str="${str}; set style histogram clustered gap 1 title offset character 0, 0, 0"
str="${str}; set style data histograms"
str="${str}; set xtics border in scale 1,0.5 nomirror rotate by -45 offset character 0, 0, 0"
str="${str}; set format x \"%2.0te%L\""
str="${str}; plot "
file=data_iter_XX_frobenius_histogram
doIterations=$(seq 2 1 ${maxIter})
doIterations="2 6 ${maxIter}"
for k in ${doIterations}
do 
    str="${str} '${file}' using ${k}:xtic(1) ti '$k'"
    if [ "${k}" -lt "${maxIter}" ]; then str="${str}, "; fi
done
str="${str}; unset logscale y"
gnustr="${gnustr}; ${str}"

# Plot angles differences
max_radius=`cat .max_radius`
xtics=`echo  "0.5 * $max_radius" | bc -l`
yposition=`echo  "0.9 * $max_radius" | bc -l`
max_radius=`echo  "1.1 * $max_radius" | bc -l`
title="Angles difference"
str="set title '${title}'"
str="${str}; set key title ''"
str="${str}; unset border"
str="${str}; set angles degrees"
str="${str}; set grid polar 15"
str="${str}; set polar"
str="${str}; set size square"
str="${str}; set style data points"
str="${str}; set xtics ${xtics} norotate"
str="${str}; unset ytics"
str="${str}; set format x \"%2.0te%L\""
# str="${str}; set yrange [-1:1]"
# str="${str}; set xrange [-1:1]"
str="${str}; set yrange [-${max_radius}:${max_radius}]"
str="${str}; set xrange [-${max_radius}:${max_radius}]"
str="${str}; set label 'radius=sin(theta) vs. angle=phi' at -${max_radius},-${yposition}"
str="${str}; plot "
# doIterations=$(seq 2 1 ${maxIter})
doIterations="2 ${maxIter}"
for k in ${doIterations}
do 
    file=data_iter_`printf %02d $k`_angles_difference_polar
    str="${str} '${file}' using 4:3 pt 4 ps .7 ti \"CSP iter $k\", "
    str="${str} '${file}' using 7:6 pt 4 ps .7 ti \"USP iter $k\""
    if [ "${k}" -lt "$maxIter" ]; then str="${str}, "; fi
done
gnustr="${gnustr}; ${str}"

# Run gnuplot
command=`echo "${gnustr}" | gnuplot`
echo ${command}
$command 

# Obtained the axial, sagital and cornal cut for every density
# map. Then concatenate them in one big image.
if [ "$write_model_cuts" == "1" ]
then 
    tput civis
    for XSP in CSP USP
    do
        doIterations=$(seq 2 1 ${maxIter})
        for k in ${doIterations}
        do
            iter=`printf %02d $k`
            model=${expName}_${XSP}_${iter}.mrc
            echo -ne "Generating model cuts for ${model}\r"
            axial=${expName}_${XSP}_${iter}_axial
            coronal=${expName}_${XSP}_${iter}_coronal
            sagital=${expName}_${XSP}_${iter}_sagital
            
            newstack ${model} ${axial}.mrc -secs 50 -scale 0,255 -mode 0 > /dev/null
            mrc2tif ${axial}.mrc ${axial}.tif > /dev/null
            convert -normalize -crop 50x50+25+25 ${axial}.tif ${axial}.png
            
            rotatevol -input ${model} -output ${sagital}.mrc -angles 0,90,0 > /dev/null
            newstack ${sagital}.mrc ${sagital}.mrc -secs 50 -scale 0,255 -mode 0 > /dev/null
            mrc2tif ${sagital}.mrc ${sagital}.tif > /dev/null
            convert -normalize -rotate 90 -crop 50x50+25+25 -negate ${sagital}.tif ${sagital}.png
            convert -negate ${sagital}.png ${sagital}.png
            
            rotatevol -input ${model} -output ${coronal}.mrc -angles 0,90,90 > /dev/null
            newstack ${coronal}.mrc ${coronal}.mrc -secs 50 -scale 0,255 -mode 0 > /dev/null
            mrc2tif ${coronal}.mrc ${coronal}.tif > /dev/null
            convert -normalize -rotate 90 -crop 50x50+25+25 -negate ${coronal}.tif ${coronal}.png
            convert -negate ${coronal}.png ${coronal}.png
            
            rm *.mrc~
            montage -tile 3x -geometry 100x100 ${expName}_${XSP}_${iter}_{axial,coronal,sagital}.png ${XSP}_${iter}.png
        done
        montage -tile 1x -geometry 300x100 ${XSP}_*.png ${XSP}.png
        rm -f ${XSP}_*.png *_{axial,coronal,sagital}.*
    done
    montage -tile 2x1 -geometry 300x$(( 100 * numIter)) CSP.png USP.png ${CSP_FOLDER}/results/${expName}_model_cuts.png
    rm -f CSP.png USP.png
    echo ""
    tput cnorm
fi

# Move the results to the results folder
#mv ${outputFile} ../../../results
echo "Results plotted in ${outputFile}"

# Get back...
cd - > /dev/null

# Show the results
if [ "${outputtype}" == "ps" ]
then
    evince ${outputFile} &
else  
    display ${outputFile} &
fi

# convert -density 150 plot.eps plot.png
