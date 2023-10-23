#!/bin/bash
#
# Version: $Id: plot.sh 566 2011-02-27 22:25:19Z fefo $
#

expName=$1
maxIter=$2
outputtype=$3
write_model_cuts=$4
use_max_radius=$5
yerrorbars=$6
gt_map_mode=$7


if [ -n $8 ]
then
	resultsdir=../../../results/
else
	resultsdir=results
fi

#
# Do the plots.
#

# Set title and output type 
title=${expName}
font=/usr/share/fonts/truetype/ttf-dejavu/DejaVuSans.ttf
fontsize=9
width=1200
height=800
if [ "${outputtype}" == "ps" ]
then
    outputFile=${resultsdir}/${expName}_results.ps
    gnustr="set term postscript enhanced color eps font 'Helvetiva,10'"
    gnustr="${gnustr}; set output '${outputFile}'"
elif [ "${outputtype}" == "svg" ]
then
    outputFile=${resultsdir}/${expName}_results.svg
    gnustr="set term svg enhanced font '${font},${fontsize}' size ${width}, ${height}"
    gnustr="${gnustr}; set output '${outputFile}'"
else  
    outputFile=${resultsdir}/${expName}_results.png
    gnustr="set terminal png font '${font},${fontsize}' size ${width}, ${height}"
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
str="${str}; set yrange [0:1.05]"
str="${str}; plot "
doIterations=$(seq 2 1 ${maxIter})
str="${str} '${expName}_USP_01_fsc.txt' using 1:2 with lines lw 3 lc rgbcolor \"#0000FF\" ti \"ini\", "
doIterations="2 ${maxIter}"
lw=1
for k in ${doIterations}
do 
    int=$(( 255-255*${k}/${maxIter} ))
    hexx=`printf %02X $int`
	hexx=00
    #lw=$(( 3*$k/15 ))
    # if [ "$gt_map_mode" == "1" ]
    # then    
    #     rgb="#${hexx}${hexx}FF"
    #     str="${str} '${expName}_USP_`printf %02d $k`_gt_fsc.txt' using 1:2 with lines lw $lw lc rgbcolor '$rgb' "
    #     if [ "$k" == "$maxIter" ]; then 
    #         str="${str} ti \"USP vs. GT\", "; 
    #     else
    #         str="${str} ti \"\", "; 
    #     fi
    # fi
    rgb="#${hexx}FF${hexx}"
    str="${str} '${expName}_USP_`printf %02d $k`_fsc.txt' using 1:2 with lines lw $lw lc rgbcolor '$rgb' "
    if [ "$k" == "$maxIter" ]; then 
        str="${str} ti \"USP\", "; 
    else
        str="${str} ti \"\", "; 
    fi
	let lw=lw+1
done
lw=1
for k in ${doIterations}
do 
    int=$(( 255-255*${k}/${maxIter} ))
    hexx=`printf %02X $int`
	hexx=00
    #lw=$(( 3*$k/15 ))
    # if [ "$gt_map_mode" == "1" ]
    # then
    #     rgb="#FF${hexx}FF"
    #     str="${str} '${expName}_CSP_`printf %02d $k`_gt_fsc.txt' using 1:2 with lines lw $lw lc rgbcolor '$rgb' "
    #     if [ "$k" == "$maxIter" ]; then 
    #         str="${str} ti \"CSP vs. GT\", "; 
    #     else
    #         str="${str} ti \"\", "; 
    #     fi
    # fi
    rgb="#FF${hexx}${hexx}"
    str="${str} '${expName}_CSP_`printf %02d $k`_fsc.txt' using 1:2 with lines lw $lw lc rgbcolor '$rgb' "
    if [ "$k" == "$maxIter" ]; then 
        str="${str} ti \"CSP\" "; 
    else
        str="${str} ti \"\", "; 
    fi
	let lw=lw+1
done
str="${str}; set yrange [*:*]"
gnustr="${gnustr}; ${str}"

# Global PR
file=${expName}_DATA_PR
# sed -i 1d ${file} # remove fisrt line
title="Mean Phase Residual per Iteration"
xlabel=Iterations
ylabel="Mean Phase Residual"
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
file=${expName}_DATA_DPR
sed -i 1d ${file} # remove fisrt line
title="Mean Differential Phase Residual"
xlabel=Iterations
ylabel="Mean Differential Phase Residual per Iteration"
str="set title '${title}'"
str="${str}; set xlabel '${xlabel}'"
str="${str}; set ylabel '${ylabel}'"
str="${str}; set style data linespoints"
str="${str}; plot '${file}' using 1:2 pt 5 ti \"CSP\", "
str="${str} '${file}' using 1:3 pt 5 ti \"USP\""
gnustr="${gnustr}; ${str}"

# FSC curves
if [ "$gt_map_mode" == "1" ]
then    
    title="FSC curves"
    xlabel="Resolution"
    ylabel=""
    str="set title '${title}'"
    str="${str}; set xlabel '${xlabel}'"
    str="${str}; set ylabel '${ylabel}'"
    str="${str}; set yrange [0:1.05]"
    str="${str}; plot "
    doIterations=$(seq 2 1 ${maxIter})
    doIterations="2 ${maxIter}"
	lw=1
    for k in ${doIterations}
    do 
        int=$(( 255-255*${k}/${maxIter} ))
        hexx=`printf %02X $int`
		hexx=00
        #lw=$(( 3*$k/15 ))
        rgb="#${hexx}FF${hexx}"
        str="${str} '${expName}_USP_`printf %02d $k`_gt_fsc.txt' using 1:2 with lines lw $lw lc rgbcolor '$rgb' "
        if [ "$k" == "$maxIter" ]; then 
            str="${str} ti \"USP vs. GT\", "; 
        else
            str="${str} ti \"\", "; 
        fi
		let lw=lw+1
    done
    lw=1
    for k in ${doIterations}
    do 
        int=$(( 255-255*${k}/${maxIter} ))
        hexx=`printf %02X $int`
        hexx=00
		#lw=$(( 3*$k/15 ))
        rgb="#FF${hexx}${hexx}"
        str="${str} '${expName}_CSP_`printf %02d $k`_gt_fsc.txt' using 1:2 with lines lw $lw lc rgbcolor '$rgb' "
        if [ "$k" == "$maxIter" ]; then 
            str="${str} ti \"CSP vs. GT\" "; 
        else
            str="${str} ti \"\", "; 
        fi
		let lw=lw+1
    done
    str="${str}; set yrange [*:*]"
    gnustr="${gnustr}; ${str}"
fi

# Frobenius of Transformation Matrices
title="Frobenius of Difference in Transformation Matrices"
xlabel="Particle Projection Index"
ylabel=""
str="set title '${title}'"
str="${str}; set logscale y"
str="${str}; set xlabel '${xlabel}'"
str="${str}; set ylabel '${ylabel}'"
str="${str}; set style data points"
# str="${str}; set key title 'CSP iter'"
str="${str}; plot "
doIterations="2 ${maxIter}"
index=1
lw=1
for k in ${doIterations}
do 

    int=$(( 255-255*${k}/${maxIter} ))
    hexx=`printf %02X $int`
    # lw=1 #$(( 3*$k/15 ))

    titer=`printf %02d $k`
    file1=${expName}_DATA_FROB_${titer}
    file2=${expName}_DATA_FROB_SORT_${titer}

    # USP

    # # Points
    # int=$(( 255-255*${k}/${maxIter}/2 ))
    # hexx=`printf %02X $int`
    # rgb="#${hexx}FF${hexx}"
    # str="${str} '${file1}' using 1:3 pt 2 ps .5 lt ${index} lc rgbcolor '$rgb' ti \"USP@${k}\", "

    int=$(( 255-255*${k}/${maxIter} ))
    hexx=00 # `printf %02X $int`
    rgb="#${hexx}FF${hexx}"
    str="${str} '${file2}' using 1:3 with lines lw $lw lc rgbcolor '$rgb' "
    index=$(( index + 1 ))

    if [ "$k" == "$maxIter" ]; then 
        str="${str} ti \"USP\", " ; 
    else
        str="${str} ti \"\", "; 
    fi

    # CSP

    # # Points
    # int=$(( 255-255*${k}/${maxIter}/2 ))
    # hexx=`printf %02X $int`
    # rgb="#FF${hexx}${hexx}"
    # str="${str} '${file1}' using 1:2 pt 4 ps .5 lt ${index} lc rgbcolor '$rgb' ti \"CSP@${k}\", "

    int=$(( 255-255*${k}/${maxIter} ))
    hexx=00 # `printf %02X $int`
    rgb="#FF${hexx}${hexx}"
    str="${str} '${file2}' using 1:2 with lines lw $lw lc rgbcolor '$rgb' "
    # if [ "${k}" -lt "${maxIter}" ]; then str="${str}, "; fi
    index=$(( index + 1 ))
    
    if [ "$k" == "$maxIter" ]; then 
        str="${str} ti \"CSP\" " ; 
    else
        str="${str} ti \"\", "; 
    fi

    lw=$(( lw + 1 ))

done
gnustr="${gnustr}; ${str}"

# # Plot angles differences
# if [ "$use_max_radius" == "1" ]; then
#     max_radius=`cat .max_radius`
# else
#     max_radius=1
# fi
# xtics=`echo  "0.5 * $max_radius" | bc -l`
# yposition=`echo  "0.9 * $max_radius" | bc -l`
# max_radius=`echo  "1.1 * $max_radius" | bc -l`
# title="Angles difference"
# str="set title '${title}'"
# str="${str}; set key title ''"
# str="${str}; unset border"
# str="${str}; unset logscale y"
# str="${str}; set angles degrees"
# str="${str}; set grid polar 15"
# str="${str}; set polar"
# str="${str}; set size square"
# str="${str}; set style data points"
# str="${str}; set xtics ${xtics} norotate"
# str="${str}; unset ytics"
# str="${str}; set format x \"%2.0te%L\""
# str="${str}; set yrange [-${max_radius}:${max_radius}]"
# str="${str}; set xrange [-${max_radius}:${max_radius}]"
# str="${str}; set label 'radius=sin(theta) vs. angle=phi' at -${max_radius},-${yposition}"
# str="${str}; plot "
# # doIterations=$(seq 2 1 ${maxIter})
# doIterations="2 ${maxIter}"
# for k in ${doIterations}
# do 
#     titer=`printf %02d $k`
#     file=${expName}_DATA_POLARANGDIF_${titer}
#     str="${str} '${file}' using 7:6 pt 2 ps .7 ti \"USP iter $k\", "
#     str="${str} '${file}' using 4:3 pt 4 ps .7 ti \"CSP iter $k\" "
#     if [ "${k}" -lt "$maxIter" ]; then str="${str}, "; fi
# done
# gnustr="${gnustr}; ${str}"

# Run gnuplot
command=`echo "${gnustr}" | gnuplot`
# echo ${command}
$command 

# Obtain the axial, sagital and coronal cuts for all the final density
# map. Then concatenate them in one image.
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
            modellpf=${expName}_${XSP}_${iter}_lpf.mrc
			# proc3d $model $modellpf apix=1.4 lp=10
            echo -ne "Generating model cuts for ${model}\r"
            axial=${expName}_${XSP}_${iter}_axial
            coronal=${expName}_${XSP}_${iter}_coronal
            sagital=${expName}_${XSP}_${iter}_sagital
            
            newstack ${model} ${axial}.mrc -secs 124 -scale 0,255 -mode 0 > /dev/null
            mrc2tif ${axial}.mrc ${axial}.tif > /dev/null
            convert -normalize -crop 128x128+64+64 ${axial}.tif ${axial}.png
            
            rotatevol -input ${model} -output ${sagital}.mrc -angles 0,90,0 > /dev/null
            newstack ${sagital}.mrc ${sagital}.mrc -secs 128 -scale 0,255 -mode 0 > /dev/null
            mrc2tif ${sagital}.mrc ${sagital}.tif > /dev/null
            convert -normalize -rotate 90 -crop 128x128+64+64 -negate ${sagital}.tif ${sagital}.png
            convert -negate ${sagital}.png ${sagital}.png
            
            rotatevol -input ${model} -output ${coronal}.mrc -angles 0,90,90 > /dev/null
            newstack ${coronal}.mrc ${coronal}.mrc -secs 128 -scale 0,255 -mode 0 > /dev/null
            mrc2tif ${coronal}.mrc ${coronal}.tif > /dev/null
            convert -normalize -rotate 90 -crop 128x128+64+64 -negate ${coronal}.tif ${coronal}.png
            convert -negate ${coronal}.png ${coronal}.png
            
            rm -f *.mrc~ $modellpf
            montage -tile 3x -geometry 256x256 ${expName}_${XSP}_${iter}_{axial,coronal,sagital}.png ${XSP}_${iter}.png
        done
        montage -tile 1x -geometry 768x256 ${XSP}_*.png ${XSP}.png
        rm -f ${XSP}_*.png *_{axial,coronal,sagital}.*
    done
    montage -tile 2x1 -geometry 768x$(( 256 * ( maxIter - 1 ) )) CSP.png USP.png ${resultsdir}/${expName}_model_cuts.png
    rm -f CSP.png USP.png
    echo ""
    tput cnorm
fi

# Move the results to the results folder
#mv ${outputFile} ${resultsdir}
echo "Results plotted in ${outputFile}"

# Show the results
if [ "${CLUSTER_ID}" == "NONE" ]
then 
    if [ "${outputtype}" == "ps" ]
    then
        evince ${outputFile} &
    else  
        display ${outputFile} &
    fi
fi
