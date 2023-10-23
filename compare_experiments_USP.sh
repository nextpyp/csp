#!/bin/bash

dataset=$1
GTMapFile=$2
pixel_size=$3

maxIter=`ls ${dataset}_??.mrc | sed -e 's/'${dataset}'_\(.*\).mrc/\1/g' | bc -l | tail -1`
echo maxIter = $maxIter

outputFile=${dataset}.png
output=${dataset}_DATA.txt
rm -f $output
# echo -e "N\tRES\tRESGT\tRMEASURE" > $output

doIterations=$(seq 2 ${maxIter})
for n in ${doIterations}
do
	# align ground truth to current map and get FSC curve
	xsp_map=${dataset}_`printf %02d ${n}`.mrc
	if [ -f "${xsp_map}" ]; then
		echo Processing map $xsp_map
		# map is not neccesarily aligned to GT map, we need to align them before computing the FSC
		
		# align GT to current map
		~/code/ETTK/AlignRigid3D ${xsp_map} ${GTMapFile} 10 180 aligned.mrc
		# cp ${GTMapFile} aligned.mrc
		
		map=${dataset}_`printf %02d ${n}`
		fsc_file=${map}_fsc.txt
		fsc_gt_file=${map}_gt_fsc.txt
		proc3d aligned.mrc ${xsp_map} apix=${pixel_size} fsc=${fsc_gt_file} > /dev/null
		rm -f aligned.mrc

		# get FSC resolution 
		res=`${SPA_DIR}/general/fsc_cutoff.sh ${fsc_file} | awk '{print $5}'`
		resgt=`${SPA_DIR}/general/fsc_cutoff.sh ${fsc_gt_file} | awk '{print $5}'`
		rmeasure=`${SPA_DIR}/utils/rmeasure.sh ${map}.mrc | grep "Resolution at FSC = 0.5" | awk '{print $6}'`
		echo -e "$n\t$res\t$resgt\t$rmeasure" >> $output
	fi
done

width=1500
height=500
title=${dataset}
gnustr="set terminal png size ${width}, ${height}"
gnustr="${gnustr}; set output '${outputFile}'"
gnustr="${gnustr}; set multiplot layout 1, 3 title '${title} results'"

echo $gnustr

title="FSC plots"
xlabel="Resolution"
ylabel=""
str="set title '${title}'"
str="${str}; set xlabel '${xlabel}'"
str="${str}; set ylabel '${ylabel}'"
str="${str}; set autoscale y"
str="${str}; set xrange [0:.1]"
str="${str}; plot "
doIterations=$(seq 2 2 ${maxIter})
for k in ${doIterations}
do 
	let n=k+1
	str="${str} '${dataset}_`printf %02d $k`_fsc.txt' with lines lw 2 "
	if [ "$k" == "$maxIter" ]; then 
		str=${str}
	else
		str="${str} , "; 
	fi
done
gnustr="${gnustr}; ${str}"

title="FSC plots to GT"
xlabel="Resolution"
ylabel=""
str="set title '${title}'"
str="${str}; set xlabel '${xlabel}'"
str="${str}; set ylabel '${ylabel}'"
str="${str}; set autoscale y"
str="${str}; set xrange [0:.1]"
str="${str}; plot "
for k in ${doIterations}
do
	let n=k+1
	str="${str} '${dataset}_`printf %02d $k`_gt_fsc.txt' with lines lw 2 "
	if [ "$k" == "$maxIter" ]; then 
		str=${str}
	else
		str="${str} , "; 
	fi
done
gnustr="${gnustr}; ${str}"

str="set autoscale x; set yrange [0:40]; set xlabel 'Iterations'; set ylabel 'Resolution'; plot"
str="${str} '${output}' using 1:2 with linespoints ti \"Internal FSC\","
str="${str} '${output}' using 1:3 with linespoints ti \"FSC to GT\","
str="${str} '${output}' using 1:4 with linespoints ti \"RMEASURE\""
gnustr="${gnustr}; ${str}"

# Run gnuplot
#echo $gnustr > command
#command=`echo "${gnustr}" | gnuplot`
echo "${gnustr}" | gnuplot

# rm -f ${output} ${dataset}_??_gt_fsc.txt