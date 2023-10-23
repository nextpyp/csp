#!/bin/bash

parameters_file=parameters.config
experimentName=`grep ^ParFile ${parameters_file} | awk '{print $2}'`
iters=`grep ^NumberOfIterations ${parameters_file} | awk '{print $2}'`
let iters=iters+1

# run USP
cd frealign

rm -f mpirun.config

# generate starting model
mv frealign_parameters_02 frealign_parameters_02_backup
cp frealign_parameters_01 frealign_parameters_02

# make sure we are using the proper .par files
cp ${experimentName}_*_01.par maps
frealign_iter.sh 2 100

# store map
com="rm ${experimentName}_01.mrc"
echo $com; $com
com="cp maps/${experimentName}_02.mrc ${experimentName}_01.mrc"
echo $com; $com
com="cp maps/${experimentName}_02.mrc maps/${experimentName}_01.mrc"
echo $com; $com

# save PR and FSC plots
cp maps/${experimentName}_02_PR.png maps/${experimentName}_USP_01_PR.png
cp maps/${experimentName}_02_fsc.png maps/${experimentName}_USP_01_fsc.png
cp maps/${experimentName}_02_fsc.txt maps/${experimentName}_USP_01_fsc.txt

# restore parameters
com="mv frealign_parameters_02_backup frealign_parameters_02"
$com

# run usual frealign
frealign_dispatcher.sh $iters 100
cd - > /dev/null

# kill -SIGINT $$

# change all frealign filenames to _USP
cd frealign/maps
for (( i=2; i<=$iters; i++ ))
do
	string=${experimentName}_T??_`printf %02d $i`.par
	outpar=${experimentName}_USP_`printf %02d $i`.par
	rm -f $outpar
	first=1
	for k in `ls $string`
	do
		# mv ${experimentName}_`printf %02d $i`.par ${experimentName}_USP_`printf %02d $i`.par
		# mv $k `echo $k | sed -e "s|${experimentName}|${experimentName}_USP|"`
		if [ "$first" -eq 1 ]
		then
			cat $k >> $outpar
			first=0
		else
			cat $k | grep -v C >> $outpar
		fi
		# rm -f $k
	done
	# correct numbering
	/data/${USER}/frealign-workspace/frealign_scripts/frealign_fix_indexes.sh $outpar dummy
	mv dummy $outpar
	
	mv ${experimentName}_`printf %02d $i`_fsc.txt ${experimentName}_USP_`printf %02d $i`_fsc.txt
	mv ${experimentName}_`printf %02d $i`_fsc.png ${experimentName}_USP_`printf %02d $i`_fsc.png
	if [ $i -gt 2 ]
	then
		mv ${experimentName}_`printf %02d $i`_fsc_all.png ${experimentName}_USP_`printf %02d $i`_fsc_all.png
	fi
	mv ${experimentName}_`printf %02d $i`_PR.png ${experimentName}_USP_`printf %02d $i`_PR.png
	mv ${experimentName}_`printf %02d $i`.mrc ${experimentName}_USP_`printf %02d $i`.mrc
done
cd - > /dev/null

echo "Normal program termination"