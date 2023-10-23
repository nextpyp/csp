#!/bin/bash
#
# Version: $Id: run_csp_experiment.sh 498 2010-08-21 23:17:37Z fefo $
#

# The configuration file is always the same
parameters_file=parameters.config

GTparFile=$1

# # Use a log file?
# if [ -n "$1" ]; then
    # log_file=$1
# fi

# Grep the basename from the <ParFile> in the parameters file. It is
# the second argument in the line starting with 'ParFile' without '.
# experimentName=`grep ^ParFile ${parameters_file} | awk '{print $2}' | sed -e "s/\(.*\)_1.par.*/\1/"`
experimentName=`grep ^ParFile ${parameters_file} | awk '{print $2}'`
iters=`grep ^NumberOfIterations ${parameters_file} | awk '{print $2}'`

# build FREALIGN directory structure
mkdir -p frealign
mkdir -p frealign/log
mkdir -p frealign/scratch
mkdir -p frealign/swarm
mkdir -p frealign/maps

echo -e "1\t$experimentName" > frealign/${experimentName}.series

# Build the new stacks, one per micrographs and one per particle.
${CSPDIR}/create_new_stacks.sh

# # Create links in the Frealign folder to the stack, images and parx
# # files.
# for kfile in data/${experimentName}/*_1.parx
# do 
    # lnBasename=`echo $kfile | sed -e "s|data/${experimentName}/\(.*\)_1.parx|\1|"`
    # lnParFilename=${lnBasename}_1.parx
    # ln -s -f ../data/${experimentName}/${lnParFilename} frealign/
    # lnMrcFilename=${lnBasename}_1.mrc
    # ln -s -f ../data/${experimentName}_1.mrc frealign/${lnMrcFilename}
    # lnStackFilename=${lnBasename}_stack.mrc
    # ln -s -f ../data/${experimentName}/${lnStackFilename} frealign/
# done

# # Create the output folder in the Frealign folder structure
# folderName=frealign/${experimentName}
# if ! [ -d ${folderName} ]; then
    # mkdir -p ${folderName}
    # echo "=> Folder ${folderName} created."
# else
    # rm -rf ${folderName}
    # mkdir -p ${folderName}
    # echo "=> Folder ${folderName} emptied."
# fi

cp frealign/${experimentName}_01.mrc frealign/maps
ln -s `pwd`/frealign/${experimentName}_01.mrc frealign/maps/${experimentName}_CSP_01.mrc
cp frealign/${experimentName}_01.parx frealign/maps

# # Clean the log files
# echo -n "=> Clean log files:              "
# if [ -n "${log_file}" ]; then
    # rm -rf ${log_file}
# fi
rm -rf run_frealign_reconstruction.log
rm -rf run_frealign_refine.log
echo "[Done]"

# Run the CSP framework
# echo "=> Running (go for a cup of tea, or the whole kettle...)"
# time ${HOME}/scripts/SPA/csp/bin/csp parameters.config

rm -f csp.done

# launch CSP run

# submit and get job number
command="qsub -v np=40,mydir=`pwd`,myjob=${experimentName}_CSP -l nodes=10:o2600:dc -e `pwd` ${CSPDIR}/mpi_csp.sh"
echo $command
# jobnumber=`$command | sed -e 's/.biobos//'`

echo Submitted qsub job ${jobnumber}.biobos

# Launch USP run
# echo "cd frealign; frealign_dispatcher.sh $iters 100 > ../${experimentName}_USP.log; cd - > /dev/null" > ${experimentName}_USP.swarm
# swarm -f ${experimentName}_USP.swarm

# echo WARNING - NOT RUNNING FREALIGN!!
# exit 1

let iters=iters+1
cd frealign; frealign_dispatcher.sh $iters 100; cd - > /dev/null

# change all frealign filenames to _USP
cd frealign/maps
for (( i = 2; i <= $iters; i++ ))
do
	mv ${experimentName}_`printf %02d $i`.par ${experimentName}_USP_`printf %02d $i`.par
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

# $CSPDIR/compare_experiments.sh ${experimentName}

# echo Wait until job completed...
# while ! [ -f csp.done ]
# do
	# sleep 3
# done

# rm -f csp.done

# echo Job $jobnumber.biobos done.

# # save log file in current directory
# mv ~/csp.o$jobnumber ${experimentName}.o
# mv ~/csp.e$jobnumber ${experimentName}.e

# #Remove links created
# rm -f frealign/${experimentName}*{_1.parx,_1.mrc,_stack.mrc}
# rm -f frealign/${experimentName}_T*.par

# # Backup the experiment
# zipFilename=${experimentName}_`date +%Y-%m-%d-%H%M%S`.zip
# echo -n "=> Back up files:                "
# zip -q ${zipFilename} \
    # $parameters_file \
    # $log_file \
    # run_frealign_refine.log \
    # run_frealign_reconstruction.log \
    # frealign/frealign_parameters \
    # frealign/microscope_parameters \
    # ${folderName}
# echo "[Done]"
# echo "   File ${zipFilename} created."
