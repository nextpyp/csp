#!/bin/bash

# SPA_DIR is the folder where the frealign/<script.sh> are (this one `pwd`)
spaDir=`grep WorkingFolder parameters.config | awk '{print $2}'`
export SPA_DIR=${spaDir}

iterations=`grep NumberOfIterations parameters.config | awk '{print $2}'`

fparams=frealign_parameters
cd frealign

for kfile in `ls ${fparams}_normal_noise_*`
do
    echo "Running USP on ${kfile}"
	for (( i = 2; i <= $iterations; i++ ))
    # cp ${kfile} ${fparams}
	do
		cp ${kfile} ${fparams}_`printf %02d ${iterations}`
	done
    
    experimentName=`grep data_input ${fparams} | awk '{print $2}'`
    
    # Create experiment folder.
    folder=${experimentName}_USP
    rm -rf ${folder}
    mkdir ${folder}

    # Link the data files needed.
    lnBasename=${experimentName}
    lnParFilename=${lnBasename}_1.par
    ln -s -f ../data/${experimentName}.parx ${lnParFilename}
    lnMrcFilename=${lnBasename}_1.mrc
    ln -s -f ../data/${lnMrcFilename} .
    lnStackFilename=${lnBasename}_stack.mrc
    ln -s -f ../data/${lnStackFilename} .

    log_file=${folder}/${experimentName}.log
    echo "# Log file if frealign/${log_file}"

    # Run Frealign
	# ./frealign_iter_nolog.sh 2 100 > $log_file
    frealign_dispatcher.sh ${iterations} 100 > $log_file
    
    # Copy output files to experiment folder.
    cp -u {scratch,maps}/${experimentName}_{1.,2.,2_}* ${folder}
    cp -u log/${experimentName}_2_* ${folder}

    # Remove links created.
    rm -f ${experimentName}*{_1.par,_1.mrc,_stack.mrc}
    
    # Backup the experiment
    zipFilename=../${experimentName}_USP_`date +%Y-%m-%d-%H%M%S`.zip
    echo -n "=> Back up files:                "
    zip -q ${zipFilename} \
        $fparams \
        $log_file \
        frealign_parameters \
        microscope_parameters \
        ${folder}/*
    echo "[Done]"
    echo "   File ${zipFilename} created."
    
done
cd - > /dev/null
