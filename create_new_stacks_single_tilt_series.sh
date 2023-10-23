#!/bin/bash
#
# Version: $Id: create_new_stacks.sh 468 2010-08-06 17:54:50Z fefo $
#

# Test whether command-line argument is present (non-empty).
if ! [ -n "$1" ]; then 
    pfile=parameters.config
else
    pfile=$1
fi  
if ! [ -f ${pfile} ]; then
    echo "[ERROR] `basename $0`: File ${pfile} does not exist."
    exit 1
fi

dotParxFile=`grep ^ParFile ${pfile} | awk '{print $2}'`_01.parx #frealign/GroEL_02_noise_00_01.par
filesBasename=`grep ^ParFile ${pfile} | awk '{print $2}'`

stackFile=${filesBasename}_stack.mrc
if ! [ -f frealign/${stackFile} ]; then
    echo -ne "[ERROR] `basename $0`: File ${stackFile} does not exist."
    exit 1
fi

cd frealign

particleColumn=\$14
micrographColumn=\$17
tiltSeriesColumn=\$8
numberOfParticles=`awk '{print '${particleColumn}'}' ${dotParxFile} | sort -nr | head -1`
numberOfMicrographs=`awk '{print '${micrographColumn}'}' ${dotParxFile} | sort -nr | head -1`
numberOfTiltSeries=`awk '{print '${tiltSeriesColumn}'}' ${dotParxFile} | sort -nr | head -1`

echo -n "=> Build new stacks (${numberOfTiltSeries},${numberOfParticles},${numberOfMicrographs}):    "

for (( tvalue=1; tvalue<=${numberOfTiltSeries}; tvalue++ ))
do
    # echo "Processing tilt Series `printf %02d $tvalue`"
    
    # Create a folder where to store all the files 
    mkdir -p data

    # For all the particle
    column=$particleColumn
    for (( value=1; value<=${numberOfParticles}; value++ ))
    do
        # echo -n "Processing Particle `printf %04d $value`:          "
        basename=${filesBasename}_T`printf %02d $tvalue`_P`printf %04d $value`
        parFilename=${basename}_01.parx
        mrcFilename=${basename}_01.mrc
        stackFilename=${basename}_stack.mrc
        # Get the number of the lines
        lines=`awk '{if('${column}'=='${value}') print $0}' ${dotParxFile} | awk '{printf (( $1 - 1 ))","}'`
        # Create the new stack
        if ! [ -f data/${stackFilename} ]; then
            newstack -input ${stackFile} -output data/${stackFilename} -secs ${lines} > /dev/null
            #ln -sf ../data/${stackFilename} frealign/${stackFilename}
        # else 
        #    echo -ne "${stackFilename} exists.\t"
        fi
        # Write the parx file
        if ! [ -f data/${parFilename} ]; then
            awk '{if('${column}'=='${value}') print $0}' ${dotParxFile} > data/${parFilename}
            ${CSPDIR}/reformat_parx_file.sh data/${parFilename} > data/${parFilename}_temp
            mv data/${parFilename}_temp data/${parFilename}
            #ln -sf ../data/${parFilename} frealign/${parFilename}
        # else
        #    echo -ne "${parFilename} exists.\t"
        fi
        ## Link the 3D model
        #ln -sf ../data/${filesBasename}_01.mrc frealign/${mrcFilename}
        # echo "[Done]"
    done

    # For all the micrographs
    column=$micrographColumn
    for (( value=1; value<=${numberOfMicrographs}; value++ ))
    do
        # echo -n "Processing Micrograph `printf %04d $value`:        "
        basename=${filesBasename}_T`printf %02d $tvalue`_M`printf %04d $value`
        parFilename=${basename}_01.parx
        mrcFilename=${basename}_01.mrc
        stackFilename=${basename}_stack.mrc
        # Get the number of the lines
        lines=`awk '{if('${column}'=='${value}') print $0}' ${dotParxFile} | awk '{printf (( $1 - 1 ))","}'`
        # Create the new stack
        if ! [ -f data/${stackFilename} ]; then
            newstack -input ${stackFile} -output data/${stackFilename} -secs ${lines} > /dev/null
            #ln -sf ../data/${stackFilename} frealign/${stackFilename}
        # else 
        #    echo -ne "${stackFilename} exists.\t"
        fi
        # Write the parx file
        if ! [ -f data/${parFilename} ]; then
            awk '{if('${column}'=='${value}') print $0}' ${dotParxFile} > data/${parFilename}
            ${CSPDIR}/reformat_parx_file.sh data/${parFilename} > data/${parFilename}_temp
            mv data/${parFilename}_temp data/${parFilename}
            #ln -sf ../data/${parFilename} frealign/${parFilename}
        # else
        #     echo -ne "${parFilename} exists.\t"
        fi
        ## Link the 3D model
        #ln -sf ../data/${filesBasename}_01.mrc frealign/${mrcFilename}
        # echo "[Done]"
    done
    
    # Clean up
    rm -f data/*~
    echo "[Done]"

done

cd - > /dev/null

exit 0
