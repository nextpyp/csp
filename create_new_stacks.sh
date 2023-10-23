#!/bin/bash
#
# Version: $Id: create_new_stacks.sh 511 2010-09-16 23:14:02Z fefo $
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

dotParxFile=`grep ^ParFile ${pfile} | awk '{print $2}'`_01.parx #frealign/GroEL_02_noise_00_1.par
filesBasename=`grep ^ParFile ${pfile} | awk '{print $2}'`

binning=`grep ^BinningForRefinement ${pfile} | awk '{print $2}'`

particleColumn=\$14
micrographColumn=\$17
tiltSeriesColumn=\$8

cd frealign

numberOfTiltSeries=`cat ${dotParxFile} | grep -v C | awk '{print '${tiltSeriesColumn}'}' | sort -nr | head -1`
let numberOfTiltSeries=numberOfTiltSeries+1

echo "${numberOfTiltSeries} tilt series to process:"

# In order to cut the slices from the stack newstack receives the
# sections (-secs) numbered from 0 to N-1, being N the total number of
# slices in the stack. For the second (and higher) tilt series the
# numbers of images have added the total number of images of the
# previous tilt series. It is necessary to keep record of this number
# in order to make the correct "restacking". See the codeline where
# the number of the lines are obtained and prevProjection is
# subtracted.
prevProjection=0

# clear state
rm -f ${filesBasename}*.done

total=0

rm -rf data

# Create a folder where to store all the files
mkdir -p data

for (( tvalue=0; tvalue<${numberOfTiltSeries}; tvalue++ ))
do
    let tvalueindex=tvalue+1
    
    currentseries=`cat ${filesBasename}.series | head -$tvalueindex | tail -1 | awk '{print $2}'`
    currentparfile=${currentseries}_01.par

    stackFile=${currentseries}_stack.mrc

    swarmfile=${currentseries}_create_new_stacks.swarm
    rm -f swarm/$swarmfile

    if ! [ -f ${stackFile} ]; then
        echo -ne "[ERROR] `basename $0`: File ${stackFile} does not exist."
        exit 1
    fi
    echo -n "Processing $currentseries : "

    numberOfParticles=`cat ${currentparfile} | grep -v C | awk '{print '${particleColumn}'}' | sort -nr | head -1`
    numberOfParticles=`echo "scale=0; ($numberOfParticles+1)/1" | bc`

    echo -n "$numberOfParticles particle(s), "
	
    # For all the particles
    column=$particleColumn
    for (( value=0; value<${numberOfParticles}; value++ ))
    do
        # echo -n "Processing Particle `printf %04d $value`:          "
        basename=${filesBasename}_T`printf %02d $tvalue`_P`printf %04d $value`
        parFilename=${basename}_01.par
        mrcFilename=${basename}_01.mrc
        stackFilename=${basename}_stack.mrc
        binnedStackFilename=${basename}_stack_bin`printf %02d $binning`.mrc
        # Get the number of the lines
        lines=`cat ${dotParxFile} | grep -v C | awk '{if ( '${column}'=='${value}' && $8=='${tvalue}' ) print $0}' | awk '{printf (( $1 - 1 - '${prevProjection}' ))","}'`
		
        if [ "$lines" != "" ]
        then
            if [ "$binning" -gt "1" ]
            then
                    gstack="cd ..; newstack -input ${stackFile} -output data/${stackFilename} -secs ${lines}; binvol data/${stackFilename} -xbinning $binning -ybinning $binning -zbinning 1 data/${binnedStackFilename}; touch ${stackFilename}.done; cd - > /dev/null"
            else
                    gstack="cd ..; newstack -input ${stackFile} -output data/${stackFilename} -secs ${lines}; touch ${stackFilename}.done; cd - > /dev/null"
            fi
            echo $gstack >> swarm/$swarmfile
        fi
		
        continue

        # Write the parx file
        if ! [ -f data/${parFilename} ]; then
            cat ${dotParxFile} | grep -v C | awk '{if( '${column}'=='${value}' && $8=='${tvalue}' ) print $0}' > data/${parFilename}
            ${CSPDIR}/reformat_parx_file.sh data/${parFilename} > data/${parFilename}_temp
            mv data/${parFilename}_temp data/${parFilename}
        fi
    done

    # For all the micrographs in tilt series
    numberOfMicrographs=`cat ${currentparfile} | grep -v C | awk '{print '${micrographColumn}'}' | sort -nr | head -1`
    let numberOfMicrographs=numberOfMicrographs+1
    echo -n "$numberOfMicrographs micrograph(s)"; echo

    column=$micrographColumn
    for (( value=0; value<${numberOfMicrographs}; value++ ))
    do
        # echo -n "Processing Micrograph `printf %04d $value`:        "
        basename=${filesBasename}_T`printf %02d $tvalue`_M`printf %04d $value`
        parFilename=${basename}_01.par
        mrcFilename=${basename}_01.mrc
        stackFilename=${basename}_stack.mrc
        binnedStackFilename=${basename}_stack_bin`printf %02d $binning`.mrc
        # Get the number of the lines
        lines=`cat ${dotParxFile} | grep -v C |awk '{if ( '${column}'=='${value}' && $8=='${tvalue}' ) print $0}' | awk '{printf (( $1 - 1 - '${prevProjection}' ))","}'`

        if [ "$lines" == "" ]
        then
            continue
        fi

        cat ${dotParxFile} | grep -v C | awk '{if ( '${column}'=='${value}' && $8=='${tvalue}' ) print $0}' | awk '{printf (( $1 - 1 - '${prevProjection}' ))"'\\n'"}' > $basename
        split -l 100 -a 5 $basename $basename; rm -f $basename
        
        list=
        counter=0
        command=
        for i in `ls ${basename}?????`
        do
            indexes=
            for j in `cat $i`
            do
                indexes=`echo ${indexes}${j},`
            done

            # newstack parameter file
            command="${command} newstack -input ${stackFile} -output ${basename}_${counter}.mrc -secs ${indexes}; "
            list=`echo $list ${basename}_${counter}.mrc`
            #echo $com; echo; $com
            let counter=counter+1
        done

        # merge all sub-stacks
        command="${command} newstack $list data/${stackFilename}; "

        command="${command} rm -f ${basename}????? $list; "

        if [ "$binning" -gt "1" ]
        then
                gstack="cd ..; $command binvol data/${stackFilename} -xbinning $binning -ybinning $binning -zbinning 1 data/${binnedStackFilename}; touch ${stackFilename}.done; cd - > /dev/null"
        else
                gstack="cd ..; $command touch ${stackFilename}.done; cd - > /dev/null"
        fi
        echo $gstack >> swarm/$swarmfile

        continue
		
        # Write the parx file
        if ! [ -f data/${parFilename} ]
        then
            cat ${dotParxFile} | grep -v C | awk '{if( '${column}'=='${value}' && $8=='${tvalue}' ) print $0}' > data/${parFilename}
            ${CSPDIR}/reformat_parx_file.sh data/${parFilename} > data/${parFilename}_temp
            mv data/${parFilename}_temp data/${parFilename}
        fi
    done
    
    let total=total+`cat swarm/$swarmfile | wc -l | awk '{print $1}'`
	
    temp=`cat ${dotParxFile} | grep -v C | awk '{if ( $8=='${tvalue}' ) print $0}' | wc -l`
    echo Tilt series $tvalue has $temp particle projections.
    prevProjection=$(( prevProjection + temp ))
    
    # Clean up
    rm -f data/*~
    echo "[Done]"

    # submit all newstack calls to swarm
    cd swarm
    com="swarm -f $swarmfile ${AVAILABLE_NODES} -b 20"
    echo $com; $com
    cd - > /dev/null
done

echo "Waiting for $total processes to finish..."
processed=0
while [ $processed != $total ]
do
	sleep 2
	processed=`find . -type f -name "${filesBasename}*.done" -exec ls '{}' \; | wc -l | awk '{print $1}'`
	echo $processed of $total processed
done
echo All processes finished. 

find . -type f -name "${filesBasename}*.done" -exec rm '{}' \;

# merge all images in single mrc stack
totalstacks=`wc -l ${filesBasename}.series | awk '{print $1}'`
echo $totalstacks > inputfile
for i in `cat ${filesBasename}.series | awk '{print $2}'`
do
    echo ${i}_stack.mrc >> inputfile
    sections=`cat ${i}_01.par | grep -v C | wc -l | awk '{print $1}'`
    let secs=sections-1
    echo "0-"$secs >> inputfile
done

rm -f ${filesBasename}_stack.mrc

newstack -fileinlist inputfile -output ${filesBasename}_stack.mrc

if [ -f ${filesBasename}_stack.mrc ]
then
    rm inputfile
    for i in `cat ${filesBasename}.series | awk '{print $2}'`
    do
        rm ${i}_stack.mrc ${i}_01.par
    done
    # overwrite .series file
    echo -e "0 \t ${filesBasename}" > ${filesBasename}.series
fi

cd - > /dev/null