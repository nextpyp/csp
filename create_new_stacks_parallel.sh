#!/bin/bash
#
# Version: $Id: create_new_stacks_parallel.sh 511 2010-09-16 23:14:02Z fefo $
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

dotParxFile=`grep ^ParFile ${pfile} | awk '{print $2}'`_01.parx
filesBasename=`grep ^ParFile ${pfile} | awk '{print $2}'`

tiltSeriesColumn=\$8

cd frealign

numberOfTiltSeries=`cat ${dotParxFile} | grep -v C | awk '{print '${tiltSeriesColumn}'}' | sort -nr | head -1`
let numberOfTiltSeries=numberOfTiltSeries+1

echo "${numberOfTiltSeries} tilt series to process:"

# clear state
rm -f *.done
rm -rf data

# Create a folder where to store all the files
mkdir -p data

swarmfile=create_new_stacks.swarm
rm -f swarm/$swarmfile

prevProjection=0

for (( tvalue=0; tvalue<${numberOfTiltSeries}; tvalue++ ))
do
    let tvalueindex=tvalue+1
    currentseries=`cat ${filesBasename}.series | head -$tvalueindex | tail -1 | awk '{print $2}'`
    echo "${CSPDIR}/create_new_stacks_single.sh $tvalue $prevProjection > ../../log/${currentseries}_create_new_stack.log" >> swarm/$swarmfile	
    temp=`cat ${dotParxFile} | grep -v C | awk '{if ( $8=='${tvalue}' ) print $0}' | wc -l`
    echo Tilt series $tvalue has $temp particle projections.
    prevProjection=$(( prevProjection + temp ))    
done

# submit all newstack calls to swarm
cd swarm
com="swarm -f $swarmfile ${AVAILABLE_NODES}"
echo $com; $com
cd - > /dev/null

echo "Waiting for ${numberOfTiltSeries} processes to finish..."
processed=0
while [ $processed != $numberOfTiltSeries ]
do
    sleep 2
    processed=`find . -type f -name "*.done" -exec ls '{}' \; | wc -l | awk '{print $1}'`
    echo $processed of $numberOfTiltSeries processed
done
echo All processes finished. 

find . -type f -name "*.done" -exec rm '{}' \;

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

com="newstack -fileinlist inputfile -output ${filesBasename}_stack.mrc"
echo $com; $com;

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