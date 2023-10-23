#!/bin/bash
#
# Version: $Id: create_new_stacks_single.sh 511 2010-09-16 23:14:02Z fefo $
#

if ! [ -n "$2" ]
then
  echo "Usage: `basename $0` tvalue prevProjection"
  exit 1
else
  echo `basename $0` $@; echo
fi  

# go to frealign folder
cd ..

pfile=../parameters.config
if ! [ -f ${pfile} ]; then
    echo "[ERROR] `basename $0`: File ${pfile} does not exist."
    exit 1
fi

dotParxFile=`grep ^ParFile ${pfile} | awk '{print $2}'`_01.parx
filesBasename=`grep ^ParFile ${pfile} | awk '{print $2}'`

particleColumn=\$14
micrographColumn=\$17

tvalue=$1
prevProjection=$2

let tvalueindex=tvalue+1

currentseries=`cat ${filesBasename}.series | head -$tvalueindex | tail -1 | awk '{print $2}'`
currentparfile=${currentseries}_01.par

stackFile=${currentseries}_stack.mrc

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
    stackFilename=${basename}_stack.mrc
    # Get the number of the lines
    lines=`cat ${dotParxFile} | grep -v C | awk '{if ( '${column}'=='${value}' && $8=='${tvalue}' ) print $0}' | awk '{printf (( $1 - 1 - '${prevProjection}' ))","}'`
            
    if [ "$lines" != "" ]
    then
        com="${IMOD_DIR}/bin/newstack -input ${stackFile} -output data/${stackFilename} -secs ${lines}"
        echo $com; $com
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
    stackFilename=${basename}_stack.mrc
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
    #command=
    for i in `ls ${basename}?????`
    do
        indexes=
        for j in `cat $i`
        do
            indexes=`echo ${indexes}${j},`
        done

        # newstack parameter file
        # command="${command} ${IMOD_DIR}/bin/newstack -input ${stackFile} -output ${basename}_${counter}.mrc -secs ${indexes}; "
        com="${IMOD_DIR}/bin/newstack -input ${stackFile} -output ${basename}_${counter}.mrc -secs ${indexes}"
        list=`echo $list ${basename}_${counter}.mrc`
        echo $com; echo; $com
        let counter=counter+1
    done

    # merge all sub-stacks
    if [ $counter -eq 1 ]
    then
        #command="${command} mv $list data/${stackFilename}; "
        com="mv $list data/${stackFilename}"
    else
        #command="${command} ${IMOD_DIR}/bin/newstack $list data/${stackFilename}; "
        com="${IMOD_DIR}/bin/newstack $list data/${stackFilename}"
    fi
    echo $com; echo; $com

    #command="${command} rm -f ${basename}????? $list"
    #echo $command; $command
    com="rm -f ${basename}????? $list"
    echo $com; echo; $com
done

touch ${currentseries}.done
