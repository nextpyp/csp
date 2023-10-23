#!/bin/bash
#
# Process the results of an experiment, the zip file or the log file.
#
# Version: $Id: process_csp_experiment.sh 478 2010-08-13 22:47:37Z fefo $
#

function usage {
    echo "Usage: `basename $0` [-z zip_file | -l log_file ] [-h] [-s]"
}

show=0
stopOnErrors=1
outFile=
while [ "$1" != "" ]
do
    case "$1" in
        -z | --zip ) 
            shift
            zipFile=$1
            ;;
        -l | --log) 
            shift
            logFile=$1
            ;;
        -o | --output )
            shift
            outFile=$1
            ;;
        -s | --show )
            show=1
            ;;
        -nse | --no-stop-on-errors )
            stopOnErrors=0
            ;;
        -h | --help )           
            usage
            exit
            ;;
        * )                     
            usage
            exit 1
    esac
    shift
done

# Unzip the log file and move to the temp folder
if [ -e "${zipFile}" ]
then
    folder=`dirname $zipFile` #temp #_`date +%Y%m%d%H%M%S`
    logFile=`unzip -q -l $zipFile | grep output.log | awk '{print $4}'`
    unzip -q -u $zipFile -d $folder $logFile
    cd $folder
fi

if [ "${outFile}" == "" ]
then
    outFile=${logFile}_refined
fi
rm -f $outFile

# Check if there are any errors in the log file
nErrors=`grep ERROR ${logFile} | wc -l`
if [ $nErrors -gt 0 ]
then
    echo "There are errors (${nErrors}) in the log file ${logFile}"
    if [ "${stopOnErrors}" == "1" ]
    then
        exit 1
    fi
fi

refinementFile=refinedFile
grep Refin ${logFile} > $refinementFile
exec 0<"${refinementFile}"
while read -r line
do
    word=`echo $line | grep Refining | awk '{print $2}'`
    if [[ "$word" == "Refining" ]]
    then
        object=`echo $line | awk '{print $3}'`
        number=`echo $line | awk '{print $4}'`
        #echo "${object} (${number})"
        read -r line
        line1=`echo $line | tr , '\ '`
        tarray1=( $line1 )
        read -r line
        line2=`echo $line | tr , '\ '`
        tarray2=( $line2 )
        #result=`echo ${tarray1[3]} - ${tarray2[4]} | bc -l`
        #echo $result
        if [ "${object}" == "micrograph" ]; then echo -ne "0\t${number}\t" >> $outFile; fi
        if [ "${object}" == "particle" ]; then echo -ne "1\t${number}\t" >> $outFile; fi
        for ind in 9 10 11 12 13 16
        do
            echo -ne "`printf %7.3f ${tarray1[$ind]}`\t" >> $outFile
        done
        for ind in 9 10 11 12 13 16
        do
            echo -ne "`printf %7.3f ${tarray2[$ind]}`\t" >> $outFile
        done
        echo "" >> $outFile
    fi
done
rm -f $refinementFile
# if [ "${show}" == "1" ]
# then
#     cat $outFile
# fi

cd - > /dev/null

# Run SciLab for computing metrics and results processing
echo "Run SciLab"
sed -e "s|REPLACE_FILE|${folder}/${outFile}|" scilab/ProcessExperiment.sce.base > scilab/ProcessExperiment.sce
scilabLogFile=${folder}/${outFile}_scilab.log
scilab -nw -nb -f "scilab/ProcessExperiment.sce" > ${scilabLogFile}
if [ "${show}" == "1" ]
then
    cat ${scilabLogFile}
fi
