#!/bin/bash
#
# Version: $Id: create_new_stacks.sh 468 2010-08-06 17:54:50Z fefo $
#

function usage {
    echo "Usage: ./`basename $0` -G groundtruth_parX_file -U Unconstrained_SP_par_file -C Constrained_SP_par_file -o log_file [-n particle_projection_number]"
}

write_costs=0
linemode=0
while [ "$1" != "" ]
do
    case "$1" in
        -n ) 
            shift
            kline=$1
            linemode=1
            ;;
        -G) 
            shift
            iniParFile=$1
            ;;
        -U )
            shift
            USPParFile=$1
            ;;
        -C )
            shift
            CSPParFile=$1
            ;;
        -o )           
            shift
            log_file=$1
            ;;
        --output-log )           
            shift
            output_log=$1
            write_costs=1
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

# Check existence of the input files 
if ! [ -f ${iniParFile} ]; then
    echo "[ERROR] `basename $0`: File ${iniParFile} does not exist."
    exit 1
fi
if ! [ -f ${USPParFile} ]; then
    echo "[ERROR] `basename $0`: File ${USPParFile} does not exist."
    exit 1
fi
if ! [ -f ${CSPParFile} ]; then
    echo "[ERROR] `basename $0`: File ${CSPParFile} does not exist."
    exit 1
fi
if [ "$write_costs" == "1" ]
then
    if ! [ -f ${output_log} ]; then
        echo "[ERROR] `basename $0`: File ${output_log} does not exist."
        exit 1
    fi
fi
rm -f ${log_file}

particleColumn=\$14
micrographColumn=\$17
tiltSeriesColumn=\$8
numberOfParticles=`awk '{print '${particleColumn}'}' ${iniParFile} | sort -nr | head -1`
numberOfMicrographs=`awk '{print '${micrographColumn}'}' ${iniParFile} | sort -nr | head -1`
numberOfTiltSeries=`awk '{print '${tiltSeriesColumn}'}' ${iniParFile} | sort -nr | head -1`

if [ "$linemode" == "1" ]; then  
    awk '{if($1 == '$kline') printf "GTH@%d:\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", $1, $2, $3, $4, $5, $6 }' ${iniParFile}
    awk '{if($1 == '$kline') printf "USP@%d:\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", $1, $2, $3, $4, $5, $6 }' ${USPParFile}
    awk '{if($1 == '$kline') printf "CSP@%d:\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", $1, $2, $3, $4, $5, $6 }' ${CSPParFile}
else
    #echo -e "part_proj_index Groundtruth(psi theta phi shx shy) USP(psi theta phi shx shy) CSP(psi theta phi shx shy)"
    column=$micrographColumn
    for (( value=1; value<=${numberOfMicrographs}; value++ ))
    do
        lines=`awk '{if('${column}'=='${value}') print $0}' ${iniParFile} | awk '{print $1}'`
        # echo "$value "
        for k in $lines
        do
            # awk '{if($1 == '$k') printf "%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t\t", $1, $2, $3, $4, $5, $6 }' ${iniParFile}
            # awk '{if($1 == '$k') printf "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t\t", $2, $3, $4, $5, $6 }' ${USPParFile}
            # awk '{if($1 == '$k') printf "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t", $2, $3, $4, $5, $6 }' ${CSPParFile}
            awk '{if($1 == '$k') printf "%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t\t", $1, $2, $3, $4, $5, $6 }' ${iniParFile} >> ${log_file}
            awk '{if($1 == '$k') printf "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t\t", $2, $3, $4, $5, $6 }' ${USPParFile} >> ${log_file}
            awk '{if($1 == '$k') printf "%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t", $2, $3, $4, $5, $6 }' ${CSPParFile} >> ${log_file}
            if [ "$write_costs" == "1" ]; then
                costos=`grep "Refined .* ( 0 $(printf "%02d" $value)" ${output_log} | awk '{print $17}'`
                # echo ${costos} | awk '{printf "%f\t%f\n", $1, $2}'
                echo ${costos} | awk '{printf "%f\t%f\n", $1, $2}' >> ${log_file}
            else
                echo -e "0\t0"  >> ${log_file}
            fi
        done
    done

    # Run SciLab for computing metrics and results processing
    echo "Run SciLab"
    sed -e "s|REPLACE_IN_FILE|${log_file}|" scilab/CompareResults.sce.base > scilab/CompareResults.sce
    sed -i -e "s|REPLACE_OUT_FILE|${log_file}.scilab|" scilab/CompareResults.sce
    scilab -nw -nb -f "scilab/CompareResults.sce" 
fi

