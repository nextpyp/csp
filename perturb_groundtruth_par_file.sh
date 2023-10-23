#!/bin/bash
#
# Version: $Id: perturb_groundtruth_par_file.sh 511 2010-09-16 23:14:02Z fefo $
#

function usage {
    echo ""
    echo "NAME"
    echo -e "\t`basename $0`\n"
    echo -e "SYNOPSIS"
    echo -e "\t./`basename $0` -i input_par_file -m A B C D E -o output_par_file\n"
    echo "DESCRIPTION"
    echo -e "\tGenerates a perturbed version of the parx input file. The perturbation"
    echo -e "\tcould be in any of the CSP parameters (tilt angle, tilt axis angle,"
    echo -e "\tpsi, theta, phi) randomly perturbed by a Gaussian noise of standard"
    echo -e "\tdeviation controlled by the binary mask "A B C D E".\n"
    echo -e "\tDepends on the function <randn>. ./randn mean sigma seed This function"
    echo -e "\tcould be compiled from the gaussian_distribution.cpp file.\n"
}


if  [ ${#@} == 0 ]; then usage; exit; fi
 
save_file=0
prta=1
prtxa=1
prn0=1
prn1=1
prn2=1
while test -n "$1"
do
    case "$1" in
        -i )
            shift
            iniParFile=$1
            ;;
        -o )
            shift
            outParFile=$1
            save_file=1
            ;;
        -m )
            shift
            prta=$1 
            shift
            prtxa=$1
            shift
            prn0=$1
            shift
            prn1=$1
            shift
            prn2=$1
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
if [ -f ${outParFile} ]; then
    rm -f ${outParFile}
fi

cat ${iniParFile} | grep C > ${outParFile}

# Turn off cursor
tput civis

particleColumn=\$14
micrographColumn=\$17
numberOfParticles=`cat ${iniParFile} | grep -v C | awk '{print '${particleColumn}'}' | sort -nr | head -1`
let numberOfParticles=numberOfParticles+1
numberOfMicrographs=`cat ${iniParFile} | grep -v C |awk '{print '${micrographColumn}'}' | sort -nr | head -1`
let numberOfMicrographs=numberOfMicrographs+1

echo Processing $numberOfParticles particles and $numberOfMicrographs micrographs

# Draw the random angles to add to the micrographs' angles
for (( k=1; k<=$numberOfMicrographs; k++ ))
do
    ind1=$(( $k * 2 - 1 ))
    ind2=$(( $k * 2 ))
    # rmgp[$ind1]=$(${CSPDIR}/bin/gaussian_distribution 0 $tilt_angle_noise_std `date +%N`) # $ind1
    # rmgp[$ind2]=$(${CSPDIR}/bin/gaussian_distribution 0 $noise_std `date +%N`) # $ind2
	rmgp[$ind1]=`echo "($RANDOM/32767-.5)*2*${prta}" | bc -l`
	rmgp[$ind2]=`echo "($RANDOM/32767-.5)*2*${prtxa}" | bc -l`
    # echo -e ${rmgp[$ind1]}"\t"${rmgp[$ind2]}
done

# Draw the random angles to add to the particles' euler angles
for (( k=1; k<=$numberOfParticles; k++ ))
do
    ind1=$(( $k * 3 - 2 ))
    ind2=$(( $k * 3 - 1 ))
    ind3=$(( $k * 3 ))
    # rprt[$ind1]=$(${CSPDIR}/bin/gaussian_distribution 0 $noise_std `date +%N`) # $ind1
    # rprt[$ind2]=$(${CSPDIR}/bin/gaussian_distribution 0 $noise_std `date +%N`) # $ind2
    # rprt[$ind3]=$(${CSPDIR}/bin/gaussian_distribution 0 $noise_std `date +%N`) # $ind3
	rprt[$ind1]=`echo "($RANDOM/32767-.5)*2*${prn0}" | bc -l`
	rprt[$ind2]=`echo "($RANDOM/32767-.5)*2*${prn1}" | bc -l`
	rprt[$ind3]=`echo "($RANDOM/32767-.5)*2*${prn2}" | bc -l`
    # echo -e ${rprt[$ind1]}"\t"${rprt[$ind2]}"\t"${rprt[$ind3]}
done

number_lines=`cat $iniParFile | wc -l `
number=1

# Write an extra header explaining the generation of this file.
# echo -ne "C\nC Perturbed version of ${iniParFile} with mask ${prta}${prtxa}${prn0}${prn1}${prn2} and sigma ${noise_std}.\nC `date +%Y-%m-%d-%H%M%S`\nC `echo $USER`\nC\n" >> ${outParFile}

exec 0<"${iniParFile}"
while read -r line
do
    # echo ${line:0:60}
    echo -ne "$number/$number_lines\r"
    number=$((number+1))
	if [ ${line:0:1} == 'C' ]
    then
        echo $line > /dev/null
	else
        
        # IFS=':'
        # echo "     ${line:0:500}"
        # IFS=' '

        tarray=( $line )
        img_index=${tarray[0]}
        psi=${tarray[1]}
        theta=${tarray[2]}
        phi=${tarray[3]}
        shift_x=${tarray[4]}
        shift_y=${tarray[5]}
        mag=${tarray[6]}
        film=${tarray[7]}
        df1=${tarray[8]}
        df2=${tarray[9]}
        angast=${tarray[10]}
        presa=${tarray[11]}
        dpres=${tarray[12]}
        ptl_index=${tarray[13]}
        tilt_angle=${tarray[14]}
        dose=${tarray[15]}
        scan_order=${tarray[16]}
        confidence=${tarray[17]}
        ptl_CCX=${tarray[18]}
        tilt_axis_angle=${tarray[19]}        
        normal_0=${tarray[20]}; normal_1=${tarray[21]}; normal_2=${tarray[22]};
        matrix_00=${tarray[23]}; matrix_01=${tarray[24]}; matrix_02=${tarray[25]}; matrix_03=${tarray[26]};
        matrix_10=${tarray[27]}; matrix_11=${tarray[28]}; matrix_12=${tarray[29]}; matrix_13=${tarray[30]};
        matrix_20=${tarray[31]}; matrix_21=${tarray[32]}; matrix_22=${tarray[33]}; matrix_23=${tarray[34]};
        matrix_30=${tarray[35]}; matrix_31=${tarray[36]}; matrix_32=${tarray[37]}; matrix_33=${tarray[38]}

        let micrograph=$scan_order+1
        let particle=$ptl_index+1
        indm1=$(( $micrograph * 2 - 1 ))
        indm2=$(( $micrograph * 2 ))
        indp1=$(( $particle * 3 - 2 ))
        indp2=$(( $particle * 3 - 1 ))
        indp3=$(( $particle * 3 ))

        # Generate random angles to add
        rta=0
        rtxa=0
        rn0=0
        rn1=0
        rn2=0
        rta=${rmgp[$indm1]}
        rtxa=${rmgp[$indm2]}
        rn0=${rprt[$indp1]}
        rn1=${rprt[$indp2]}
        rn2=${rprt[$indp3]}
        
        # echo "$img_index $micrograph $particle $rta $rtxa $rn0 $rn1 $rn2"
        
        # Compute new angles in CSP manifold.
        tilt_angle=`echo $tilt_angle + $rta | bc -l`
        tilt_axis_angle=`echo $tilt_axis_angle + $rtxa | bc -l`
        normal_0=`echo $normal_0 + $rn0 | bc -l`
        normal_1=`echo $normal_1 + $rn1 | bc -l`
        normal_2=`echo $normal_2 + $rn2 | bc -l`
                
        # Map and get the new angles in the USP manifold.
        mapped_angles=`${CSPDIR}/bin/SPA_euler_angles $tilt_angle $tilt_axis_angle 0 0 0 $normal_0 $normal_1 $normal_2 $matrix_00 $matrix_01 $matrix_02 $matrix_03 $matrix_10 $matrix_11 $matrix_12 $matrix_13 $matrix_20 $matrix_21 $matrix_22 $matrix_23 $matrix_30 $matrix_31 $matrix_32 $matrix_33 | tail -n 1`
        psi=`echo $mapped_angles | awk '{print $1}'`
        theta=`echo $mapped_angles | awk '{print $2}'`
        phi=`echo $mapped_angles | awk '{print $3}'`

        # Format and write the new .parx line
        ${CSPDIR}/frealign/format_dotpar_line.sh $img_index $psi $theta $phi $shift_x $shift_y $mag $film $df1 $df2 $angast $presa $dpres $ptl_index $tilt_angle $dose $scan_order $confidence $ptl_CCX $tilt_axis_angle $normal_0 $normal_1 $normal_2 $matrix_00 $matrix_01 $matrix_02 $matrix_03 $matrix_10 $matrix_11 $matrix_12 $matrix_13 $matrix_20 $matrix_21 $matrix_22 $matrix_23 $matrix_30 $matrix_31 $matrix_32 $matrix_33 >> ${outParFile}
        
	fi
done

# Turn on cursor
tput cnorm
    
