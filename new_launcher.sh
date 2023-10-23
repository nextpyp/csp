#!/bin/bash

# check arguments
if ! [ -n "$2" ] 
then
  echo "Usage: `basename $0` experiment_name starting_model.mrc [with full path]"
  exit 1
else
  echo `basename $0` $@; echo
fi  

parameters_file=parameters.config
if ! [ -f ${parameters_file} ]
then
	echo ERROR - File ${parameters_file} not found.
	exit 1
fi

experiment=$1
starting_model=$2

uspswarmfile=${experiment}_USP.swarm
cspswarmfile=${experiment}_CSP.swarm
mpiswarmfile=${experiment}_CSP_MPI.swarm
spswarmfile=${experiment}.swarm
rm -f ${uspswarmfile} ${cspswarmfile} ${mpiswarmfile} ${spswarmfile}

#for (( km=0; km<=10; km+=1 ))
for (( km=7; km<9; km+=1 ))
do
	# create perturbed .parx file
	for i in `ls ../frealign/${experiment}_*_01.par`
	do
		# perturbed_file=`echo $i | sed -e "s|../frealign/||" | sed -e "s|_01.par||"`_ini_`printf %02d $km`_01.parx
		ts=`echo $i | sed -e "s|../frealign/${experiment}_||" | sed -e "s|_01.par||"`
		perturbed_file=${experiment}_ini_`printf %02d $km`_T`printf %02d $ts`_01.par
		if ! [ -f ${perturbed_file} ]
		then
			echo Generating file ${perturbed_file}
			# the perturbation in the tilt angle has to preserve the ordering of the images in the tilt series.
			# Therefore, if the perturbation is bigger than half the initial tilt separation, we keep limit the perturbation
			# to comply with the ordering constraint.
			kmf=`echo "3*$km" | bc -l`
			if [ `echo "$kmf<1" | bc -l` -eq 1 ]
			then
				kaf=$kmf
			else
				kaf=1.0
			fi
			com="${CSPDIR}/perturb_groundtruth_par_file.sh -i ${i} -m $kaf $kmf $kmf $kmf $kmf -o ${perturbed_file}"
			echo $com; $com
		else
			echo File ${perturbed_file} exists.
		fi
	done

	continue
	
	for (( nn=0; nn<=10; nn++ ))
	do
		# create stacks with noise
		for i in `ls ../frealign/${experiment}_*_stack.mrc`
		do
			# noisystack=`echo $i | sed -e "s|../frealign/||" | sed -e "s|_stack.mrc||"`_noise_`printf %02d $nn`_stack.mrc
			ts=`echo $i | sed -e "s|../frealign/${experiment}_||" | sed -e "s|_stack.mrc||"`
			noisystack=${experiment}_noise_`printf %02d $nn`_T`printf %02d $ts`_stack.mrc
			if ! [ -f $noisystack ]
			then
				if (( $nn > 0 ))
				then
					effectivenoise=`echo ".125*2^($nn-1)" | bc -l`
				else
					effectivenoise=$nn
				fi
				effectivenoise=`echo "2.0*5.0/4.0*$nn"| bc -l`
				effectivenoise=`echo "1.0*$nn"| bc -l`
				com="proc3d $i $noisystack noise=$effectivenoise"
				echo $com; $com
			fi
		done

		# create folder for current experiment
		current_experiment=${experiment}_ini_`printf %02d $km`_noise_`printf %02d $nn`
		mkdir -p ${current_experiment}

		cd ${current_experiment}

		# create frealign directory to store data and parameter files
		mkdir -p frealign
		
		iters=-1
		for i in `cd ..; ls frealign_parameters_??; cd - > /dev/null`
		do
			cat ../$i | sed 's|'EXPERIMENT'|'${current_experiment}'|g' > frealign/$i	
			let iters=iters+1
		done
		
		cat ../${parameters_file} | sed 's|'EXPERIMENT'|'${current_experiment}'|g' | sed 's|'ITERATIONS'|'${iters}'|g' > ${parameters_file}
		
		if ! [ -f microscope_parameters ]
		then
			rm -f microscope_parameters
			ln -s ../../microscope_parameters .
		fi
		
		# create links to image stacks
		imagestacks=../${experiment}_noise_`printf %02d $nn`_T??_stack.mrc
		for i in `ls $imagestacks`
		do
			tail=_noise_`printf %02d $nn`
			newtail=_ini_`printf %02d $km`$tail
			current_stack=frealign/`echo $i | sed -e "s|../||" | sed -e "s|$tail|$newtail|"`
			com="rm -f ${current_stack}"
			$com
			com="ln -s `pwd`/$i ${current_stack}"
			$com
		done
		
		# create links to noiseless image stacks
		imagestacks=../${experiment}_noise_00_T??_stack.mrc
		for i in `ls $imagestacks`
		do
			tail=_noise_00
			newtail=_ini_`printf %02d $km`_noise_`printf %02d $nn`
			current_stack=frealign/`echo $i | sed -e "s|../||" | sed -e "s|$tail|$newtail|"`.noiseless
			com="rm -f ${current_stack}"
			$com
			com="ln -s `pwd`/$i ${current_stack}"
			$com
		done

		# create links to .par files
		parfiles=../${experiment}_ini_`printf %02d $km`_T??_01.par
		for i in `ls $parfiles`
		do
			tail=_ini_`printf %02d $km`
			newtail=${tail}_noise_`printf %02d $nn`
			current_stack=frealign/`echo $i | sed -e "s|../||" | sed -e "s|$tail|$newtail|"`
			com="rm -f ${current_stack}"
			$com
			com="ln -s `pwd`/$i ${current_stack}"
			$com
			# com="rm -f ${current_stack}x"
			# $com
			# com="ln -s `pwd`/${current_stack} ${current_stack}x"
			# $com
		done
		
		# generate link to reference file
#		if ! [ -f frealign/${current_experiment}_01.mrc ]
#		then
			rm -f frealign/${current_experiment}_01.mrc
			ln -s ${starting_model} frealign/${current_experiment}_01.mrc
#		fi
		
		jobname=i_`printf %02d $km`_n_`printf %02d $nn`
		
		# run usp & csp
		echo -e "cd ${current_experiment}; ${CSPDIR}/run_usp.sh > ${current_experiment}_USP.log; qsub -N ${jobname} -m ae -v np=48,mydir=`pwd` -l nodes=6:c8 -e `pwd` ${HOME}/scripts/SPA/CSP/run_csp_mpi.sh; cd - > /dev/null" >> ../$uspswarmfile
		echo -e "cd ${current_experiment}; ${CSPDIR}/run_csp.sh init > ${current_experiment}_CSP.log; cd - > /dev/null" >> ../$cspswarmfile
		echo -e "cd ${current_experiment}; qsub -N ${jobname} -m ae -v np=8,mydir=`pwd` -l nodes=1:e2666 -e `pwd` ${HOME}/scripts/SPA/CSP/run_csp_mpi.sh; cd - > /dev/null" >> ../$mpiswarmfile

		echo -e "cd ${current_experiment}; ${CSPDIR}/run_usp.sh > ${current_experiment}_USP.log; ${CSPDIR}/run_csp.sh init > ${current_experiment}_CSP.log; ${CSPDIR}/run_csp.sh iter >> ${current_experiment}_CSP.log; cd - > /dev/null" >> ../$spswarmfile

		cd ..
	done
done

# swarm -f ${uspswarmfile}
# swarm -f ${cspswarmfile}
