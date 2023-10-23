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

swarmfile=${experiment}.swarm
rm -f ${swarmfile}

# for expk in $experiments
# do
#for k in 0 1 3 5
#for km in 0 1 3 5
for km in 0
do
#	for kp in 220 200 20 002
	for kp in 002
	do
		# create perturbed .parx file
		perturbed_file=${experiment}_ini_M_`printf %02d $km`_P_`printf %03d $kp`.parx
		#perturbed_file=${experiment}_ini_`printf %02d $km`.parx
		if ! [ -f ${perturbed_file} ]
		then
			echo Generating file ${perturbed_file} 
			${CSPDIR}/perturb_groundtruth_par_file.sh -i ${experiment}_01.parx -m $km $km $kp $kp $kp -o ${perturbed_file}
		else
			echo File ${perturbed_file} exists.
		fi

#		for nn in 0 15 30 45 60
		for nn in 0
		do
			# create folder for current experiment
			current_experiment=${experiment}_ini_M_`printf %02d $km`_P_`printf %03d $kp`_noise_`printf %02d $nn`
			#current_experiment=${experiment}_ini_`printf %02d $km`_noise_`printf %02d $nn`
			mkdir -p ${current_experiment}

			# cd ${current_experiment}/frealign/maps
			# for iter in 2 3 4 5 6 7 8 9 10
			# do
				# ${CSPDIR}/concatenate_par_files.sh ${current_experiment}_CSP_${iter}.parx ${current_experiment}_CSP_${iter}.parxt test.par
				# mv test.par ${current_experiment}_CSP_${iter}.parx
				# rm ${current_experiment}_CSP_${iter}.parxt
			# done
			# cd - > /dev/null
			# continue		

			cd ${current_experiment}

			# create frealign directory to store data and parameter files
			mkdir -p frealign
			
			iters=0
			for i in `cd ..; ls frealign_parameters_??; cd - > /dev/null`
			do
				cat ../$i | sed 's|'EXPERIMENT'|'${current_experiment}'|g' > frealign/$i	
				let iters=iters+1
			done
			
			cat ../${parameters_file} | sed 's|'EXPERIMENT'|'${current_experiment}'|g' | sed 's|'ITERATIONS'|'${iters}'|g' > ${parameters_file}
			
			if ! [ -f microscope_parameters ]
			then
				ln -s `pwd`/../microscope_parameters .
			fi
			
			# create link to image stack
			if ! [ -f frealign/${current_experiment}_stack.mrc ]
			then
				com="ln -s `pwd`/../${experiment}_noise_`printf %02d $nn`_stack.mrc frealign/${current_experiment}_stack.mrc"
				echo $com; $com
			else
				echo File frealign/${current_experiment}_stack.mrc already exists.
			fi
			
			# create link to perturbed .parx file
			rm -f frealign/${current_experiment}_01.par*
			cp `pwd`/../${perturbed_file} frealign/${current_experiment}_01.parx
			ln -s `pwd`/frealign/${current_experiment}_01.parx frealign/${current_experiment}_01.par
			# else
				# echo File frealign/${current_experiment}_01.parx already exisits.
			# fi
			# # copy .parx file from original run
			# cp ../../GroEL_03/${current_experiment}/frealign/${current_experiment}_01.parx frealign/${current_experiment}_01.parx
			# rm -f frealign/${current_experiment}_01.par
			# ln -s `pwd`/frealign/${current_experiment}_01.parx frealign/${current_experiment}_01.par
			
			# # generate link to reference file
			# if ! [ -f frealign/${current_experiment}_01.mrc ]
			# then
				rm -f frealign/${current_experiment}_01.mrc
				ln -s ${starting_model} frealign/${current_experiment}_01.mrc
			# fi
			
			# run usp
			# echo "cd frealign; frealign_dispatcher.sh $iters 100 > ../${current_experiment}_USP.log; cd - > /dev/null" > ${current_experiment}.swarm
			echo -e "cd ${current_experiment}; ${CSPDIR}/run_csp_experiment.sh ${experiment}_01.parx > ${current_experiment}_CSP.log; cd - > /dev/null" >> ../$swarmfile
			# swarm -f ${current_experiment}.swarm

			cd ..
			# run csp
			# run_csp_experiment.sh ${log_file}
			

			# expk=parameters.config.GroEL_02_noise_${nn}_ini_${k}
			# basename=`echo $expk | sed -e "s|${parameters_file}.||"`
			# echo "[`date +%Y-%m-%d@%H:%M:%S`] Running experiment: ${basename} ($expk)"
			# cp $expk ${parameters_file}
			# log_file=${basename}_output.log
			# ./run_csp_experiment.sh $log_file
			# echo ""
			# echo ""
		done
	done
done

# swarm -f ${swarmfile}
