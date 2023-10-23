#!/bin/bash
# This file is multpich3.sh
#
#PBS -k oe
#PBS -N csp
#
# To execute: qsub -v np=160,mydir=`pwd` -l nodes=40:dc:centos multipich3.sh
#

cd $mydir
echo "Running on directory: $mydir"

mpdboot -f $PBS_NODEFILE -n `cat $PBS_NODEFILE | wc -l`

# com="/usr/local/mpich/bin/mpirun -machinefile $PBS_NODEFILE -np $np /usr/local/bin/multirun -m $mpirunfile"
# echo $com; $com

cp $PBS_NODEFILE frealign/mpirun.mynodes

#!/bin/bash
#
# Constrained Single Particle refinement
# Requirements:
#	- parameters.config
#	- dataset_01.mrc
#	- dataset_*{_stack.mrc,.par}

# The configuration file is always the same
parameters_file=parameters.config

if ! [ -f ${parameters_file} ]
then
	echo ERROR - File ${parameters_file} not found.
	exit 1
fi

# Grep the basename from the <ParFile> in the parameters file. It is
# the second argument in the line starting with 'ParFile' without '.
experimentName=`grep ^ParFile ${parameters_file} | awk '{print $2}'`
iters=`grep ^NumberOfIterations ${parameters_file} | awk '{print $2}'`

inner_iters=3

cp frealign/${experimentName}_01.parx frealign/maps/${experimentName}_CSP_01.parx

# Run the CSP framework using MPI+swarm combination for alignment and reconstruction, respectively.

for (( i = 14; i < 19; i++ ))
do
	let iteration=i+2
	iformat=`printf %02d $iteration`

	# cleanup
	rm -f csp.done
	rm -f ${experimentName}_CSP_${iformat}.log

	for (( inner = 0; inner < 3; inner++ ))
	do
	
	# do independent runs for micrographs and particles
	for (( mode = 0; mode < 2; mode++ ))
	do
	#	mode=2;
		
		if (( ( i == 0 ) && ( mode == 0 ) && ( inner == 0 ) ))
		then
			rm -f frealign/scratch/${experimentName}_MPI_state.bin
		fi
		
		command="qsub -v np=24,mydir=`pwd`,myjob=${experimentName}_CSP_${iformat},iter=$i,mode=$mode -l nodes=6:dc -e `pwd` ${CSPDIR}/mpi_csp.sh"
		command="mpiexec -n $np ${CSPDIR}/bin/csp ${parameters_file} $i $mode"
		echo $command; $command

		# Backup experiment logs
		echo Deleting log files ...
		#com="tar cvfz ${experimentName}_CSP_${iformat}_mode_`printf %02d ${mode}`_msearch_n.tgz $listoffiles > /dev/null"
		find frealign/log -type f -name "${experimentName}_CSP_*${iformat}_msearch_n.log" -exec rm '{}' \;

		echo Clearing scratch ...
		find frealign/scratch -type f -name "${experimentName}_CSP_*.parx" -exec rm '{}' \;

		# concatenate current .par files into .parx file
		cd frealign/maps
		let j=i+1
		jformat=`printf %02d $j`
		first=1
		for tilt_series in `cat ../${experimentName}.series | awk '{print $2}'`
		do
			if [ $first -eq 1 ]
			then
				cat ${tilt_series}_${jformat}.par > ${experimentName}_CSP_${jformat}.parx
				first=0
			else
				cat ${tilt_series}_${jformat}.par | grep -v C >> ${experimentName}_CSP_${jformat}.parx
			fi
		done
		${SPA_DIR}/utils/frealign_fix_indexes.sh ${experimentName}_CSP_${jformat}.parx dummy
		mv dummy ${experimentName}_CSP_${jformat}.parx
		cd - > /dev/null
		
	done
	done
	
	# run FREALIGN (refine shifts and do reconstruction)
	cd frealign;
	let prev=iteration-1
	prevf=`printf %02d $prev`
	
	# # figure out index of last valid reconstruction
	# last_valid_reconstruction=`echo "1+(${i}/${inner_iters})*${inner_iters}" | bc`

	# use proper reference (last valid reconstruction)
	# cp maps/${experimentName}_CSP_`printf %02d ${last_valid_reconstruction}`.mrc maps/${experimentName}_${prevf}.mrc
	cp maps/${experimentName}_CSP_${prevf}.mrc maps/${experimentName}_${prevf}.mrc
	echo 1 > .last_iteration_completed
	
	# echo FREALIGN iteration $iteration using reference maps/${experimentName}_CSP_`printf %02d ${last_valid_reconstruction}`.mrc

	# override FREALIGN's MASK to only do shift refinement
	mv frealign_parameters_`printf %02d $iteration` frealign_parameters_tmp
	oldmask=`cat frealign_parameters_tmp | grep MASK`
	newmask="MASK                    0,0,0,1,1"
	cat frealign_parameters_tmp | sed -e "s|${oldmask}|${newmask}|" > frealign_parameters_`printf %02d $iteration`

	echo "$np `pwd`/mpirun.mynodes" > mpirun.config
	
	com="frealign_iter.sh $iteration 0 100"
	echo Running frealign reconstruction ...
	echo $com; $com
	
	# rm -f mpirun.config
	
	# wait until reconstruction is done
	last_iteration_completed=1
	while [ "${last_iteration_completed}" -ne "$iteration" ]
	do
		sleep 3
		last_iteration_completed=`cat .last_iteration_completed`
	done
	
	# restore frealign parameter file
	mv frealign_parameters_tmp frealign_parameters_`printf %02d $iteration`
	
	# if (( ${last_valid_reconstruction} == ${i} - 1 ))
	# then
		echo Saving current reconstruction
		mv maps/${experimentName}_${iformat}.mrc maps/${experimentName}_CSP_${iformat}.mrc
	# else
		# echo Ignoring this reconstruction and using last valid one: ${experimentName}_CSP_`printf %02d ${last_valid_reconstruction}`.mrc
		# rm -f maps/${experimentName}_${iformat}.mrc
		# rm -f maps/${experimentName}_CSP_${iformat}.mrc
		# ln -s `pwd`/maps/${experimentName}_CSP_`printf %02d ${last_valid_reconstruction}`.mrc `pwd`/maps/${experimentName}_CSP_${iformat}.mrc
	# fi
	rm -f maps/${experimentName}_${prevf}.mrc
		
	cd - > /dev/null	
	
	# substitute newly computed shifts into extended .parx files
	let prev=iteration-1
	cd frealign/maps
	first=1
	for tilt_series in `cat ../${experimentName}.series | awk '{print $2}'`
	do
		$CSPDIR/concatenate_par_files.sh ${tilt_series}_`printf %02d $prev`.par ${tilt_series}_${iformat}.par test.par
		mv test.par ${tilt_series}_${iformat}.par
		if [ $first -eq 1 ]
		then
			cat ${tilt_series}_${iformat}.par > ${experimentName}_CSP_${iformat}.parx
			first=0
		else
			cat ${tilt_series}_${iformat}.par | grep -v C >> ${experimentName}_CSP_${iformat}.parx
		fi
	done
	${SPA_DIR}/utils/frealign_fix_indexes.sh ${experimentName}_CSP_${iformat}.parx dummy
	mv dummy ${experimentName}_CSP_${iformat}.parx
	
	mv ${experimentName}_${iformat}_fsc.txt ${experimentName}_CSP_${iformat}_fsc.txt
	mv ${experimentName}_${iformat}_fsc.png ${experimentName}_CSP_${iformat}_fsc.png
	mv ${experimentName}_${iformat}_fsc_all.png ${experimentName}_CSP_${iformat}_fsc_all.png
	mv ${experimentName}_${iformat}_PR.png ${experimentName}_CSP_${iformat}_PR.png

	cd - > /dev/null
	
done

# cleanup
find frealign/scratch -type f -exec rm '{}' \;
	
rm -f frealign/mpirun.mynodes

echo Normal program termination.

mpdallexit

# move log files to current directory
mv ~/csp.o`echo ${PBS_JOBID} | sed -e "s/.biobos//"` ${experimentName}_CSP.log
mv ~/csp.e`echo ${PBS_JOBID} | sed -e "s/.biobos//"` ${experimentName}_CSP.error
