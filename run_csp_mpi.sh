#!/bin/bash
# This file is multpich3.sh
#
#PBS -k oe
#PBS -N csp
#
# To execute: qsub -v np=160,mydir=`pwd`,cutoff=100 -l nodes=40:dc:centos multipich3.sh
#

cd $mydir
echo "Running on directory: $mydir"

mpdboot -f $PBS_NODEFILE -n `cat $PBS_NODEFILE | wc -l`

# com="/usr/local/mpich/bin/mpirun -machinefile $PBS_NODEFILE -np $np /usr/local/bin/multirun -m $mpirunfile"
# echo $com; $com

# mpirun can only handle a limited number of processes
# if [ $np -lt 161 ]
# then
	cp $PBS_NODEFILE frealign/mpirun.mynodes
	# let npmpirun=np-1
# else
	# cat $PBS_NODEFILE | head -20 > frealign/mpirun.mynodes
	# let npmpirun=160
# fi

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

starting_iters=`grep ^StartingIteration ${parameters_file} | awk '{print $2}'`
if [ -z $starting_iters ]
then
	starting_iters=0
fi

# induce number of inner iterations from configuratuion file
inner_iters=`grep ^BinningForRefinement ${parameters_file} | sed -e "s|BinningForRefinement||" | wc -w | awk '{print $1}'`
if [ -z $inner_iters ]
then
	inner_iters=1
fi
echo Using $inner_iters inner iterations.

if [ ${starting_iters} -eq 0 ]
then
	cp frealign/${experimentName}_01.parx frealign/maps/${experimentName}_CSP_01.parx
fi

# get total number of particle projections for PR threshold calculation
particles=`cat frealign/maps/${experimentName}_CSP_01.parx | grep -v C | wc -l`
prcutoff=`echo "${particles}*${cutoff}/100" | bc`

echo Processing $particles total particles and using only $prcutoff.

# start using all particle projections
UseImagesForRefinementMin=`grep ^UseImagesForRefinementMin ${parameters_file} | awk '{print $2}'`
UseImagesForRefinementMax=`grep ^UseImagesForRefinementMax ${parameters_file} | awk '{print $2}'`
UseImagesForReconstructionMin=`grep ^UseImagesForReconstructionMin ${parameters_file} | awk '{print $2}'`
UseImagesForReconstructionMax=`grep ^UseImagesForReconstructionMax ${parameters_file} | awk '{print $2}'`
# thresh=`cat frealign/${experimentName}_01.parx | grep -v C | awk '{print $12}' | sort -n | head -$prcutoff | tail -1`

echo Using [ $UseImagesForRefinementMin, $UseImagesForRefinementMax ] scan_order range for refinement.
echo Using [ $UseImagesForReconstructionMin, $UseImagesForReconstructionMax ] scan_order range for reconstruction.

# Run the CSP framework using MPI+swarm combination for alignment and reconstruction, respectively.

for (( i = ${starting_iters}; i < $iters; i++ ))
do
	let iteration=i+2
	iformat=`printf %02d $iteration`

	# cleanup
	rm -f csp.done
	rm -f ${experimentName}_CSP_${iformat}.log

	for (( inner = 0; inner < $inner_iters; inner++ ))
	do

	# store inner iteration number
	inneri=`echo "$inner+2" | bc`
	echo $inneri > .csp_current_inner_iteration	

	# get binning for this inner iteration
	binning=`grep ^BinningForRefinement ${parameters_file} | awk '{print $'${inneri}'}'`
	if [ -z $binning ]
	then
		binning=1
	fi
	echo Global iteration $i, Frealign iteration $iteration, Inner iteration $inner, Binning $binning

	# do independent runs for micrographs and particles
	for (( mode = 0; mode < 2; mode++ ))
	do
	#	mode=2;
		
		if (( ( i == 0 ) && ( mode == 0 ) && ( inner == 0 ) ))
		then
			rm -f frealign/scratch/${experimentName}_MPI_state.bin
		fi

		command="qsub -v np=24,mydir=`pwd`,myjob=${experimentName}_CSP_${iformat},iter=$i,mode=$mode -l nodes=6:dc -e `pwd` ${CSPDIR}/mpi_csp.sh"
		command="mpiexec -n $np ${CSPDIR}/bin/csp ${parameters_file} $i $mode ${UseImagesForRefinementMin} ${UseImagesForRefinementMax}"
		echo $command; $command
                
		# Backup experiment logs
		echo Deleting log files ...
		#com="tar cvfz ${experimentName}_CSP_${iformat}_mode_`printf %02d ${mode}`_msearch_n.tgz $listoffiles > /dev/null"
		find frealign/log -maxdepth 1 -type f -name "${experimentName}_CSP_*${iformat}_msearch_n.log" -exec rm '{}' \;

		echo Clearing scratch ...
		find frealign/scratch  -maxdepth 1 -type f -name "${experimentName}_CSP_*.parx" -exec rm '{}' \;

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
		cp ${experimentName}_CSP_${jformat}.parx ${experimentName}_CSP_${jformat}_`printf %02d $inner`.parx

		cd - > /dev/null
		
	done
	done
	
	# retrieve first binning value
	binning=`grep ^BinningForRefinement ${parameters_file} | awk '{print $2}'`
	
	# run FREALIGN (refine shifts and do reconstruction)
	cd frealign;
	let prev=iteration-1
	prevf=`printf %02d $prev`
	
	# # figure out index of last valid reconstruction
	# last_valid_reconstruction=`echo "1+(${i}/${inner_iters})*${inner_iters}" | bc`

	# use proper reference (last valid reconstruction)
	# cp maps/${experimentName}_CSP_`printf %02d ${last_valid_reconstruction}`.mrc maps/${experimentName}_${prevf}.mrc
	cp maps/${experimentName}_CSP_${prevf}.mrc maps/${experimentName}_${prevf}.mrc
	
	# echo FREALIGN iteration $iteration using reference maps/${experimentName}_CSP_`printf %02d ${last_valid_reconstruction}`.mrc

	# override FREALIGN's MASK to only do shift refinement
	mv frealign_parameters_`printf %02d $iteration` frealign_parameters_tmp
	oldmask=`cat frealign_parameters_tmp | grep MASK`
	newmask="MASK                    0,0,0,1,1"
	cat frealign_parameters_tmp | sed -e "s|${oldmask}|${newmask}|" > frealign_parameters_`printf %02d $iteration`

# 	if [ $np -gt 201 ]
# 	then
# 		echo "201 `pwd`/mpirun.mynodes" > mpirun.config
# 	else
		echo "$np `pwd`/mpirun.mynodes" > mpirun.config
# 	fi
	
	# do multiple frealign iterations without changing the orientations
	max_frealign_iterations=2
	for (( k = 0; k < max_frealign_iterations; k++ ))
	do
		echo 1 > .last_iteration_completed
		echo "minscanor $UseImagesForReconstructionMin" > .scanor_limits
		echo "maxscanor $UseImagesForReconstructionMax" >> .scanor_limits
		com="frealign_iter.sh $iteration 0 $cutoff"
		echo "Running FREALIGN reconstruction round $k of $max_frealign_iterations ..."
		echo $com; $com
                rm -f .scanor_limits
		
		# wait until reconstruction is done
		last_iteration_completed=1
		while [ "${last_iteration_completed}" -ne "$iteration" ]
		do
			sleep 3
			last_iteration_completed=`cat .last_iteration_completed`
		done
		
		# replace reference with most up-to-date
		cp maps/${experimentName}_${iformat}.mrc maps/${experimentName}_${prevf}.mrc
		
		# replace parameter files with most up-to-date

		# keep the information in the extended parameter files
		echo Generating extended .parx files ...
		mpirunfile="concatenate_par_files.mpirunfile"
		rm -f $mpirunfile
		counter=1
		for tilt_series in `cat ${experimentName}.series | awk '{print $2}'`
		do
			# generate mpirun script
			if [ $counter -eq 1 ]
			then
				echo \#\!/bin/bash > $mpirunfile
				echo "case \$MP_CHILD in" >> $mpirunfile
			fi
			let countern=counter-1
			echo "${countern})" >> $mpirunfile
			echo "$CSPDIR/concatenate_par_files.sh maps/${tilt_series}_${prevf}.par maps/${tilt_series}_${iformat}.par maps/${tilt_series}_${iformat}_test.par" >> $mpirunfile
			echo "mv maps/${tilt_series}_${iformat}_test.par maps/${tilt_series}_${prevf}.par" >> $mpirunfile
			echo ";;" >> $mpirunfile
			let counter=counter+1
		done
		echo "esac" >> $mpirunfile
		
		com="chmod u+x $mpirunfile"
		echo $com; $com
		let npmpirun=counter-1
		com="/usr/local/mpich2/bin/mpirun -machinefile $PBS_NODEFILE -np $npmpirun /data/${USER}/multirun/multirun -m `pwd`/$mpirunfile > /dev/null"
		echo $com; $com	
		rm -f $mpirunfile
	done
	
	# restore frealign parameter file
	mv frealign_parameters_tmp frealign_parameters_`printf %02d $iteration`
	rm -f mpirun.config	
	
	# if (( ${last_valid_reconstruction} == ${i} - 1 ))
	# then
		echo Saving current reconstruction
		mv maps/${experimentName}_${iformat}.mrc maps/${experimentName}_CSP_${iformat}.mrc
		if [ "$binning" -gt "1" ]
		then
			echo Generate binned reference for alignment maps/${experimentName}_CSP_${iformat}_bin`printf %02d $binning`.mrc
			binvol maps/${experimentName}_CSP_${iformat}.mrc -bin $binning maps/${experimentName}_CSP_${iformat}_bin`printf %02d $binning`.mrc > /dev/null
			# binvol maps/${experimentName}_CSP_01.mrc -bin $binning maps/${experimentName}_CSP_${iformat}_bin`printf %02d $binning`.mrc
		fi

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

	# branch out
	echo Generating extended .parx files ...
	mpirunfile="concatenate_par_files.mpirunfile"
	rm -f $mpirunfile
	counter=1
	for tilt_series in `cat ../${experimentName}.series | awk '{print $2}'`
	do
		# generate mpirun script
		if [ $counter -eq 1 ]
		then
			echo \#\!/bin/bash > $mpirunfile
			echo "case \$MP_CHILD in" >> $mpirunfile
		fi
		let countern=counter-1
		echo "${countern})" >> $mpirunfile
		echo "$CSPDIR/concatenate_par_files.sh ${tilt_series}_`printf %02d $prev`.par ${tilt_series}_${iformat}.par ${tilt_series}_${iformat}_test.par" >> $mpirunfile
		echo "mv ${tilt_series}_${iformat}_test.par ${tilt_series}_${iformat}.par" >> $mpirunfile
		echo ";;" >> $mpirunfile
		let counter=counter+1
	done
	echo "esac" >> $mpirunfile
	
	com="chmod u+x $mpirunfile"
	echo $com; $com 
	let npmpirun=counter-1
	com="/usr/local/mpich2/bin/mpirun -machinefile $PBS_NODEFILE -np $npmpirun /data/${USER}/multirun/multirun -m `pwd`/$mpirunfile > /dev/null"
	echo $com; $com	
	rm -f $mpirunfile

	first=1
	for tilt_series in `cat ../${experimentName}.series | awk '{print $2}'`
	do
		# $CSPDIR/concatenate_par_files.sh ${tilt_series}_`printf %02d $prev`.par ${tilt_series}_${iformat}.par test.par
		# mv test.par ${tilt_series}_${iformat}.par
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

	# compute new tilt_angle_thresh
	# tilt_angle_thresh=`cat ${experimentName}_CSP_${iformat}.parx | grep -v C | awk '{print $12}' | sort -n | head -$prcutoff | tail -1`
	
	mv ${experimentName}_${iformat}_fsc.txt ${experimentName}_CSP_${iformat}_fsc.txt
	mv ${experimentName}_${iformat}_fsc.png ${experimentName}_CSP_${iformat}_fsc.png
	mv ${experimentName}_${iformat}_fsc_all.png ${experimentName}_CSP_${iformat}_fsc_all.png
	mv ${experimentName}_${iformat}_PR.png ${experimentName}_CSP_${iformat}_PR.png

	cd - > /dev/null
	
	# update number of iterations from file
	iters=`grep ^NumberOfIterations ${parameters_file} | awk '{print $2}'`
done

# cleanup
echo Clearing scratch ...
find frealign/scratch -type f -exec rm '{}' \;
	
rm -f frealign/mpirun.mynodes

echo Normal program termination.

mpdallexit

# move log files to current directory
dts=$(date '+%Y%m%d_%H%M%S')
mv ~/${PBS_JOBNAME}.o`echo ${PBS_JOBID} | sed -e "s/.biobos//"` ${experimentName}_CSP_${dts}.log
mv ~/${PBS_JOBNAME}.e`echo ${PBS_JOBID} | sed -e "s/.biobos//"` ${experimentName}_CSP_${dts}.error
