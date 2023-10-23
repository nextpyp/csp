#!/bin/bash
#
# Constrained Single Particle refinement
# Requirements:
#	- parameters.config
#	- dataset_01.mrc
#	- dataset_*{_stack.mrc,.par}

if ! [ -n "$1" ] 
then
  echo "Usage: `basename $0` execmode (init,iter)"
  exit 1
fi

execmode=$1

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
binning=`grep ^BinningForRefinement ${parameters_file} | awk '{print $2}'`
if [ -z $binning ]
then
	binning=1
fi

if [ "$execmode" == "init" ]
then
	# build FREALIGN directory structure
	mkdir -p frealign
	mkdir -p frealign/log
	mkdir -p frealign/scratch
	mkdir -p frealign/maps
	mkdir -p frealign/swarm

	# Data preparation (produces single .parx file for dataset and collects list of tilt series)
	${CSPDIR}/prepare_for_csp.sh ${experimentName}

	# Build the new stacks, one per micrograph and one per particle.
	${CSPDIR}/create_new_stacks_parallel.sh

	# Initial model
	starting_model=frealign/${experimentName}_01.mrc
	if ! [ -f ${starting_model} ]
	then
		echo ERROR - File ${starting_model} does not exist.
		exit 1
	else
		# cp ${starting_model} frealign/maps/${experimentName}_CSP_01.mrc
		#if ! [ -f frealign/maps/${experimentName}_CSP_01.mrc ]
		#then
		#	ln -s `pwd`/${starting_model} frealign/maps/${experimentName}_CSP_01.mrc
		#fi
		rm -f frealign/maps/${experimentName}_CSP_01.mrc
		ln -s `pwd`/frealign/${experimentName}_01.mrc frealign/maps/${experimentName}_CSP_01.mrc
		if [ "$binning" -gt "1" ]
		then
		    binvol frealign/maps/${experimentName}_CSP_01.mrc -bin $binning frealign/maps/${experimentName}_CSP_01_bin`printf %02d $binning`.mrc
		fi
	fi

	# Initial alignments
	starting_alignments=frealign/${experimentName}_01.parx
	if ! [ -f ${starting_alignments} ]
	then
		echo ERROR - File ${starting_alignments} does not exist.
		exit 1
	else
		cp ${starting_alignments} frealign/maps/${experimentName}_CSP_01.parx
		cp ${starting_alignments} frealign/maps/${experimentName}_USP_01.par
		cp ${starting_alignments} frealign/${experimentName}_01.par
	fi
	
else

	# Run the CSP framework using MPI+swarm combination for alignment and reconstruction, respectively.

	for (( i = 0; i < $iters; i++ ))
	do
		let iteration=i+2
		iformat=`printf %02d $iteration`
	
		# cleanup
		rm -f csp.done
		rm -f ${experimentName}_CSP_${iformat}.log

		# do independent runs for micrographs and particles
		for (( mode = 0; mode < 2; mode++ ))
		do
#			mode=2;
			
			# figure out the number of processes needed
			totalNumberOfMicrographs=0
			totalNumberOfParticles=0
			for ts in `ls frealign/*_stack.mrc`
			do
				currentTS=`echo $ts | sed -e 's/_stack.mrc//'`
				numberOfMicrographs=`cat ${currentTS}_01.par | grep -v C | awk '{print $17}' | sort -n | tail -1`
				numberOfParticles=`cat ${currentTS}_01.par | grep -v C | awk '{print $14}' | sort -n | tail -1`
				totalNumberOfMicrographs=`echo "$totalNumberOfMicrographs+$numberOfMicrographs+1" | bc -l`
				totalNumberOfParticles=`echo "$totalNumberOfParticles+$numberOfParticles+1" | bc -l`
			done
			if [ "$mode" == "0" ]
			then
				totalNumberOfJobs=80
			else
				totalNumberOfJobs=$totalNumberOfParticles
			fi
			
#			totalNumberOfJobs=$totalNumberOfMicrographs

			# determine largest number of avaialble nodes
			c8_nodes_available=`freen | head -15 | grep e2666 | awk '{split($6,a,"/"); print a[1];}'`
			c4_nodes_available=`freen | head -15 | grep :dc | awk '{split($6,a,"/"); print a[1];}' | sort -n | tail -1`
			c2_nodes_available=`freen | head -15 | grep "o2800 " | awk '{split($6,a,"/"); print a[1];}' | sort -n | tail -1`
			c2_m4_nodes_available=`freen | head -15 | grep "o2800 " | awk '{split($5,a,"/"); print a[1];}' | sort -n | tail -1`

			let c8_process_available=c8_nodes_available*8
			let c4_process_available=c4_nodes_available*4
			let c2_process_available=c2_nodes_available*2

			# decide which nodes to use based on availability
			if (( ( $c8_process_available > 0 ) && ( $c8_process_available >= $c4_process_available ) && ( $c8_process_available >= $c2_process_available ) ))
			then
				total_procs_available=$c8_process_available
				procs_per_node=8
				pool=c8
			else
				if (( ( $c4_process_available > 0 ) && ( $c4_process_available >= $c2_process_available ) ))
				then
					total_procs_available=$c4_process_available
					procs_per_node=4
					pool=dc
				else
					total_procs_available=$c2_process_available
					procs_per_node=8
					if [ "${c2_m4_nodes_available}" -eq "${c2_nodes_available}" ]
					then
						pool=c8
					else
						pool=c8
					fi
				fi
			fi

			echo Nodes available [ c8 = $c8_nodes_available, dc = $c4_nodes_available, o2800 = $c2_nodes_available ], using $pool
			
			maxNumberOfJobs=`batchlim | grep "norm " | awk '{print $2}'`
			userRemainingJobs=`jobload $USER | grep "Core Total" | awk '{print $3}'`
			
			# let maxNumberOfJobs=maxNumberOfJobs-userRemainingJobs
			maxNumberOfJobs=`echo "scale=0; $maxNumberOfJobs/8" | bc`
			echo $maxNumberOfJobs maximum number of jobs allowed

			if (( ( $maxNumberOfJobs > $total_procs_available ) && ( $total_procs_available > 0 ) ))
			then
				let maxNumberOfJobs=total_procs_available
			fi
			
			if [ "$totalNumberOfJobs" -lt "$procs_per_node" ]
			then
				nodes=1
			else
				if [ "$totalNumberOfJobs" -lt "$maxNumberOfJobs" ]
				then
					nodes=`echo "($totalNumberOfJobs / $procs_per_node + .5 )" | bc -l | awk '{split($1,a,"."); print a[1]}'`
				else
					nodes=`echo "($maxNumberOfJobs / $procs_per_node + .5 )" | bc -l | awk '{split($1,a,"."); print a[1]}'`
				fi
			fi
			
			let np=nodes*procs_per_node

			if (( ( i == 0 ) && ( mode == 0 ) ))
			then
				rm -f frealign/scratch/${experimentName}_MPI_state.bin
			fi
			
			echo $totalNumberOfJobs jobs to run, processes per node $procs_per_node, $pool nodes requested = $nodes

			# if $nodes is empty, assign a fix value
			if [ "$nodes" -eq "" ]
			then
				nodes=10
			fi
			
			command="qsub -v np=24,mydir=`pwd`,myjob=${experimentName}_CSP_${iformat},iter=$i,mode=$mode -l nodes=6:dc -e `pwd` ${CSPDIR}/mpi_csp.sh"
			echo $command

			# if (( $i > 0 ))
			# then
				jobnumber=`$command | sed -e 's/.biobos//'`
				echo Iteration $iformat, mode $mode submitted ${jobnumber}.biobos
			# fi
	
			# exit 1
			
			echo Wait until job completed ...
			while ! [ -f csp.done ]
			do
				sleep 3
			done
			rm -f csp.done

			# Backup experiment logs
			echo Deleting log files ...
			#com="tar cvfz ${experimentName}_CSP_${iformat}_mode_`printf %02d ${mode}`_msearch_n.tgz $listoffiles > /dev/null"
			find frealign/log -type f -name "${experimentName}_CSP_*${iformat}_msearch_n.log" -exec rm '{}' \;

			echo Clearing scratch ...
			find frealign/scratch -type f -name "${experimentName}_CSP_*.parx" -exec rm '{}' \;
		
		done
	
		# submit refinement and get job number
			
		# run FREALIGN (refine shifts and do reconstruction)
		cd frealign;
		let prev=iteration-1
		prevf=`printf %02d $prev`
		cp maps/${experimentName}_CSP_${prevf}.mrc maps/${experimentName}_${prevf}.mrc
		echo 1 > .last_iteration_completed

		# override FREALIGN's MASK to only do shift refinement
		mv frealign_parameters_`printf %02d $iteration` frealign_parameters_tmp
		oldmask=`cat frealign_parameters_tmp | grep MASK`
		newmask="MASK                    0,0,0,1,1"
		cat frealign_parameters_tmp | sed -e "s|${oldmask}|${newmask}|" > frealign_parameters_`printf %02d $iteration`
		
		com="frealign_iter.sh $iteration 100"
		echo Running frealign reconstruction ...
		echo $com; $com
		
		# wait until reconstruction is done
		last_iteration_completed=1
		while [ "${last_iteration_completed}" -ne "$iteration" ]
		do
			sleep 3
			last_iteration_completed=`cat .last_iteration_completed`
		done
		
		# restore frealign parameter file
		mv frealign_parameters_tmp frealign_parameters_`printf %02d $iteration`
		
		echo Saving reconstruction results ...
		mv maps/${experimentName}_${iformat}.mrc maps/${experimentName}_CSP_${iformat}.mrc
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
		
		# cd maps
		# # rename output files
		# iformat=`printf %02d $i`
		# mv ${experimentName}_${iformat}.par ${experimentName}_CSP_${iformat}.par
		# mv ${experimentName}_${iformat}_fsc.txt ${experimentName}_CSP_${iformat}_fsc.txt
		# mv ${experimentName}_${iformat}_fsc.png ${experimentName}_CSP_${iformat}_fsc.png
		# mv ${experimentName}_${iformat}_fsc_all.png ${experimentName}_CSP_${iformat}_fsc_all.png
		# mv ${experimentName}_${iformat}_PR.png ${experimentName}_CSP_${iformat}_PR.png
		# mv ${experimentName}_${iformat}.mrc ${experimentName}_CSP_${iformat}.mrc
		# cd .. # maps
		# cd .. # frealign
	done

	# cleanup
	# rm -f frealign/scratch/*
	find frealign/scratch -type f -exec rm '{}' \;
	
fi

echo Normal program termination.