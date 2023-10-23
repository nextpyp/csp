#!/bin/bash
#
# $Id: frealign_iter.sh 142 2009-12-16 03:38:21Z fefo $

# Test whether command-line argument is present (non-empty).
if ! [ -n "$2" ]; then
    echo "Usage: `basename $0` iteration_number cutoff_percentage"
    exit 1
else
	echo "[`date '+%Y%m%d-%H%M%S'`] `basename $0` $@"
fi

if [ $1 -le 1 ]; then
    echo "[ERROR] `basename $0`: Iteration number must be greater than 1."
    exit
fi

# Create directories if needed
mkdir -p scratch
mkdir -p maps
mkdir -p swarm
mkdir -p log

swarmfile=frealign_iter.swarm
dts=`date '+%Y%m%d_%H%M%S'`
logfile=log/${dts}_frealign_iter_`printf %02d $1`.log
echo ${SPA_DIR}/frealign/frealign_iter_implement.sh $1 $2 '>&' ../${logfile} > swarm/${swarmfile}

# Call the job dispatcher. Based in CLUSTER_ID, will call swarm with
# the correct parameters or just run the selected scripts.
cd swarm
${SPA_DIR}/frealign/frealign_jobs_dispatcher.sh -f ${swarmfile} -l nodes=1:o2600:dc
cd - > /dev/null

