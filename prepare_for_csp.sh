#!/bin/bash

if ! [ -n "$1" ] 
then
  echo "Usage: `basename $0` dataset"
  exit 1
#else
#  echo `basename $0` $1; echo
fi

experiment=$1

cd frealign

# concatenate .par from each tilt series into single file

dataset=${experiment}.series
datasetpar=${experiment}_01.parx

rm -f $dataset
rm -f $datasetpar

counter=0

for i in `cat ../micrograph_list | sed -e "s|_aliavg||"`
#for i in `ls *_stack.mrc`
do
        currentTS=${experiment}_${i}
	#currentTS=`echo $i | sed -e 's/_stack.mrc//'`
	# set correct tilt series index
        if [ -f ${currentTS}_stack.mrc ]
        then
            com="$CSPDIR/set_tilt_series_index.sh ${currentTS}_01.par dummy $counter"
            $com
            if [ $counter -eq 0 ]
            then
                    cat dummy >> ${datasetpar}
            else
                    cat dummy | grep -v C >> ${datasetpar}
            fi
            echo -e "$counter \t $currentTS" >> ${dataset}
            let counter=counter+1
        else
            echo ${currentTS}_stack.mrc does not exist.
            rm -f ${currentTS}_01.par
        fi
done

${SPA_DIR}/utils/frealign_fix_indexes.sh ${datasetpar} dummy

mv dummy ${datasetpar}

echo $counter tilt series found.
echo
echo File $datasetpar contains `cat $datasetpar | grep -v C | wc -l` particles.

cd - > /dev/null