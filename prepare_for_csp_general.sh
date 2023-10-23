#!/bin/bash
# combines all .par files into single extended .parx file
if ! [ -n "$3" ] 
then
  echo "Usage: `basename $0` dataset maps_folder iteration"
  exit 1
#else
#  echo `basename $0` $1; echo
fi

experiment=$1
maps=$2
iter=$3

cd frealign

# concatenate .par from each tilt series into single file
iterformat=`printf %02d $iter`
dataset=${experiment}_${maps}_${iterformat}.series
datasetpar=${maps}/${experiment}_USP_${iterformat}.parx

rm -f $dataset
rm -f $datasetpar

counter=0

for i in `ls *_stack.mrc`
do
	currentTS=`echo $i | sed -e 's/_stack.mrc//'`
	# set correct tilt series index
	com="$CSPDIR/concatenate_par_files.sh ${currentTS}_01.par ${maps}/${currentTS}_${iterformat}.par dummy_1_${maps}_${iterformat}.par"
	$com
	com="$CSPDIR/set_tilt_series_index.sh dummy_1_${maps}_${iterformat}.par dummy_2_${maps}_${iterformat}.par $counter"
	$com
	if [ $counter -eq 0 ]
	then
		cat dummy_2_${maps}_${iterformat}.par | grep C >> ${datasetpar}
		cat ${currentTS}_01.par | grep C >> ${datasetpar}
		cat dummy_2_${maps}_${iterformat}.par | grep -v C >> ${datasetpar}
	else
		cat dummy_2_${maps}_${iterformat}.par | grep -v C >> ${datasetpar}
	fi
	echo -e "$counter \t $currentTS" >> ${dataset}
	let counter=counter+1
done

${SPA_DIR}/utils/frealign_fix_indexes.sh ${datasetpar} dummy_${maps}_${iterformat}

mv dummy_${maps}_${iterformat} ${datasetpar}

rm dummy_1_${maps}_${iterformat}.par dummy_2_${maps}_${iterformat}.par

rm $dataset

echo $counter tilt series found.
echo
echo File $datasetpar contains `cat $datasetpar | grep -v C | wc -l` particles.

cd - > /dev/null