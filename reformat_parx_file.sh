#!/bin/bash

IFS=';'
img_index=1
exec 0<"$1"
while read -r line
do
	if [ ${line:0:1} == 'C' ]
    then
        echo $line
	else
        echo -n "`printf %7d $img_index`"
        echo -n ${line:7}
        img_index=`expr $img_index + 1`
        echo ${line:0:7}
	fi
done

#./frealign/format_dotpar_line.sh 
