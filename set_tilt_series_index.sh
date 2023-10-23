#!/bin/bash

if ! [ -n "$3" ] 
then
  echo "Usage: `basename $0` input.par output.par number"
  exit 1
#else
#  echo `basename $0` $1 $2 $3; echo
fi  

index=$3

file=tmp_${1}

cat $1 | grep C > $2
cat $1 | grep -v C > $file

awk -v counter=$index '
  /^[ \t]*$/ { print; next }
  /^[ \t]*!/ { print; next }
  {
    printf("%7d%8.2f%8.2f%8.2f%8.2f%8.2f%8.0f%6d%9.1f%9.1f%8.2f%7.2f%8.2f%9d%9.2f%9.2f%9d%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f%9.2f\n",
             $1,$2,$3,$4,$5,$6,$7,counter,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39);
  }
' ${file} >> $2

rm $file
