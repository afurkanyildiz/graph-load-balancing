#!/bin/bash

if [ "$1" == "--help" ] || [ "$#" -ne 2 ]
then
  echo -e "parameter #1:"
  echo -e "report name"
  echo -e "parameter #2:"
  echo -e "output file name for new level sizes\n"
fi

start=`grep -n new $1 | head -1 | cut -d':' -f1`
end=`grep -n "rewritten rows:" $1 | head -1 | cut -d':' -f1` ; end=$(( $end -1 ))
sed -ne "${start},${end}p" $1 >> $2
