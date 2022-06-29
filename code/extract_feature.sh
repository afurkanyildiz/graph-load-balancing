#!/bin/bash

if [ "$1" == "--help" ] || [ "$#" -ne 1 ]
then
  echo "\$1: pattern to look for in all the reports"
fi

while read matrix
do
  name=`echo $matrix | cut -d',' -f2`
  echo "matrix: $name"
  grep "$1" report."$name"
  echo  ""
#done < matrix.ID.name
done < torso

