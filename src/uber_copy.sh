#!/bin/bash

#input="/scratch/c.mcbsd1/nextflow/testing_uber/trace.txt"
input=$1
#workDir="/scratch/c.mcbsd1/nextflow/testing_uber/work/"
workDir=$2
#workingDir="/scratch/c.mcbsd1/nextflow/testing_uber/output"
workingDir=$3
linkDir=$4
keyword=$5

rest1="*/*"
while IFS= read -r line
do
  #var1=`awk '$4 == $keyword { print $2 }'`
  var1=`grep -n "$keyword" | grep "COMPLETED" | cut -f2`
done < "$input"

array=($var1)

for i in "${array[@]}"
do
  cmd1=`realpath $workDir/$i$rest1`
  cmd2=`ln -fs $cmd1 $linkDir`
  if ! [[ -L "$file" ]]
  then
    echo "$cmd2"
  fi
  #cmd1=`ln -s $workDir$i $linkDir`
done
