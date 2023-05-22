#!/bin/bash

set -eu

if [[ "$1" == "--help" ]] || [[ $# -lt 2 ]]
    then
    echo "Transfer files within a directory"
    echo "To execute:"
    echo "transfer-data.sh <directory-name> <destination-dir>"
    exit 0
fi

dirName=$1
dirDestination=$2
for file in `find ./$dirName/ -type f`; do
    echo transferring $file
    filename=`basename $file`
    pod=nostalgic-nightingale
    namespace=nextflow
    cat $file | kubectl exec -i $pod -n $namespace -- tee $dirDestination/$filename >/dev/null
done
