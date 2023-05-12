#!/bin/bash

args=("$@")
VAR1=$(ls -d "${args[0]}"/*f*q.gz | head -1)
SEQLEN=$(gzip -cd ${VAR1} | head -2 | tail -1 | wc -c)

echo "$(($SEQLEN - 2))"
