#!/bin/bash

n=500

OUTPUT_FILE="many_seeds.txt"

> $OUTPUT_FILE

#for i in $(seq 1 $n); do
#    NUMBER=$RANDOM
#    echo "$RANDOM" >> $OUTPUT_FILE
#done

shuf -i 1-32767 -n 500 > $OUTPUT_FILE

echo "Generated $n random numbers in $OUTPUT_FILE"
