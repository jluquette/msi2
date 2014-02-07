#!/bin/bash

for i in cortex_bulk neuron_2; do
    bsub -q long -W 72:00 -n 2 -o $i.rmdup.nomapq0.log "gunzip -c ../bams/$i.rmdup.nomapq0.fa.gz | ./sputnik_V6_len5 /dev/stdin | gzip - > $i.rmdup.nomapq0.sputnik.txt.gz"
done
