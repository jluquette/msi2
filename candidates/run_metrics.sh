#!/bin/bash

for i in */*.str_summary.txt; do
    echo $i
    bsub -q short -W 2:00 -o /dev/null "./metrics_msitools.py $i > ${i/.str_summary.txt/.metrics.txt}"
    bsub -q short -W 2:00 -o /dev/null "./metrics_msitools.py --mapq60 $i > ${i/.str_summary.txt/.mapq60_metrics.txt}"
done
