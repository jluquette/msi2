for i in ../msitools/*/*.str_summary.txt; do
    echo $i
    #bsub -q short -W 4:00 "./metrics_msitools.py $i > ${i/.str_summary.txt/}.metrics.txt"
    bsub -q short -W 4:00 "./metrics_msitools.py --min-mapq 36 --min-units 4 --min-supp 15 --max-ref-diff 80 $i > ${i/.str_summary.txt/}.minmapq36_minunits4_minsupp15_maxrefdiff80.metrics.txt"
    bsub -q short -W 4:00 "./candidate_loci.py --min-mapq 36 --min-units 4 --min-supp 15 --max-ref-diff 80 $i > ${i/.str_summary.txt/}.minmapq36_minunits4_minsupp15_maxrefdiff80.candidates.txt"
done
