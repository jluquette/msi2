for i in ../sputnik/*.sputnik.txt.gz; do
    x=$(basename $i .rmdup.nomapq0.sputnik.txt.gz)
    bsub -q short -W 12:00 -o $x/$x.log -R 'rusage[mem=52000]' "gunzip -c $i | ./msitools.pl --input /dev/stdin --outprefix $x/$x.flank10 --flank_bp 10"
done
