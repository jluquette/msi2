MSITools is a set of perl scripts developed by Tae-min for mapping Sputnik's
results to a library of known repeat regions in the genome and for summarizing
the mapping.

Changes to Tae-min's code:
    * fix left-flank loop-to-beginning of read bug
    * change the flank_size requirement from 2bp->10bp
    * msitools now creates a file of all reads that were used to support one
      of its calls in a file named 'SUPPORT_reads_<filename>.txt'
    * many, many changes were made to msitools.  it is almost a complete
      rewrite at this point.  the above features are more or less correct,
      but do not reflect all the changes.

Tae-min's database contains two weird STR loci:
    * chr2    242858713   242858723   intergenic      mo
    * chr17   64722990    64722999    intergenic      d
BOTH of these loci are the last loci in their respective chromosomes.  Maybe
there was a corruption or bug when Tae-min was generating the database?  I've
removed both lines from my version of the DB.

NOTE: the analysis was run using a version of msitools that miscomputes
mean searches per read; it divides by the total number of lines processed
instead of 10,000, which is the chunk each mean is computed for.

DEPENDENCY!
http://cpansearch.perl.org/src/BENBOOTH/Set-IntervalTree-0.07/lib/Set/IntervalTree.pm

perl -MCPAN -e shell
install Set::IntervalTree

Must be installed in this directory under 'perllib'.
