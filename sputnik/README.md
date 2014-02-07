Results of Sputnik on the FASTAs in ../bams.  Running Sputnik (and gzipping
the input and output) takes a while (~3 days).

When copying data from /files, two of the Sputnik results files reported
md5sum: Input/output error.

cortex_bulk.rmdup.nomapq0.sputnik.txt.gz
neuron_2.rmdup.nomapq0.sputnik.txt.gz

rerun.sh was used to rerun these samples, and is also the exact way the
uncorrupt samples were processed.


The version of Sputnik used is: sputnik_V6_len5.  This corresponds to a
specific setting of the length parameter in the Sputnik code.  I don't know
if the version of the code in this repository has that variable set to 5.
