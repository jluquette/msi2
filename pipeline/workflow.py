from cosmos.contrib.ezflow.dag import DAG, add_, sequence_, map_, reduce_, split_
from cosmos.contrib.ezflow.tool import INPUT

import tools

def gen_inputs(inputs):
    """Create a list of INPUT()s from the command line arguments."""
    return [ INPUT(path=filename, name='fastq.gz', fmt='gz', tags=tags)
             for filename, tags in inputs.iteritems() ]

def build_dag(input_fastqs, chrs=[]):
    """Build a DAG to handle a FASTQ->polymorphic STR loci."""

    return DAG().sequence_(
        add_(gen_inputs(input_fastqs)),
        reduce_([ 'sample', 'readgroup', 'chunk' ], tools.align),
        reduce_([ 'sample' ], tools.remdup),
        map_(tools.index_bam),
        split_([ ('chrom', chrs) ], tools.sputnik),
        map_(tools.index_bam, stage_name='Index Sputnik BAMs'),
        map_(tools.msitools),
        reduce_([ 'sample' ], tools.combine_summaries),
        #map_(tools.genotyper),
        #reduce_([], tools.join_loci)
    )
