from cosmos.contrib.ezflow.tool import Tool

class align(Tool):
    inputs = [ 'fastq.gz' ]
    outputs = [ 'bam' ]
    mem_req = 4096 + 4*1024
    cpu_req = 2
    time_req = 60
    name = 'BWA MEM paired alignment'

    def cmd(self, i, s, p):
        return """set -o pipefail ;
                  {s[bwa_binary]} mem
                    -M
                    -v 2
                    -R '@RG\tID:{p[readgroup]}\tSM:{p[sample]}\tPL:ILLUMINA'
                    {s[reference_genome]}
                    {i[fastq.gz][0]}
                    {i[fastq.gz][1]}
                    | {s[java_binary]}
                        -Xmx4g
                        -Djava.io.tmpdir={s[tmpdir]}
                        -jar {s[picard_home]}/SortSam.jar
                        SORT_ORDER=coordinate
                        I=/dev/stdin
                        O=$OUT.bam"""


class remdup(Tool):
    inputs = [ 'bam' ]
    outputs = [ 'bam', 'markdup_metrics.txt' ]
    mem_req = 20480
    cpu_req = 1
    time_req = 48*60
    name = 'Remove PCR duplicates'

    def cmd(self, i, s, p):
        inlist = ' '.join('I=' + str(f) for f in i['bam'])
        return """{s[java_binary]}
                    -Xmx20g
                    -Djava.io.tmpdir={s[tmpdir]}
                    -jar {s[picard_home]}/MarkDuplicates.jar
                    VALIDATION_STRINGENCY=LENIENT
                    REMOVE_DUPLICATES=true
                    O=$OUT.bam
                    %s
                    METRICS_FILE=$OUT.markdup_metrics.txt""" % inlist


class index_bam(Tool):
    inputs = [ 'bam' ]
    forward_input = True
    mem_req = 1024
    cpu_req = 1
    time_req = 12*60
    name = "Index BAM"

    def cmd(self, i, s, p):
        return """{s[samtools_binary]} index {i[bam][0]} {i[bam][0]}.bai"""


class sputnik(Tool):
    inputs = [ 'bam' ]
    outputs = [ 'bam' ]
    mem_req = 1024
    cpu_req = 1
    time_req = 12*60
    name = 'sputnik'

    # FIXME: sputnik_wrapper uses the global version of samtools
    def cmd(self, i, s, p):
        return """{s[sputnik_wrapper]}
                    {i[bam][0]}
                    $OUT.bam
                    {p[chrom]}"""


class msitools(Tool):
    inputs = [ 'bam' ]
    outputs = [ 'bam', 'str_summary.txt' ]
    mem_req = 20480
    cpu_req = 1
    time_req = 12*60
    name = 'msitools'

    def cmd(self, i, s, p):
        return """set -o pipefail;
                  {s[samtools_binary]} view -h {i[bam][0]} {p[chrom]}
                  | {s[msitools_script]}
                    --resource_path {s[resource_path]}
                    --summary $OUT.str_summary.txt
                    --flank_bp 10
                    --chr {p[chrom]}
                  | {s[samtools_binary]} view -bSo $OUT.bam -"""


class combine_summaries(Tool):
    """Combine per-chromosome summary files in reference order."""
    inputs = [ 'str_summary.txt' ]
    outputs = [ 'str_summary.txt' ]
    mem_req = 128
    cpu_req = 1
    time_req = 60
    name = 'combine summaries'

    def cmd(self, i, s, p):
        """Assumes s[chr] is in reference order."""
        # Order the files by s['chr'].  To do this, look at the 'chrom' tag
        # of the task that produced each file.
        chr_to_input = {}
        for p in self.parents:
            sums = [ tf for tf in p.output_files if tf.name == 'str_summary.txt' ]
            chr_to_input[p.tags['chrom']] = str(sums[0])

        ilist = ' '.join(chr_to_input[chrom] for chrom in s['chr'])

        # Headers should all match; use the first one.
        return """(head -1 {i[str_summary.txt][0]};
                   tail --quiet -n +2 %s) > $OUT.str_summary.txt""" % ilist


class genotyper(Tool):
    inputs = [ 'str_summary.txt' ]
    outputs = [ 'calls.txt', 'filter_metrics.txt', 'error_profile.txt']
    mem_req = 2048
    cpu_req = 1
    time_req = 60
    name = 'genotyper'

    def cmd(self, i, s, p):
        single_cell = '--single-cell' if p['sample'] in s['single_cell'] else ''
        return """{s[genotyper_script]}
                    --error-distn-file $OUT.error_profile.txt
                    --filter-metrics-file $OUT.filter_metrics.txt
                    --min-mapq 30
                    --min-depth 10
                    %s
                    {i[str_summary.txt][0]} > $OUT.calls.txt""" % single_cell


class join_loci(Tool):
    inputs = [ 'calls.txt' ]
    outputs = [ 'calls.txt' ]
    mem_req = 30720
    cpu_req = 1
    time_req = 12*60
    name = 'join_loci'

    def cmd(self, i, s, p):
        # Input files must be specified as file,sampleID on the command line
        ilist = ''
        for p in self.parents:
            calls = [ tf for tf in p.output_files if tf.name == 'calls.txt' ]
            ilist += '%s,%s' % (str(calls[0]), p.tags['sample'])
        return """Rscript {s[join_script]}
                    $OUT.calls.txt
                    %s""" % ilist
