#!/usr/bin/env python

"""
Build a JSON file containing pipeline parameters and an input file manifest.
The manifest file should work as an input to start the pipeline in Cosmos.

The input FASTQ file names must specify:
    1. a sample identifier
    2. a *single* read group (i.e., FASTQ files must already be broken into
       read groups)
    3. a mate identifier: it doesn't matter how these are specified as long
       as it is consistent within a single sample.
    4. an optional chunk identifier: FASTQs should be broken into reasonable
       chunks (5-10M reads each) so that alignment can be parallelized.
"""

import os
import argparse
import sys
import json
import os.path
from collections import defaultdict

def get_path_dependencies(config):
    """Return a dictionary of files and paths that depend on the user
    specified base, resource, annovar and result paths.  If any of these
    paths or files do not exist, throw an error."""

    base_path = os.path.abspath(config['base_path'])
    resource_path = os.path.abspath(config['resource_path'])
    result_path = os.path.abspath(config['result_path'])
    dep_path = os.path.join(config['dep_path'])

    deps = {
        # Some base paths
        'base_path': base_path,
        'resource_path': resource_path,
        'result_path': result_path,
        'dep_path': dep_path,

        # Files that must be in the resource path
        'reference_genome': os.path.join(resource_path, 'hg19.fa'),
        'reference_dict': os.path.join(resource_path, 'hg19.dict'),

        # Scripts and compiled programs used by various tools
        'java_binary': '/usr/bin/java',
        'bwa_binary': os.path.join(dep_path, 'bwa-0.7.5a', 'bwa'),
        'samtools_binary': os.path.join(dep_path, 'samtools-0.1.19', 'samtools'),
        'sputnik_wrapper': os.path.join(base_path, 'sputnik', 'wrap_sputnik.sh'),
        'msitools_script': os.path.join(base_path, 'msitools', 'msitools.pl'),
        'genotyper_script': os.path.join(base_path, 'genotyper', 'genotyper.py'),
        'join_script': os.path.join(base_path, 'candidates', 'join_loci.R'),

        # Miscellany
        'picard_home': os.path.join(dep_path, 'picard-tools-1.95'),
        'tmpdir': os.path.join(result_path, 'tmp'),
    }

    missing = [ (k, v) for k, v in deps.items() if not os.path.exists(v) ]
    if missing:
        raise RuntimeError('missing required files/paths: ' + str(missing))

    return deps


# Files were split into 4M read chunks and renamed by a custom python script.
# The file names are in the following format:
#    MDA-2_ACAGTG_L003_R1_C007.fastq.gz
#    <sampleID>_<barcode>_<lane>_<mate>_<chunk>.fastq.gz
# The various IDs may contain "-" but may not contain underscores.
# All sample IDs that are MDA amplified must begin with "MDA-".
# All sample IDs not derived from single cells must contain "bulk"
# FASTQ files must be contained in a directory named by the unique machine
# run or flowcell ID.
# min_file_size in bytes
def handle_input(fastq_manifest, min_file_size=100000):
    with open(fastq_manifest, 'r') as f:
        fastq_files = [ x.strip() for x in f.readlines() ]

    result = {}
    # get rid of "small" files (param: min_file_size), but keep
    # in mind that if one FASTQ file is removed due to size, its mated file
    # must also be removed
    too_small = defaultdict(bool)  # default value of bool() is False
    for f in fastq_files:
        f = os.path.abspath(f)
        if not f.lower().endswith('.fastq.gz'):
            raise RuntimeError('expected files to end in ".fastq.gz", got ' + f)
        if not os.path.exists(f):
            raise RuntimeError('FASTQ file "%s" does not exist' % f)

        run_id = os.path.split(os.path.split(f)[0])[1]
        fields = os.path.basename(f).rstrip('.gz').rstrip('.fastq').split('_')
        try:
            #sample, readgroup, mate, chunk = [ x.strip() for x in fields ]

            sample, barcode, lane, mate, chunk = [ x.strip() for x in fields ]
            # This is a crappy way to do things, but easy
            readgroup = "_".join([run_id, lane, barcode])
            sc = False if 'bulk' in sample else True
            amp = "MDA" if sample.startswith("MDA-") else "PCR"
            small = os.path.getsize(f) < min_file_size
            result[f] = { 'sample': sample,
                          'readgroup': readgroup,
                          'run_id': run_id,
                          'barcode': barcode,
                          'lane': lane,
                          'mate': mate,
                          'chunk': chunk,
                          'single_cell': sc,
                          'amplification': amp,
                          'too_small': small,
                          'size': os.path.getsize(f) }
            # OR the values: because one mate may be too small and one may
            # be fine.  if assign is used, the larger mate may overwrite the
            # previous smaller mate's decision.
            too_small[sample + readgroup + chunk] |= small
        except KeyError:
            sys.stderr.out('Error parsing FASTQ file name:\n')
            sys.stderr.out('fields=' + str(fields))
            exit(1)

    '''
    for k,v in too_small.items():
        if v:
            print(k + "\t" + str(v))
    '''

    # Second pass to determine total sample size, in bytes.  Drop the FASTQs
    # failing the file size cutoff
    sample_size = defaultdict(int)
    for k, v in result.items():
        if too_small[v['sample'] + v['readgroup'] + v['chunk']]:
            del(result[k])
        else:
            sample_size[v['sample']] += v['size']

    # Now all FASTQs failing the cutoff should be removed.
    for k, v in result.items():
        v['sample_size'] = sample_size[v['sample']]

    return result


def is_int(x):
    try:
        int(x)
        return True
    except ValueError:
        return False


def split_region_by(reg_len, unit):
    """Split a region of length `reg_len` into roughly equal segments of
    size `unit`.  The unit is used to determine the number of desired regions
    and is then modified so that each region is as equally sized as possible."""
    num_regions = int(round(float(reg_len) / unit))
    if num_regions == 0:
        return [(1, reg_len)]
    equi_unit, rem = divmod(reg_len, num_regions)
    return [ (1 + k*equi_unit + min(k, rem), (k+1)*equi_unit + min(k+1, rem))
             for k in range(num_regions) ]
 

def get_chrs(config, split_chr, chr22):
    """Sets config['chr'] and config['split_chr']."""

    # Get sequence names and lengths from the reference dictionary.
    with open(config['reference_dict'], 'r') as refdict:
        seqs = []
        for line in refdict:
            if line.strip().startswith('@SQ'):
                fields = dict([ x.strip().split(':')[0:2]
                                for x in line.split('\t') if ':' in x ])
                chrom = fields['SN'].replace('chr', '').upper()
                if is_int(chrom) or chrom in [ 'X', 'Y' ]:
                    if not chr22 or chrom == '22':
                        seqs += [ (fields['SN'], int(fields['LN'])) ]

    config['chr'] = [ x[0] for x in seqs ]

    # Use full chromosomes as the minimum unit for parallelization of genotyping
    if split_chr.lower().startswith('chr'):
        config['split_chr'] = \
            dict(zip(config['chr'], [ [k] for k in config['chr'] ]))
    else:
        regs = defaultdict(list)
        for contig, length in seqs:
            this_regions = split_region_by(length, int(split_chr))
            for start, stop in this_regions:
                regs[contig] += [ '%s:%d-%d' % (contig, start, stop) ]
        config['split_chr'] = regs
 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description= \
        "Prepare a JSON configuration file describing an analysis.  This " \
        "config file can be supplied as the sole argument to Cosmos to run " \
        "an analysis.")
    #parser.add_argument('FASTQ_files', metavar='file', nargs='+', type=str,
        #help='Gzipped FASTQ input files')
    parser.add_argument('fastq_manifest', metavar='file',
        help='File containing FASTQ files to analyze, one per line.');
    parser.add_argument('--name', metavar='string', type=str, required=True,
        help='A name for this analysis, unique amongst analyses in the ' \
             'Cosmos database.')
    parser.add_argument('--base-path', metavar='path', type=str,
        default=os.getcwd(),
        help='Path containing the pipeline code, resources, etc.')
    parser.add_argument('--resource-path', metavar='path', type=str,
        default=os.path.join(os.getcwd(), 'resources'),
        help='Path to a resource directory, which must contain the human '\
             'reference FASTAs, indexes for BWA, and the STR locus database.')
    parser.add_argument('--dep-path', metavar='path', type=str,
        default=os.path.join(os.getcwd(), 'dependencies'),
        help='Path containing dependencies such as BWA, Picard, etc.')
    parser.add_argument('--result-path', metavar='path', type=str,
        default=os.path.join(os.getcwd(), 'results'),
        help='Directory to store pipeline results.')
    parser.add_argument('--split-chr', type=str, metavar='INT', default='chr',
        help='Potentially for parallelizing at a sub-chrom level.  ' \
             'Set to "chr" to disable.  Currently ignored.')
    parser.add_argument('--chr22', action='store_true', dest='chr22',
        default=False, help='Use only chromosome 22 for a quick test run.')
    parser.add_argument('--save-config', metavar='path', type=str,
        default='/dev/stdout', help='Write the config file to this path.')
    args = parser.parse_args()


    # Get a list of all files/paths that must exist in specified paths
    config = get_path_dependencies(vars(args))

    config['name'] = args.name

    # We only have an STR database for hg19
    config['genome_version'] = 'hg19'

    # Parse input files and determine additional characteristics
    config['input'] = handle_input(args.fastq_manifest)
    config['samples'] = list(set(v['sample'] for v in config['input'].values()))
    config['single_cell'] = list(set(v['sample'] for v in config['input'].values() if v['single_cell']))

    # Get a list of chromosomes for parallelization.  Also a list of subchrom
    # segments if any part of the analysis ever wants a subchrom split.
    get_chrs(config, args.split_chr, args.chr22)

    # Write config file to output if specified.  Do not overwite a previous
    # config file if it exists--require the user to delete it.
    output = args.save_config
    if os.path.exists(output) and output != '/dev/stdout':
        print("ERROR: Configuration file already exists: " + output)
        exit(1)

    with open(output, 'w') as fout:
        json.dump(config, fout, sort_keys=True, ensure_ascii=True,
                  indent=4, encoding='ascii')
