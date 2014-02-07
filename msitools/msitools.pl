#!/usr/bin/env perl

# The previous version of this program split its original input into one
# file per chromosome and worked on a single chromosome at a time.  It also
# required the input to be sorted for efficient lookup of repeat records.
# This version is almost completely rewritten by Joe to neither split the
# input nor to require the input to be sorted.

# Using the ~8M repeat locus database produced by Tae-min, this program
# requires ~11GB of RAM to store the full database and human reference
# sequences in memory.  Typical RAM usage for analyzing the data in this
# experiment is ~50G.  One sample required as much as 52G.  It is usually
# more practical to split your analysis by chromosome, using the --chr
# option.

# TODO: TEST! compare against previous output.  will not be identical
# due to fixing the negative flanking index bug as well as removing all
# lines for repeats that have no supporting reads.

# NOTE: This program standardizes on chromosome strings withOUT "chr".

use strict;
use warnings;
use Getopt::Long;
use Cwd qw(getcwd abs_path);
use List::Util qw(max);

my $script_dir;
BEGIN {
    $script_dir = abs_path($0);
    $script_dir =~ s/msitools.pl//g;
}
use lib $script_dir . 'perllib';
use Set::IntervalTree;

my $cmdline = "$0 " . join(" ", @ARGV);

my $flank_bp = 2;
my $short_test = 0;
my $debug = 0;
my $resource_path = getcwd;
my @chrarray = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y);
my $inputfile = '/dev/stdin';
undef my $summary;
GetOptions("flank_bp=i" => \$flank_bp,
           "input=s" => \$inputfile,
           "summary=s" => \$summary,
           "resource_path=s" => \$resource_path,
           "chr=s" => sub { @chrarray = map { s/chr//; $_; } split(",", $_[1]) },
           "debug" => \$debug)  # EXTREMELY verbose.  Do not use on large input
    or die("error parsing arguments");

die "--input is required" if not defined $inputfile;
die "--summary is required" if not defined $summary;
die "--resource_path is required" if not defined $resource_path;

print STDERR "reading SAM records from $inputfile\n";
print STDERR "flank_bp=$flank_bp\n";
print STDERR "summary file=$summary\n";

# Build a hash of accepted chromosomes
print STDERR "chroms=" . join(", ", @chrarray) . "\n";
my %chrhash;
foreach my $chr (@chrarray) {
    $chrhash{$chr} = 1;
}
print STDERR "Using " . scalar(keys %chrhash) . " chromosomes...\n";



# Reads the full set of human reference sequences into memory
print STDERR "Reading chromosome sequences into RAM..\n";
my %chrSeq;
foreach my $chr (@chrarray) {
    print STDERR "reading chromosome $chr.. ";
    open (F, "<", "$resource_path/chr$chr.fa") or die;
    my $null = <F>;  # strip off the first >chrXXX line
    my $nbytes = 0;
    while (<F>) {
        my $line = $_;
        chomp $line;
        $chrSeq{$chr} .= $line;
        $nbytes += length($line);
    }
    print STDERR "$nbytes bytes.\n";
    close F; 
}
print STDERR "done.\n";



# Build a record of all known repeats on this chromosome
# This used to be done one chrom at a time.  Changed by Joe.
print STDERR "Building repeat index..\n";
my %repeatdb;
my %repeatdb_stats;
open (F, "<", "$resource_path/WGRef_7892585MS_withGENEcategory_FINAL_withMSTYPE_hg19.txt") or die;
while (<F>) {
    chomp;
    my @element = split("\t");
    my $chr = $element[0];
    $chr =~ s/chr//g;
    # Only load repeats in the specified chromosome list
    if (exists $chrhash{$chr}) {
        $element[5] =~ s/\s//g;
        # There was a significant performance improvement for using an array
        # instead of a hash for these records, albeit at a cost to code
        # understandability.
        my $rec = [
            $chr,         # chr
            $element[1],  # start
            $element[2],  # end
            $element[3],  # region (exonic, intronic, intergenic, 3utr, 5utr)
            $element[4],  # seq
            $element[5],  # unit (mono, di, tri, tetra)
            uc(substr($chrSeq{$chr}, $element[1] - $flank_bp - 1, $flank_bp)),
            uc(substr($chrSeq{$chr}, $element[2] , $flank_bp)),  # flank2
            [],           # passing lens
            [],           # passing strands
            [],           # passing mapqs
            [],           # failing lens
            [],           # failing strands
            [],           # failing mapqs
            []            # failing reasons
        ];

        # The library doesn't make these values easily queryable.
        # track the number of records and the max record high value.
        $repeatdb_stats{$chr} = [ 0, 0, 0 ] if not defined $repeatdb_stats{$chr};
        ++$repeatdb_stats{$chr}[0];
        $repeatdb_stats{$chr}[1] = max($repeatdb_stats{$chr}[1], $element[2]+1);

        # If either flank sequence is a subsequence of the repeat, discard the
        # locus.  We have to do extra work to correctly assign supporting reads
        # to these loci.  For now it's easier to just ignore them.
        if ($rec->[4] =~ $rec->[6] or $rec->[4] =~ $rec->[7]) {
            if ($debug and 0) {
                print STDERR "flank matched repeat sequence, discarding locus:\n";
                print STDERR "chr=$rec->[0], start=$rec->[1], end=$rec->[2]\n";
                print STDERR "flank1=$rec->[6]\n";
                print STDERR "flank2=$rec->[7]\n";
                print STDERR "seq=$rec->[4]\n";
                print STDERR "marked=$rec->[6]|$rec->[4]|$rec->[7]\n";
                print STDERR "complete seq=";
                print STDERR uc(substr($chrSeq{$chr},
                                       $element[1] - $flank_bp - 1,
                                       $element[2] - $element[1] + 1 + 2*$flank_bp));
                print STDERR "\n";
            }
            ++$repeatdb_stats{$chr}[2];
        } else {
            # Intervals are half open: [a, b).
            $repeatdb{$chr} = Set::IntervalTree->new if not defined $repeatdb{$chr};
            $repeatdb{$chr}->insert($rec, $element[1], $element[2]+1);
        }

    }
}
close F;

while (my ($k, $v) = each(%repeatdb)) {
    print STDERR "chr$k: $repeatdb_stats{$k}[0] repeat records in [0, $repeatdb_stats{$k}[1]), $repeatdb_stats{$k}[2] loci dropped since they contained flanking sequence\n";
}
print STDERR "done.\n";



my $sputnik_in_header = 0;
my $inheader = 1;
my $nread = 0;
my $nsupp = 0;
my $num_searches = 0;
my $num_notindb = 0;
print STDERR "Reading input data..\n";
open (F, "<", $inputfile) or die;
while (<F>) {
    # Always print out SAM header lines
    if (/^@/) {
        print $_;
        if (/PN:sputnik/) {
            $sputnik_in_header = 1;
        }
        next;
    } elsif ($inheader) {
        if (!$sputnik_in_header) {
            print STDERR "ERROR: This SAM does not have a Sputnik \@PG tag ";
            print STDERR "in its header.  Aborting.\n";
            exit 1;
        }
        # The first non-header line.  Before possibly printing anything else,
        # add a @PG tag for msitools.
        $inheader = 0;
        my @md5sum = split " ", `md5sum $0`;
        print "\@PG\tID:msitools\tPN:msitools\tCL:$cmdline\tVN:$md5sum[0]\n";
        print "\@CO\tTG:ls\tTY:i\tDS:Reference repeat locus start coordinate\n";
        print "\@CO\tTG:le\tTY:i\tDS:Reference repeat locus end coordinate\n";
        print "\@CO\tTG:lr\tTY:Z\tDS:Reference repeat locus region: one of intergenic, intronic, 3utr, 5utr, or exonic\n";
        print "\@CO\tTG:f1\tTY:Z\tDS:Left flanking sequence, length depends on --flank_bp parameter to msitools\n";
        print "\@CO\tTG:f2\tTY:Z\tDS:Right flanking sequence, length depends on --flank_bp parameter to msitools\n";
        print "\@CO\tTG:ln\tTY:Z\tDS:Reference repeat locus nucleotide sequence\n";
        print "\@CO\tTG:fr\tTY:Z\tDS:Reason for filtering the read.  One of: pass, cigar, unit, flank1, flank2\n";
    }

    ++$nread;
    if ($nread % 10000 == 0) {
        print STDERR "$nread lines processed, $nsupp supporting reads, ";
        print STDERR $num_searches / 10000 . " mean searches per read.\n";
        $num_searches = 0;
    }

    chomp;
    my @element = split("\t"); 
    my $chrom = $element[2];
    $chrom =~ s/chr//g;
    next if not exists $chrhash{$chrom};  # Skip if not in chr list

    # Parse the optional SAM tags (which start at element 11) for the
    # Sputnik repeat tags.  All 3 tags must be present.
    undef my $repunit;
    undef my $repstart;
    undef my $repend;
    for (@element[11 .. $#element]) {
        if (/^ru:Z:(.+)/) {
            $repunit = $1;
        } elsif (/^rs:i:(.+)/) {
            $repstart = $1;
        } elsif (/^re:i:(.+)/) {
            $repend = $1;
        }
    }
    die("Required tag 'ru' not found in SAM record") unless defined $repunit;
    die("Required tag 'rs' not found in SAM record") unless defined $repstart;
    die("Required tag 're' not found in SAM record") unless defined $repend;

    my $replen = scalar($repend - $repstart + 1);
    my $read_start = $element[3];
    my $startRepeat = $read_start + $repstart - 1;
    my $endRepeat = $read_start + $repend - 1;
    my $strand = (int $element[1] & 0x10) == 0 ? '+' : '-';
    my $mapq = int $element[4];
    my $cigar = $element[5];
    my $this_flank1 = substr($element[9], $repstart-$flank_bp-1, $flank_bp);
    my $this_flank2 = substr($element[9], $repend, $flank_bp);
    
    my $overlapping_repeats = $repeatdb{$chrom}->fetch($startRepeat, $endRepeat);
    if (scalar(@$overlapping_repeats) == 0) {
        if ($debug) {
            print STDERR "DEBUG: not in db: $chrom:[$startRepeat, $endRepeat)\n";
        }
        ++$num_notindb;
    }
    for my $record (@$overlapping_repeats) {
        my ($chr, $start, $end, $region, $seq, $unit, $flank1, $flank2, $passing_lens, $passing_strands, $passing_mapqs, $failing_lens, $failing_strands, $failing_mapqs, $failing_reasons) = @$record;

        my $reflen = $end - $start + 1;
        my $refdiff = $replen - $reflen;

        ++$num_searches;
        # Assert that an insertion or deletion is reported in the CIGAR string.
        # E.g., in a stretch like TTTATTTTTTTTTT; an A->T polymorphism or reads
        # with a common Illumina homopolymer sequencing error at the A would
        # support a +(T)6 insertion.  However, the local alignment against ref
        # should reveal that this is a 1bp mismatch, not a 6bp insertion.
        # CIGAR format: (number)(letter)[(number)(letter)...]
        my $cigar_support = 0;
        my $cigar_ref = "", my $cigar_ref_at_repeat = "";
        my $cigar_read = "", my $cigar_read_at_repeat = "";
        if (not $cigar eq "*") {  # Don't try to parse unaligned reads
            my @cigar_elements = split(qr/([A-Z])/, $cigar);
            for (my $i = 0; $i < scalar(@cigar_elements); $i = $i + 2) {
                my $bps  = $cigar_elements[$i];
                my $mode = $cigar_elements[$i+1];
                if ($mode eq "M" or $mode eq "S" or $mode eq "H") {
                    # M=Match OR mismatch, S=softclip, H=hardclip
                    $cigar_ref = $cigar_ref . $mode x $bps;
                    $cigar_read = $cigar_read .+ $mode x $bps;
                } elsif ($mode eq "D") {
                    # Deletion relative to ref
                    $cigar_ref = $cigar_ref . $mode x $bps;
                    $cigar_read = $cigar_read . "-" x $bps;
                } elsif ($mode eq "I") {
                    # Insertion relative to ref
                    $cigar_ref = $cigar_ref . "-" x $bps;
                    $cigar_read = $cigar_read . $mode x $bps;
                } else {
                    die("unsupported CIGAR mode $mode in CIGAR=$cigar");
                }
            }
    
            # XXX: also require 10M for each flanking region?
            my $complen = $reflen > $replen ? $reflen : $replen;  # max
            $cigar_ref_at_repeat = substr($cigar_ref, $start - $read_start, $complen);
            $cigar_read_at_repeat = substr($cigar_read, $repstart - 1, $complen);
            # Ms should match.  The reason they might not be equal is the substr()
            # calls.  We could get a wacky replen (e.g., the A->T above) which may
            # cause complen to exceed the length of either CIGAR expansion.
            # sum(-1 for D, +1 for I) should match refdiff
            my $ref_ms = $cigar_ref_at_repeat =~ tr/M//;
            my $read_ms = $cigar_read_at_repeat =~ tr/M//;
            my $ds = $cigar_ref_at_repeat =~ tr/D//;
            my $is = $cigar_read_at_repeat =~ tr/I//;
            if ($ref_ms == $read_ms and $is - $ds == $refdiff) {
                $cigar_support = 1;
            }
        }

        my $filter_reason = "";
        if (not $repunit eq $unit) {
            $filter_reason = "unit";
        # joe: fix the case where $element[9]-$flank_bp-1 is negative,
        # which causes substr to select the end of the read
        } elsif (not $flank1 eq $this_flank1 or $repstart - $flank_bp - 1 < 0) {
            $filter_reason = "flank1";
        } elsif (not $flank2 eq $this_flank2) {
            $filter_reason = "flank2";
        # I want flank failures before cigar, because insufficient flank seq
        # will also cause a cigar error when the substr is taken.
        } elsif (not $cigar_support) {
            $filter_reason = "cigar";
        } else {
            $filter_reason = "pass";
        }

        if ($filter_reason eq "pass") {
            push @$passing_lens, $replen;
            push @$passing_strands, $strand;
            push @$passing_mapqs, $mapq;
            ++$nsupp;
        } else {
            push @$failing_lens, $replen;
            push @$failing_strands, $strand;
            push @$failing_mapqs, $mapq;
            push @$failing_reasons, $filter_reason;
        }

        # Print out all BAM records incl. fails with the failing filter reason
        # Print out a valid SAM record with 7 new tags describing the
        # repeat locus matching this read.
        print "$_\tls:i:$start\tle:i:$end\tlr:Z:$region\t";
        print "f1:Z:$flank1\tf2:Z:$flank2\tln:Z:$seq\tfr:Z:$filter_reason\n";
        if ($debug) {
            my $seq_context = uc(substr($chrSeq{$chrom},
                                        $start - $flank_bp - 1,
                                        $end - $start + 2*$flank_bp + 1));
            print STDERR "read: $chrom:$element[3], STR region: $startRepeat-$endRepeat\n";
            print STDERR "record: $start\t$end\t$replen\t@element\n";
            print STDERR "$flank1\t$flank2\n";
            print STDERR "$this_flank1\t$this_flank2\n";
            print STDERR "repeatSeq=" . " " x $flank_bp . "$seq\n";
            print STDERR "chrSeq=   $seq_context\n";
            print STDERR "CIGAR=$cigar\n";
            print STDERR "cigar_ref=  $cigar_ref\n";
            print STDERR "cigar_read= $cigar_read\n";
            print STDERR "cigar_ref_at_repeat=  $cigar_ref_at_repeat\n";
            print STDERR "cigar_read_at_repeat= $cigar_read_at_repeat\n";
            print STDERR "PASS lens here: @$passing_lens\n";
            print STDERR "PASS strands here: @$passing_strands\n";
            print STDERR "PASS mapqs here: @$passing_mapqs\n";
            print STDERR "FAIL lens here: @$failing_lens\n";
            print STDERR "FAIL strands here: @$failing_strands\n";
            print STDERR "FAIL mapqs here: @$failing_mapqs\n";
            print STDERR "FAIL reasons here: @$failing_reasons\n";
            print STDERR "FILTER REASON: $filter_reason\n";
            print STDERR "-" x 80 . "\n";
        }
        # Don't allow one repeat stretch in a read to match multiple repeat
        # loci; one READ can still support multiple loci, however, if it has
        # multiple repeat stretches.
        last;
    }
}
close F;
print STDERR "done.\n";
print STDERR "$nread total reads, $nsupp placed, $num_notindb not found in repeat db.\n";



print STDERR "Writing summary to $summary.. ";
open (OUTPUT, ">$summary");
print OUTPUT "chr\tstart\tend\tunit\tregion\tflank1\tflank2\tsequence\tpassRepArray\tpassStrandArray\tpassMapQArray\tfailRepArray\tfailStrandArray\tfailMapQArray\tfailReasonArray\n";
# Don't loop over keys %repeatdb, this way preserves chrom order
for my $chr (@chrarray) {  
    # Hack.  The IntervalTree library does not export a method for iterating
    # over its nodes, so we do a search that returns all possible nodes.
    # Recall that repeatdb_stats[1] is the max upper bound.
    my $all_recs = $repeatdb{$chr}->fetch(0, $repeatdb_stats{$chr}[1]);

    # Sort output records by start and end position (all recs have the same chr)
    for my $record (sort { $a->[1] <=> $b->[1] or $a->[2] <=> $b->[2] } @$all_recs) {
        my ($chr, $start, $end, $region, $seq, $unit, $flank1, $flank2, $passing_lens, $passing_strands, $passing_mapqs, $failing_lens, $failing_strands, $failing_mapqs, $failing_reasons) = @$record;
        print OUTPUT "$chr\t$start\t$end\t$unit\t$region\t$flank1\t$flank2\t$seq";
        print OUTPUT "\t" . join(",", @$passing_lens);
        print OUTPUT "\t" . join(",", @$passing_strands);
        print OUTPUT "\t" . join(",", @$passing_mapqs);
        print OUTPUT "\t" . join(",", @$failing_lens);
        print OUTPUT "\t" . join(",", @$failing_strands);
        print OUTPUT "\t" . join(",", @$failing_mapqs);
        print OUTPUT "\t" . join(",", @$failing_reasons) . "\n";
    }
}
close OUTPUT;
print STDERR "done.\n";
