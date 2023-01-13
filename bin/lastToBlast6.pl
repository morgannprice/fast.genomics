#!/usr/bin/perl -w
use strict;
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use neighbor qw{parseLast};

die <<END
Run as a filter on the maf file from lastal.  Writes tab-delimited
"BLAST6" format with field queryId, subjectId, percent identity,
alignment length, number of mismatches, number of gap opens, query
start, query end, subject start, subject end, evalue, bits. All
positions are 1-based.
END
  unless @ARGV == 0;

my $hits = parseLast(\*STDIN);
foreach my $hit (@$hits) {
  # BLAST rounds identity to 2 and bits to 1 decimal point, so do the same
  print join("\t", $hit->{query}, $hit->{subject},
             sprintf("%.2f",$hit->{identity}), $hit->{alnLength},
             $hit->{nMismatch}, $hit->{nGapOpen},
             $hit->{qBegin}, $hit->{qEnd}, $hit->{sBegin}, $hit->{sEnd},
             $hit->{evalue}, sprintf("%.1f", $hit->{bits}))."\n";
}
