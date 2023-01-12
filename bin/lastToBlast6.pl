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
  foreach my $field (qw{query subject identity nMismatch nGapOpen qBegin qEnd sBegin sEnd evalue bits}) {
    die "No value for $field in hit" unless defined $hit->{$field};
  }
  # BLAST and usearch round identity and bits to 1 decimal point, so do the same
  print join("\t", $hit->{query}, $hit->{subject},
             sprintf("%.1f",$hit->{identity}), $hit->{nMismatch}, $hit->{nGapOpen},
             $hit->{qBegin}, $hit->{qEnd}, $hit->{sBegin}, $hit->{sEnd},
             $hit->{evalue}, sprintf("%.1f", $hit->{bits}))."\n";
}
