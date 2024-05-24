#!/usr/bin/perl -w
use strict;

die "Usage: seekToFasta.pl fastaFile seqPos\n"
  unless @ARGV == 2
  && $ARGV[1] =~ m/^\d+$/;

my ($inFile, $at) = @ARGV;
open(my $fh, "<", $inFile) || die "Cannot read $inFile\n";
seek($fh, $at, 0) # 0 means absolute seek
  || die "Cannot seek to $at\n";
my $header = <$fh> || die "Cannot read header from $inFile at $at\n";
die "Not a header at $at\n" unless $header =~ m/^>/;
print $header;
my $seq = "";
while (my $line = <$fh>) {
  if ($line =~ m/^>/) {
    last;
  } else {
    print $line;
  }
}
close($fh) || die "Error reading $inFile\n";
