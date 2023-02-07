#!/usr/bin/perl -w
use strict;
use FindBin qw{$RealBin};

my $usage = "Usage: touchSliced.pl prefix\n"
  . "  where prefix is the stub of the slice names\n";
die $usage unless @ARGV == 1;
my ($pre) = @ARGV;

open (my $fh, "<", "$pre.nSlices") || die "Cannot read $pre.nSlices\n";
my $nSlices = <$fh>;
close($fh) || die "Error reading $pre.nSlices\n";
chomp $nSlices;
die "Invalid #slices in $pre.nSlices\n"
  unless $nSlices >= 1 && $nSlices <= 100;

my $mmseqs = "$RealBin/mmseqs";
system("$mmseqs -h > /dev/null") == 0
  || die "Cannot run $mmseqs\n";

for (my $i = 0; $i < $nSlices; $i++) {
  my $mmdb = $pre.$i;
  my $index = "$mmdb.idx";
  foreach my $file ($mmdb,$index) {
    die "No such file: $pre$i\n"
      unless -e $file;
  }
  my $cmd = "$mmseqs touchdb $mmdb";
  system("$cmd > /dev/null") == 0
    || die "mmseqs failed: $!\n$cmd\n";
}
print STDERR "Touched all $nSlices database slices for $pre\n";
