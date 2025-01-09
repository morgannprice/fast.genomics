#!/usr/bin/perl -w

use strict;
die "Run as a filter\n" unless @ARGV == 0;

my $name = undef;
my $seq = undef;

while(<STDIN>) {
  chomp;
  if (m/^>(.*)$/) {
    print "$name\t$seq\n" if defined $name;
    $name = $1;
    $seq = ""
  } else {
    $seq .= $_;
  }
}
print "$name\t$seq\n" if defined $name;
