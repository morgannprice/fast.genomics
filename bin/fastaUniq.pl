#!/usr/bin/perl -w
use strict;
use Digest::MD5 qw{md5_hex};
use lib "/usr2/people/mprice/src/PaperBLAST/lib";
use pbutils qw{ReadFastaEntry};

die "Run as a filter -- silently skips any repeated sequences\n"
  if @ARGV > 0;

my $state ={};
my %seen = ();
my %seen2 = ();
my $nSkip = 0;
my $nWritten = 0;
while (my ($header,$sequence) = ReadFastaEntry(\*STDIN, $state)) {
  my $md5 = md5_hex(uc($sequence));
  my $key2 = md5_hex(uc($sequence) . length($sequence));
  if (exists $seen{$md5} && exists $seen2{$key2}) {
    $nSkip++;
    next;
  }
  $seen{$md5} = 1;
  $seen2{$key2} = 1;
  print ">$header\n$sequence\n";
  $nWritten++;
}
print STDERR "Wrote $nWritten sequences and skipped $nSkip\n";
