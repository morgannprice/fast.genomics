#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib "$RealBin/../../PaperBLAST/lib";
use pbutils qw{ReadFasta reverseComplement};

my $usage = <<END
get16S.pl [ -prefix GCF_00351045.1 ] -in ind/refseq_GCF_003351045.1 ~/ind/refseq_GCA_021372535.1

  Extracts the nucleotide sequences for annotated 16S sequences, using
  in.fna and in.features.tab

  Optionally adds prefix: to the front of each extracted sequence's name

  Accepts multiple inputs
END
  ;

my @in;
my $prefix;
die $usage
  unless GetOptions('prefix=s' => \$prefix, 'in=s{1,}' => \@in)
  && @ARGV == 0
  && @in > 0;
foreach my $in (@in) {
  foreach my $suffix ("fna", "features.tab") {
    die "No such file: $in.$suffix\n" unless -e "$in.$suffix";
  }
}
foreach my $in (@in) {
  my $ntseq = ReadFasta("$in.fna");
  open(my $fh, "<", "$in.features.tab")
    || die "Cannot read $in.features.tab";
  while (my $line = <$fh>) {
    next unless $line =~ m/^rRNA\t/ && $line =~ m/16S|SSU|small/;
    chomp $line;
    my @F = split /\t/, $line;
    my $scaffold = $F[6];
    my $begin = $F[7];
    my $end = $F[8];
    my $strand = $F[9];
    my $locusTag = $F[16];
    unless ($begin =~ m/^\d+/ && $end =~ m/^\d+$/) {
      print STDERR "Skipping $F[16] with non-standard coordinates at $begin $end\n";
      next;
    }
    die "$locusTag on nonexisting scaffold $scaffold\n" unless exists $ntseq->{$scaffold};
    my $seq = $ntseq->{$scaffold};
    die "Invalid begin end for $locusTag\n"
      unless $begin >= 1 && $begin < $end && $end <= length($seq);
    my $subseq = substr($seq, $begin-1, $end-$begin+1);
    $subseq = reverseComplement($subseq) if $strand eq "-";
    my $header = $locusTag;
    $header = $prefix . ":" . $locusTag if defined $prefix;
    print ">" . $header . "\n" . $subseq . "\n";
  }
}


