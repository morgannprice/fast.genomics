#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $kmerLen = 15;
my $usage = <<END
Usage: fastaToSeeks.pl -in in.fasta -out prefix [ -k $kmerLen ]

Writes tab-delimited tables prefix_firstk.tsv and prefix_lastk.tsv
which indexes the 1st or last kmer of each sequence to its seek
position in the input fasta file. Sequences of kmer length or shorter
are included as a single kmer in the first file. Leading Ms are ignored
in the input sequences.
END
;

my ($inFile, $outPre);
die $usage
  unless GetOptions('in=s' => \$inFile,
                    'out=s' => \$outPre,
                    'k=i' => \$kmerLen)
  && @ARGV == 0
  && defined $inFile && defined $outPre;
die "Invalid kmer length\n" unless $kmerLen > 1 && $kmerLen <= 50;

open(my $fhIn, "<", $inFile) || die "Cannot read $inFile\n";
my $out1 = $outPre."_firstk.tsv";
my $out2 = $outPre."_lastk.tsv";
open (my $fhOut1, ">", $out1) || die "Cannot write to $out1\n";
open (my $fhOut2, ">", $out2) || die "Cannot write to $out2\n";

my $at = 0;
my $header = undef;
my $seq = "";

for(;;) {
  my $lineAt = tell($fhIn);
  my $line = <$fhIn>;
  if (!defined $line || $line =~ m/^>/) {
    if (defined $header) {
      die "No sequence for $header\n" if $seq eq "";
      $seq =~ s/^[mM]+//; # ignore leading methionines
      my $kmer1 = substr($seq, 0, $kmerLen);
      print $fhOut1 join("\t", $kmer1, $at)."\n";
      if (length($seq) > $kmerLen) {
        my $kmer2 = substr($seq, length($seq) - $kmerLen);
        print $fhOut2 join("\t", $kmer2, $at)."\n";
      }
    }
  }
  if (!defined $line) {
    last;
  } else {
    chomp $line;
    if ($line =~ m/^>/) {
      $header = $line;
      $header =~ s/^>//;
      die "Empty header" if $header eq "";
      $seq = "";
      $at = $lineAt; # seek for beginning of this sequence
    } else {
      die "No header yet\n" unless defined $header;
      die "Invalid sequence line $line\n" unless $line =~ m/^[a-zA-Z*]*$/;
      $line =~ s/[*]//g; # ignore * characters (usually from translation problems)
      $line = uc($line);
      $seq .= $line;
    }
  }
} # end loop over lines
close($fhIn) || die "Error reading $inFile\n";
close($fhOut1) || die "Error writing to $out1\n";
close($fhOut2) || die "Error writing to $out2\n";
