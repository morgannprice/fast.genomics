#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib "$RealBin/../../PaperBLAST/lib";
use pbutils qw{ReadFastaEntry};

my $k = 7;

my $usage = <<END
Usage: buildSliced.pl -in in.faa -out db -slices 6

Given a large fasta file of protein sequences, split it into the
requested number of slices and format the mmseqs databases. The
resulting files will be named based on the output argument,
i.e. db0.mmseqs, db1.mmseqs, ..., and associated index files.
These can be queried with searchSliced.pl. The #slices is recorded in
db.nSlices and db will be an empty file.

Optional argument:
  -k $k -- kmer size for each slice (or 0 to let MMseqs2 choose)
END
;

my ($inFaa, $outPre, $nSlices, $debug);
die $usage
  unless GetOptions('in=s' => \$inFaa, 'out=s' => \$outPre, 'slices=i' => \$nSlices,
                    'k=s' => \$k,
                    'debug' => \$debug)
  && @ARGV == 0
  && defined $inFaa && defined $outPre && defined $nSlices;
die "Invalid nSlices -- must be positive\n" unless $nSlices > 0;
die "Invalid nSlices -- limited to 100\n" unless $nSlices <= 100;

my $mmseqs = "$RealBin/mmseqs";
system("$mmseqs -h > /dev/null") == 0
  || die "Cannot run $mmseqs\n";

my $tmp = $ENV{TMPDIR} || "/tmp";
die "Not a directory: $tmp\n" unless -d $tmp;
my $tmpDir = "$tmp/buildSliced.$$";
mkdir($tmpDir) || die "Cannot mkdir $tmpDir";

open(my $fhIn, "<", $inFaa) || die "Cannot read $inFaa\n";

my @sliceFile = ();
my @fhOut = ();
for (my $i = 0; $i < $nSlices; $i++) {
  $sliceFile[$i] = "$tmpDir/faa$i";
  open($fhOut[$i], ">", $sliceFile[$i]) || die "Cannot write to $fhOut[$i]\n";
}

my $state = {};
my $nSeq = 0;
while (my ($header, $seq) = ReadFastaEntry($fhIn, $state)) {
  die "Invalid sequence $seq\n" unless $seq =~ m/^[a-zA-Z*]+$/;
  $seq =~ s/[*]//g;
  $seq = uc($seq);
  my $i = $nSeq++ % $nSlices;
  $fhOut[$i]->print(">", $header, "\n", $seq, "\n");
}
close($fhIn) || die "Error reading $inFaa\n";

for (my $i = 0; $i < $nSlices; $i++) {
  close($fhOut[$i]) || die "Error writing to $sliceFile[$i]";
}

for (my $i = 0; $i < $nSlices; $i++) {
  my $mmDb = $outPre.$i;
  my @cmds = ("$mmseqs createdb $sliceFile[$i] $mmDb",
              "$mmseqs createindex $mmDb $tmpDir -k $k");
  foreach my $cmd (@cmds) {
    print STDERR "Running $cmd\n" if defined $debug;
    system("$cmd > $tmpDir/mmseqs.log") == 0
      || die "mmseqs failed -- $! -- command line\n$cmd\nsee $tmpDir/mmseqs.log\n";
  }
  print STDERR "Created and indexed $mmDb\n" if defined $debug;
}
print STDERR "Created and indexed\n${outPre}0\nthrough\n$outPre".($nSlices-1)."\nwith k-mer size $k\n";

open (my $fhS, ">", "$outPre.nSlices") || die "Cannot write to $outPre.nSlices";
print $fhS $nSlices,"\n";
close($fhS) || die "Error writing to $outPre.nSlices\n";

open (my $fh, ">", $outPre) || die "Cannot write to $outPre";
close($fh) || die "Error writing to $outPre\n";

system("rm -Rf $tmpDir");
