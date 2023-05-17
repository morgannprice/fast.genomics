#!/usr/bin/perl -w
use strict;
use DBI;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use lib "$RealBin/../../PaperBLAST/lib";
use clusterBLASTp;
use pbutils qw{ReadFasta};

my $maxHits1 = 500;
my $maxHits2 = 500;
my $scale = 1;
my $eValue = 1e-3;
my $dataDir = "$RealBin/../data";
my $nCPUs = 12;

my $usage = <<END
clusteredBLASTp.pl -subdb order -in sequence -out hits
  The input file should have just a single query.

Optional arguments:
 -dataDir -- default $dataDir
 -scale (of evalues) -- defaults to $scale
 -evalue (maximum) -- default $eValue
 -maxHits1 -- defaults to $maxHits1
 -maxHits2 -- defaults to $maxHits2
 -nCPUs -- defaults to $nCPUs
 -quiet
END
;

my ($subdb, $inFile, $outFile, $quiet);
die $usage
  unless GetOptions('in=s' => \$inFile,
                    'out=s' => \$outFile,
                    'subdb=s' => \$subdb,
                    'dataDir=s' => \$dataDir,
                    'scale=f' => \$scale,
                    'evalue=f' => \$eValue,
                    'maxHits1=i' => \$maxHits1,
                    'maxHits2=i' => \$maxHits2,
                    'nCPUs=i' => \$nCPUs,
                    'quiet' => \$quiet)
  && @ARGV == 0
  && defined $inFile && defined $outFile && defined $subdb;

my $subdbDir = "$dataDir/$subdb";
die "No such directory: $subdbDir\n" unless -d $subdbDir;
my $subdbFile = "$subdbDir/sub.db";
die "No such file: $subdbFile\n" unless -e $subdbFile;
my $dbh = DBI->connect("dbi:SQLite:dbname=$subdbFile", "", "", { RaiseError => 1 }) || die $DBI::errstr;
my $clusterDb = "$subdbDir/cluster.faa";
die "No such file: $clusterDb.plusdb.pin" unless -e "$clusterDb.plusdb.pin";

my $seqs = ReadFasta($inFile);
die "Must have 1 input sequence\n" unless scalar(keys %$seqs) == 1;
my ($name, $seq) = %$seqs;

my $hits = clusteredBLASTp('query' => $seq, 'clusterDb' => $clusterDb,
                           'dbh' => $dbh,
                           'nCPUs' => $nCPUs,
                           'bin' => $RealBin,
                           'maxHits' => [ $maxHits1, $maxHits2 ],
                           'scale' => $scale,
                           'eValue' => $eValue,
                           'quiet' => defined $quiet);
open(my $fhOut, ">", $outFile) || die "Cannot write to $outFile\n";
foreach my $hit (@$hits) {
  print $fhOut join("\t", $name, $hit->{subject}, $hit->{identity}*100,
                    map $hit->{$_}, qw{alnLength nMismatch nGapOpens qBegin qEnd sBegin sEnd eValue bits})
    . "\n";
}
close($fhOut) || die "Error writing to $outFile\n";
