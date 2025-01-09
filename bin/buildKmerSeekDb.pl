#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use IO::Handle; # for autoflush

my $kmerLen = 15;

my $usage = <<END
Usage: buildKmerSeekDb.pl -in in.fasta -out kmer.db
  Given a fasta file, build a database for searching by first or last
  kmer.

Optional arguments:
-k $kmerLen -- the kmer size
END
;

my $tmpDir = $ENV{TMPDIR} || "/tmp";
my $tmpPre = "$tmpDir/buildKmerSeekDB.$$";

my ($inFile, $dbFile);
die $usage
  unless GetOptions('in=s' => \$inFile,
                    'out=s' => \$dbFile,
                    'k=i' => \$kmerLen)
  && @ARGV ==0
  && defined $inFile
  && defined $dbFile;
die "No such file: $inFile\n" unless -e $inFile;

my $sqlFile = "$RealBin/../lib/kmerSeek.sql";
die "No such file: $sqlFile\n" unless -e $sqlFile;

my $cmd = "$RealBin/fastaToSeeks.pl -in $inFile -out $tmpPre -k $kmerLen";
system($cmd) == 0 || die "Failed command\n$cmd\n$!";

unlink($dbFile);
system("sqlite3 $dbFile < $sqlFile") == 0|| die $!;

open(SQLITE, "|-", "sqlite3", $dbFile) || die "Cannot run sqlite3 on $dbFile";
autoflush SQLITE 1;
print SQLITE ".mode tabs\n";
print STDERR "Loading table KmerFirst\n";
print SQLITE ".import ${tmpPre}_firstk.tsv KmerFirst\n";
print STDERR "Loading table KmerLast\n";
print SQLITE ".import ${tmpPre}_lastk.tsv KmerLast\n";
print SQLITE <<END
  SELECT 'nFirst', COUNT(*) FROM KmerFirst;
  SELECT 'nLast', COUNT(*) FROM KmerLast;
END
;
close(SQLITE) || die "Error running sqlite3 import commands\n";

unlink("${tmpPre}_firstk.tsv");
unlink("${tmpPre}_lastk.tsv");
