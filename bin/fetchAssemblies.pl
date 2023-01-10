#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib "$RealBin/../../PaperBLAST/lib";
use FetchAssembly qw{CacheAssembly};
use pbutils qw{ReadTable};

my $usage = <<END
fetchAssemblies.pl -table genomes.tsv -dir dir

Fetches the assemblies in the fetch field of the table. This field
should have values like GCF_000092985.1 or GCA_000069185.1. For each
assembly, creates files dir/refseq_assemblyId.{faa,fna,features.tab}.
The assemblyId may not always match the fetch field, so
these are reported in a dir/assemblies.tsv table.
END
;

my ($tableFile, $dir);
die $usage
  unless GetOptions('table=s' => \$tableFile, 'dir=s' => \$dir)
  && @ARGV == 0
  && defined $tableFile && defined $dir;

my @rows = ReadTable($tableFile, ["fetch"]);
my %seen = ();
foreach my $row (@rows) {
  my $fetch = $row->{fetch};
  die "Duplicate fetch for $fetch\n" if exists $seen{$fetch};
  $seen{$fetch} = 1;
  die "Invalid assembly specifier in fetch field: $fetch\n"
    unless $fetch =~ m/^[A-Z]+_\d+[.]\d+$/;
}

my $listFile = "$dir/assemblies.tsv";
open(my $fhList, ">", $listFile) || die "Cannot write to $listFile";
print $fhList "fetch\tgid\n";

FetchAssembly::setFailMode("warning");
my @fetched = ();
foreach my $row (@rows) {
  my $fetch = $row->{fetch};
  my $assembly = CacheAssembly("NCBI", $fetch, $dir);
  if (defined $assembly) {
    push @fetched, $fetch;
    my $gid = $assembly->{gid};
    print $fhList "${fetch}\t${gid}\n";
  }
}
close($fhList) || die "Error writing to $listFile";
print STDERR "Found " . scalar(@fetched) . " of " . scalar(@rows) . " assemblies\n";
