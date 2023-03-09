#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use neighbor qw{featuresToGenes};
use lib "$RealBin/../../PaperBLAST/lib";
use FetchAssembly qw{CacheAssembly ParseNCBIFeatureFile};
use pbutils qw{ReadTable};

my $usage = <<END
fetchAssemblies.pl -table genomes.tsv -dir dir

Fetches (if needed) the assemblies in the fetch field of the
table. This field should have values like GCF_000092985.1 or
GCA_000069185.1. For each assembly, creates files
dir/refseq_assemblyId.{faa,fna,features.tab}.  The assemblyId may not
always match the fetch field, so these are reported to the fetched
table (by default, genomes.tsv.fetched).  The fetched table also
reports how many protein-coding genes and total genes the assembly
has. Assemblies that cannot be fetched (for instance, that do not have
protein annotations) are omitted from the fetched table.

Optional arguments:
-outTable genomes.tsv.fetched
END
;

my ($tableFile, $outTable, $dir);
die $usage
  unless GetOptions('table=s' => \$tableFile, 'dir=s' => \$dir, '-outTable' => \$outTable)
  && @ARGV == 0
  && defined $tableFile && defined $dir;
die "Not a directory: $dir\n" unless -d $dir;
$outTable = "$tableFile.fetched" if !defined $outTable;

my @rows = ReadTable($tableFile, ["fetch"]);
my %seen = ();
foreach my $row (@rows) {
  my $fetch = $row->{fetch};
  die "Duplicate fetch for $fetch\n" if exists $seen{$fetch};
  $seen{$fetch} = 1;
  die "Invalid assembly specifier in fetch field: $fetch\n"
    unless $fetch =~ m/^[A-Z]+_\d+[.]\d+$/;
}

open(my $fhOut, ">", $outTable) || die "Cannot write to $outTable";
print $fhOut join("\t", qw{fetch gid nGenes nProteinGenes})."\n";

FetchAssembly::setFailMode("warning");
my $nFetched = 0;
foreach my $row (@rows) {
  my $fetch = $row->{fetch};
  my $assembly = CacheAssembly("NCBI", $fetch, $dir);
  if (defined $assembly) {
    my $gid = $assembly->{gid} || die;
    my $featureFile = "$dir/refseq_${gid}.features.tab";
    die "No feature table for $fetch: $featureFile\n" unless -e $featureFile;
    my $features = ParseNCBIFeatureFile($featureFile);
    # Convert features to a list of genes; in particular, combine gene
    # and CDS entries into protein-coding genes
    my ($genes, $warnings) = featuresToGenes($features);
    foreach my $warning (@$warnings) {
      print STDERR join("\t", "Warning", $gid, $warning)."\n";
    }
    my @proteinGenes = grep { $_->{proteinId} ne "" } @$genes;
    print $fhOut join("\t", $fetch, $gid, scalar(@$genes), scalar(@proteinGenes))."\n";
    $nFetched++;
  }
}
close($fhOut) || die "Error writing to $outTable";
print STDERR "Found $nFetched of " . scalar(@rows) . " assemblies\n";
print STDERR "Wrote $outTable\n";
