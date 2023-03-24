#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use lib "$RealBin/../lib";
use neighbor qw{parseTaxString featuresToGenes orderToSubName};
use lib "$RealBin/../../PaperBLAST/lib";
use FetchAssembly qw{ParseNCBIFeatureFile};
use pbutils qw{ReadFastaEntry ReadTable};
use clusterBLASTp;
my $nCPUs = $ENV{MC_CORES} || (`egrep -c '^processor' /proc/cpuinfo`);

my $minCoverage = 0.9;
my $minIdentity = 0.7;
my $nSlices = 4;

my $usage = <<END
./buildSubDbs.pl -genomes genomes.tsv [ -fetched genomes.tsv.fetched ] -in indir -out outdir

Creates the sub databases for each order, in subdirectories of the
output directory. Also writes to outdir/SubDbs.tsv

The genomes table must include the fields fetch, (gtdb) accession,
ncbi_assembly_name, gtdb_taxonomy, ncbi_taxonomy, and
ncbi_strain_identifiers. The fetched table must include the fields
fetch and gid. The input directory should contain all the fetched
files.

More arguments (optional):
-nCPUs $nCPUs -- number of CPUs to use for cd-hit
-minCoverage $minCoverage -- minimum alignment coverage (both ways)
  for the cluster seed and the sequence
-minIdentity $minIdentity -- minimum fractional identity to cluster
-nSlices -- number of slices to use for mmseqs databases
END
;

my ($genomeFile, $fetchedFile, $inDir, $outDir);
die $usage
  unless GetOptions('genomes=s' => \$genomeFile,
                    'fetched=s' => \$fetchedFile,
                    'in=s' => \$inDir,
                    'out=s' => \$outDir,
                    'nCPUs=i' => \$nCPUs,
                    'minCoverage=f' => \$minCoverage,
                    'minIdentity=f' => \$minIdentity,
                    'slices=i' => \$nSlices)
  && defined $genomeFile
  && defined $inDir
  && defined $outDir
  && @ARGV == 0;
$fetchedFile = "$genomeFile.fetched" if !defined $fetchedFile;

die "minCoverage is out of range\n" unless $minCoverage <= 1;
die "minIdentity is out of range\n" unless $minIdentity >= 0.5 && $minIdentity <= 1;
die "slices is out of range\n" unless $nSlices >= 1 && $nSlices <= 100;

my $cdhit = "$RealBin/cd-hit";
die "No such executable: $cdhit\n" unless -x $cdhit;
my $formatdb = "$RealBin/blast/formatdb";
die "No such executable: $formatdb\n" unless -x $formatdb;

die "No such directory: $inDir\n" unless -d $inDir;
die "No such directory: $outDir\n" unless -d $outDir;
my @genomeFields = qw(fetch accession ncbi_assembly_name gtdb_taxonomy ncbi_taxonomy ncbi_strain_identifiers);
my @genomes = ReadTable($genomeFile, \@genomeFields);
print STDERR "Read " . scalar(@genomes) . " genomes\n";
my @fetched = ReadTable($fetchedFile, qw[fetch gid]);
my %fetchToGid = map { $_->{fetch} => $_->{gid} } @fetched;
@genomes = grep exists $fetchToGid{$_->{fetch}}, @genomes;
print STDERR "Genomes with fetched assemblies: " . scalar(@genomes) . "\n";

die "No genomes to load\n" if @genomes == 0;

my %orderToGenomes = ();
foreach my $genome (@genomes) {
  my $gtdbTax = parseTaxString($genome->{gtdb_taxonomy});
  die "Cannot parse taxonomy from " . $genome->{gtdb_taxonomy}
    unless defined $gtdbTax;
  my $order = $gtdbTax->{o};
  die unless defined $order;
  die "No order in " . $genome->{gtdb_taxonomy} if $order eq "";
  push @{ $orderToGenomes{$order} }, $genome;
}

print STDERR "Found " . scalar(keys %orderToGenomes) . " orders\n";

my $tabFile = "$outDir/SubDbs.tsv";
open(my $fhTab, ">", $tabFile) || die "Cannot write to $tabFile";
print $fhTab join("\t", qw{taxon level prefix nGenomes nProteins nClusters})."\n";

foreach my $order (sort keys %orderToGenomes) {
  my $orderGenomes = $orderToGenomes{$order};
  next if @$orderGenomes < 2;
  my $orderDir = "$outDir/" . orderToSubName($order);
  if (! -d $orderDir) {
    mkdir($orderDir) || die "mkdir $orderDir failed: $!\n";
  }

  # Run build.pl in test mode -- this will make the subdb and the fasta file
  my $tmp = $ENV{TMPDIR} || "/tmp";
  my $tmpPre = "$tmp/buildSubDbs.$$";
  my $tmpGenomes = "$tmpPre.genomes";
  my $tmpFetched = "$tmpPre.genomes.fetched";
  open(my $fhGenomes, ">", $tmpGenomes) || die "Cannot write to $tmpGenomes\n";
  open(my $fhFetched, ">", $tmpFetched) || die "Cannot write to $tmpFetched\n";
  print $fhGenomes join("\t", @genomeFields)."\n";
  print $fhFetched join("\t", qw{fetch gid})."\n";
  foreach my $genome (@$orderGenomes) {
    print $fhGenomes join("\t", map $genome->{$_}, @genomeFields)."\n";
    print $fhFetched join("\t", $genome->{fetch}, $fetchToGid{$genome->{fetch}})."\n";
  }
  close($fhGenomes) || die "Error writing to $tmpGenomes\n";
  close($fhFetched) || die "Error writing to $tmpFetched\n";
  print STDERR "Parsing " . scalar(@$orderGenomes) . " genomes for $order\n";
  my $buildCmd = "$RealBin/build.pl -genomes $tmpGenomes -in $inDir -out $orderDir -test -quiet";
  system($buildCmd) == 0 || die "Error for $order :\n$buildCmd\n-- $!\n";
  unlink($tmpGenomes);
  unlink($tmpFetched);


  die unless -e "$orderDir/neighbor.faa";
  my $proteinFaa = "$orderDir/sub.faa";
  rename("$orderDir/neighbor.faa", $proteinFaa) || die "Rename to $proteinFaa failed";
  my $clusterPre = "$orderDir/cluster.faa";
  print STDERR "Clustering $order\n";
  my $clusters = cluster('cdhit' => $cdhit, 'faa' => $proteinFaa, 'out' => $clusterPre,
                         'minIdentity' => $minIdentity, 'minCoverage' => $minCoverage,
                         'nCPUs' => $nCPUs);

  my $cInfoTab = "$orderDir/neighbor_ClusteringInfo.tab";
  open(my $fhCI, ">", $cInfoTab) || die "Cannot write to $cInfoTab\n";
  my $nProteins = 0;
  foreach my $cluster (values %$clusters) {
    $nProteins += scalar(@$cluster);
  }
  my $nClusters = scalar(keys %$clusters);
  my $nClusteredAA = proteinDbSize($clusterPre);
  print $fhCI join("\t", $nProteins, $nClusters, $nClusteredAA)."\n";
  close($fhCI) || die "Error writing to $cInfoTab\n";
  unlink("$clusterPre.clstr");
  formatBLASTp($formatdb, $clusterPre);

  # Build the sql database
  print STDERR "Building the sqlite3 database $orderDir/sub.db\n";
  my $dbFile = "$orderDir/sub.db";
  unlink($dbFile); # no old data
  my $sqlFile = "$RealBin/../lib/neighbor.sql";
  system("sqlite3 $dbFile < $sqlFile") == 0 || die $!;
  open(my $fhSql, "|-", "sqlite3", "$dbFile") || die "Cannot run sqlite3 on $dbFile";
  print $fhSql ".mode tabs\n";
  foreach my $table (qw{Genome Gene Protein ClusteringInfo Taxon}) {
    print $fhSql ".import $orderDir/neighbor_${table}.tab $table\n";
  }
  # Remove higher level taxa from Taxon
  print $fhSql qq{DELETE FROM Taxon WHERE level NOT IN ("order","family","genus","species");\n};
  close($fhSql) || die "Error loading database: $!";
  # Add the clusters ot the database
  clusteringIntoDb($clusters, $dbFile);
  foreach my $table (qw{Genome Gene Protein ClusterProtein ClusteringInfo Taxon}) {
    unlink("$orderDir/neighbor_${table}.tab");
  }
  system("gzip --force $clusterPre") == 0 || die "gzip failed on $clusterPre";
  system("gzip --force $proteinFaa") == 0 || die "gzip failed on $proteinFaa";
  print $fhTab join("\t", $order, "order", orderToSubName($order),
                    scalar(@$orderGenomes), $nProteins, $nClusters)."\n";
}
close($fhTab) || die "Error writing to $tabFile";
print STDERR "Wrote $tabFile\n";
print STDERR <<END
To load it into neighbor.db use
.mode tabs
DELETE * from SubDb;
.import $tabFile SubDb
END
  ;
print STDERR "buildSubDbs.pl done\n";

