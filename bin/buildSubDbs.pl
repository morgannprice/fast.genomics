#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use IO::Handle; # for autoflush
use lib "$RealBin/../lib";
use neighbor qw{parseTaxString featuresToGenes orderToSubName};
use lib "$RealBin/../../PaperBLAST/lib";
use FetchAssembly qw{ParseNCBIFeatureFile};
use pbutils qw{ReadFastaEntry ReadTable};
my $nCPUs = $ENV{MC_CORES} || (`egrep -c '^processor' /proc/cpuinfo`);

my $minCoverage = 0.8;
my $minIdentity = 0.8;
my $nSlices = 4;

my $usage = <<END
./buildSubDbs.pl -genomes genomes.tsv [ -fetched genomes.tsv.fetched ] -in indir -out outdir

Creates the sub databases for each order, in subdirectories of the
output directory.

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

my $iCluster = 0; # global cluster number
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

  # Run cd-hit
  my $proteinFaa = "$orderDir/neighbor.faa";
  die unless -e $proteinFaa;
  my $clusterPre = "$orderDir/cluster.faa";
  # cluster at 80% identity, no memory limit, use 5-mers
  # -d 0 means use sequence name in fasta header up to first white space
  # -aS and -aL set the coverage of the sequence and the cluster representative
  print STDERR "Clustering $order\n";
  my $clusterCmd = "$cdhit -i $proteinFaa -o $clusterPre -c $minIdentity -M 0 -T $nCPUs -n 5 -d 0"
    . " -aS $minCoverage -aL $minCoverage >& $clusterPre.log";
  system($clusterCmd) == 0 || die "cd-hit failed:\n$clusterCmd\n$!\n";

  # Parse the .clstr file
  my @clusters = ();
  open (my $fhClust, "<", "$clusterPre.clstr")
    || die "Cannot read $clusterPre.clstr";
  while (my $line = <$fhClust>) {
    if ($line =~ m/^>/) {
      push @clusters, []; # start an emtpy cluster
    } else {
      $line =~ m/>(\S+)[.][.][.] / || die "Cannot parse cluster line\n$line\n";
      my $proteinId = $1;
      die "Cluster member before cluster number\n" if @clusters == 0;
      push @{ $clusters[-1] }, $proteinId;
    }
  }
  close($fhClust) || die "Error reading $clusterPre.clstr";
  my $clusterTab = "$orderDir/neighbor_ClusterProtein.tab";
  open (my $fhCT, ">", $clusterTab)
    || die "Cannot write to $clusterTab\n";
  foreach my $cluster (@clusters) {
    $iCluster++;
    foreach my $proteinId (@$cluster) {
      print $fhCT join("\t", $iCluster, $proteinId)."\n";
    }
  }
  close($fhCT) || die "Error writing to $clusterTab\n";

  my $cInfoTab = "$orderDir/neighbor_ClusteringInfo.tab";
  open(my $fhCI, ">", $cInfoTab) || die "Cannot write to $cInfoTab\n";
  my $nProteins = 0;
  foreach my $cluster (@clusters) {
    $nProteins += scalar(@$cluster);
  }
  my $nClusters = scalar(@clusters);
  print $fhCI "$nProteins\t$nClusters\n";
  close($fhCI) || die "Error writing to $cInfoTab\n";

  unlink("$clusterPre.clstr");

  my $formatCmd = "$formatdb -p T -i $clusterPre";
  system($formatCmd) == 0 || die "formatdb failed:\n$formatCmd\n$!";

  # Build the sql database
  print STDERR "Building the sqlite3 database $orderDir/sub.db\n";
  my $dbFile = "$orderDir/sub.db";
  unlink($dbFile);
  my $sqlFile = "$RealBin/../lib/neighbor.sql";
  system("sqlite3 $dbFile < $sqlFile") == 0 || die $!;
  open(SQLITE, "|-", "sqlite3", "$dbFile") || die "Cannot run sqlite3 on $dbFile";
  autoflush SQLITE 1;
  print SQLITE ".mode tabs\n";
  foreach my $table (qw{Genome Gene Protein ClusterProtein ClusteringInfo Taxon}) {
    print SQLITE ".import $orderDir/neighbor_${table}.tab $table\n";
  }
  # Remove higher level taxa from Taxon
  print SQLITE qq{DELETE FROM Taxon WHERE level NOT IN ("order","family","genus","species");\n};
  close(SQLITE) || die "Error loading database: $!";
  foreach my $table (qw{Genome Gene Protein ClusterProtein ClusteringInfo Taxon}) {
    unlink("$orderDir/neighbor_${table}.tab");
  }
  system("gzip --force $clusterPre") == 0 || die "gzip failed on $clusterPre";
  system("gzip --force $proteinFaa") == 0 || die "gzip failed on $proteinFaa";
}
print STDERR "buildSubDbs.pl done\n";


