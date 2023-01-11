#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use IO::Handle;                 # for autoflush
use lib "$RealBin/../lib";
use neighbor qw{parseTaxString};
use lib "$RealBin/../../PaperBLAST/lib";
use FetchAssembly qw{ParseNCBIFeatureFile};
use pbutils qw{ReadFastaEntry ReadTable};

my $usage = <<END
build.pl -genomes genomes.tsv -fetched assemblies.tsv -dir dir

Creates the neighbor database, including the mysql database
(neighbor.db), a big fasta file (neighbor.faa), and an mmseqs index
(neigbor.mmseqs). The genomes table must include the fields fetch,
(gtdb) accession, ncbi_assembly_name, gtdb_taxonomy, ncbi_taxonomy,
and ncbi_strain_identifiers. The fetched table must include the fields
fetch and gid.

Optional arguments:
-out dir -- where to put the output files
-tmp dir/tmp -- the temporary directory for mmseqs
-mmseqs mmseqs_executable -- by default, assumes mmseqs is on the path
END
;

# sqlite3 expects CVS format, not exactly tab delimited format
# So, need to replace any " with "" and surround the field with quotes.
sub csv_quote($);

my ($genomeFile, $fetchedFile, $inDir, $outDir, $tmpDir);
my $mmseqs = "mmseqs";
die $usage
  unless GetOptions('genomes=s' => \$genomeFile,
                    'fetched=s' => \$fetchedFile,
                    'dir=s' => \$inDir,
                    'out=s' => \$outDir,
                    'tmp=s' => \$tmpDir,
                    'mmseqs=s' => \$mmseqs)
  && defined $genomeFile
  && defined $fetchedFile
  && defined $inDir
  && @ARGV == 0;

die "No such directory: $inDir\n" unless -d $inDir;
$outDir = $inDir unless defined $outDir;
die "No such directory: $outDir\n" unless -d $outDir;
if (!defined $tmpDir) {
  $tmpDir = "$inDir/tmp";
  mkdir($tmpDir) unless -d $tmpDir;
}
die "No such directory: $tmpDir\n" unless -d $tmpDir;

system("$mmseqs -h > /dev/null") == 0
  || die "Cannot run $mmseqs\n";

my @genomes = ReadTable($genomeFile,
  qw[fetch accession ncbi_assembly_name gtdb_taxonomy ncbi_taxonomy ncbi_strain_identifiers]);
print STDERR "Read " . scalar(@genomes) . " genomes\n";
my @assemblyRows = ReadTable($fetchedFile, qw[fetch gid]);
my %fetchToGid = map { $_->{fetch} => $_->{gid} } @assemblyRows;
@genomes = grep exists $fetchToGid{$_->{fetch}}, @genomes;
print STDERR "Genomes with fetched assemblies: " . scalar(@genomes) . "\n";

die "No genomes to load\n" if @genomes == 0;

my $dbFile = "$outDir/neighbor.db";
my $mmDb = "$outDir/neighbor.mmseqs";
unlink($dbFile);
unlink($mmDb);

my $sqlFile = "$RealBin/../lib/neighbor.sql";
system("sqlite3 $dbFile < $sqlFile") == 0 || die $!;
print STDERR "Created empty database $dbFile; reading genomes.\n";

# As we read the genomes, write to the tables
my @tables = qw{Genome Gene Protein Taxon};
my %files = map { $_ => "$outDir/neighbor_" . $_ . ".tab" } @tables;
my %fh = ();
foreach my $table (@tables) {
  open($fh{$table}, ">", $files{$table})
    || die "Cannot write to $files{$table}";
}

my $faaOut = "$outDir/neighbor.faa";
open(my $fhFaaOut, ">", $faaOut) || die "Cannot write to $faaOut\n";

# The same protein might be in multiple genomes
my %protSeen = ();

# For counting across taxa
my @levels = qw{domain phylum class order family genus species};
my %charToLevel = map { substr($_, 0, 1) => $_ } @levels;
my %nGenomes = (); # level => taxon => #genomes
my %taxParent = (); # level => taxon => parent

foreach my $row (@genomes) {
  my $fetch = $row->{fetch};
  my $gid = $fetchToGid{$fetch};
  my $faaIn = "$inDir/refseq_${gid}.faa";
  die "No faa file for $gid: $faaIn\n" unless -e $faaIn;

  my $featureFile = "$inDir/refseq_${gid}.features.tab";
  die "No feature table for $gid: $featureFile\n" unless -e $featureFile;
  my $features = ParseNCBIFeatureFile($featureFile);

  open(my $fhIn, "<", $faaIn) || die "Cannot read $faaIn";
  my $state = {};
  while (my ($header,$seq) = ReadFastaEntry($fhIn, $state)) {
    $header =~ s/ .*//;
    my $proteinId = $header;
    next if exists $protSeen{$proteinId};
    # If using a filehandle in a hash, need to enclode it in a block
    print { $fh{Protein} } "$proteinId\t$seq\n";
    print $fhFaaOut ">$proteinId\n$seq\n";
    $protSeen{$proteinId} = 1;
  }
  close($fhIn) || die "Error reading $faaIn";

  # Read the feature file as pairs of rows -- first a gene entry, then a CDS (even if a pseudogene)
  # or rRNA or ncRNA type entry, which will have the description under "name".
  while(@$features > 0) {
    my $gene = shift @$features;
    die "Error parsing feature table $featureFile" unless $gene->{"# feature"} eq "gene";
    my $second = shift @$features;
    die "Error parsing feature table $featureFile" 
      unless defined $second && $second->{"# feature"} ne "gene";

    my $locusTag = $gene->{"locus_tag"};
    my $proteinId = "";
    my $desc = $second->{name};
    $desc = $desc . " (pseudogene)" if $second->{class} eq "without_protein";

    if ($gene->{class} eq "protein_coding" && $second->{class} eq "with_protein") {
      $proteinId = $second->{"product_accession"};
      if (!exists $protSeen{$proteinId}) {
        print STDERR "Warning: unknown protein id $proteinId in $featureFile for $locusTag\n";
        $proteinId = "";
      }
    }

    # Fields in the Gene table are gid, scaffoldId, begin, end, strand, locusTag, proteinId, desc
    print { $fh{Gene} } join("\t", $gid,
                             $gene->{genomic_accession}, $gene->{start}, $gene->{end}, $gene->{strand},
                             $locusTag, $proteinId, csv_quote($desc)) . "\n";

  }
  # Fields in the Genome table are gid
  # gtdbDomain Phylum Class Order Family Genus Species
  # strain gtdbAccession assemblyName ncbiTaxonomy
  my $gtdbTax = parseTaxString($row->{gtdb_taxonomy});
  if (!defined $gtdbTax) {
    print STDERR "Warning: cannot parse taxonomy string $row->{gtdb_taxonomy}\n";
    $gtdbTax = {};
  }
  my @taxa = map $gtdbTax->{$_} || "", qw{d p c o f g s};
  # Increase the count for each level
  foreach my $levelChar (keys %charToLevel) {
    my $level = $charToLevel{$levelChar};
    $nGenomes{$level}{ $gtdbTax->{$levelChar} }++
      if exists $gtdbTax->{$levelChar};
  }
  # Index by domain/phylum/etc. instead of d/p/etc.
  my %gtdbLevels = map { $charToLevel{$_} => $gtdbTax->{$_} } keys(%$gtdbTax);
  foreach my $i (0..5) {
    my $parentLevel = $levels[$i];
    my $level = $levels[$i+1];
    if (exists $gtdbLevels{$level} && exists $gtdbLevels{$parentLevel}) {
      $taxParent{$level}{ $gtdbLevels{$level} } = $gtdbLevels{$parentLevel};
    }
  }
  print { $fh{Genome} } join("\t", $gid, @taxa,
                             csv_quote($row->{ncbi_strain_identifiers}), $row->{accession},
                             $row->{ncbi_assembly_name}, $row->{ncbi_taxonomy})."\n";
}

foreach my $level (@levels) {
  foreach my $taxon (sort keys %{ $nGenomes{$level} }) {
    print { $fh{Taxon} } join("\t", csv_quote($taxon), $level,
                              csv_quote($taxParent{$level}{$taxon} || ""),
                              $nGenomes{$level}{$taxon})."\n";
  }
}

close($fhFaaOut) || die "Error writing to $faaOut\n";
foreach my $table (@tables) {
  close($fh{$table}) || die "Error writing to $files{$table}";
}

print STDERR "Read all genomes and built neighbor.faa and tables for loading\n";

open(SQLITE, "|-", "sqlite3", "$dbFile") || die "Cannot run sqlite3 on $dbFile";
autoflush SQLITE 1;
print SQLITE ".mode tabs\n";
foreach my $table (@tables) {
  print STDERR "Loading table $table\n";
  print SQLITE ".import $files{$table} $table\n";
}
print SQLITE <<END
SELECT 'nGenomes', COUNT(*) FROM Genome;
SELECT 'nGenes', COUNT(*) FROM Gene;
SELECT 'nProteinGenes', COUNT(*) FROM Gene WHERE proteinId <> "";
SELECT 'nProteins', COUNT(*) FROM Protein;
SELECT 'ProteinLength', sum(length(sequence))/1e6, 'millions' FROM Protein;
SELECT 'nGenus', COUNT(*) FROM Taxon WHERE level="genus";
END
;
close(SQLITE) || die "Error running sqlite3 import commands\n";

print STDERR "Creating and indexing mmseqs database, see $outDir/mmseqs.log\n";
system("(mmseqs createdb $faaOut $mmDb; mmseqs createindex $mmDb $tmpDir) >& $outDir/mmseqs.log") == 0
  || die "Error running mmseqs: $!\n";
print STDERR "Finished building:\nsqlite3\t$dbFile\nmmseqs2\t$mmDb\n";

# Remove tab-delimited imports
foreach my $file (values %files) {
  unlink($file);
}

sub csv_quote($) {
  my ($in) = @_;
  return $in unless $in =~ m/"/;
  $in =~ s/"/""/g;
  return '"' . $in . '"';
}
