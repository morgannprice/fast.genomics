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

my $nSlices = 8;
my $usage = <<END
build.pl -genomes genomes.tsv -fetched assemblies.tsv -dir dir -test

Creates the neighbor database, including the mysql database
(neighbor.db), a big fasta file (neighbor.faa.gz), and a sliced mmseqs
index (neighbor.sliced*). The genomes table must include the fields
fetch, (gtdb) accession, ncbi_assembly_name, gtdb_taxonomy,
ncbi_taxonomy, and ncbi_strain_identifiers. The fetched table must
include the fields fetch and gid.

Optional arguments:
-out dir -- where to put the output files
-test -- test mode: build tab-delimited tables and faa file, but
   do not delete them or build the mmseqs or sqlite3 databases.
-slices $nSlices -- how many slices to use for the mmseqs db
END
;

# sqlite3 expects CVS format, not exactly tab delimited format
# So, need to replace any " with "" and surround the field with quotes.
sub csv_quote($);

my ($genomeFile, $fetchedFile, $inDir, $outDir, $test);
die $usage
  unless GetOptions('genomes=s' => \$genomeFile,
                    'fetched=s' => \$fetchedFile,
                    'dir=s' => \$inDir,
                    'out=s' => \$outDir,
                    'test' => \$test,
                    'slices=i' => \$nSlices)
  && defined $genomeFile
  && defined $fetchedFile
  && defined $inDir
  && @ARGV == 0;

die "No such directory: $inDir\n" unless -d $inDir;
$outDir = $inDir unless defined $outDir;
die "No such directory: $outDir\n" unless -d $outDir;
die "nSlices is out of range\n" unless $nSlices >= 1 && $nSlices <= 100;

my @genomes = ReadTable($genomeFile,
  qw[fetch accession ncbi_assembly_name gtdb_taxonomy ncbi_taxonomy ncbi_strain_identifiers]);
print STDERR "Read " . scalar(@genomes) . " genomes\n";
my @assemblyRows = ReadTable($fetchedFile, qw[fetch gid]);
my %fetchToGid = map { $_->{fetch} => $_->{gid} } @assemblyRows;
@genomes = grep exists $fetchToGid{$_->{fetch}}, @genomes;
print STDERR "Genomes with fetched assemblies: " . scalar(@genomes) . "\n";

die "No genomes to load\n" if @genomes == 0;

my $dbFile = "$outDir/neighbor.db";
my $slicedDb = "$outDir/neighbor.sliced";
if (!defined $test) {
  unlink($dbFile);
  my $sqlFile = "$RealBin/../lib/neighbor.sql";
  system("sqlite3 $dbFile < $sqlFile") == 0 || die $!;
  print STDERR "Created empty database $dbFile\n";
}
print STDERR "Reading genomes\n";

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

  # The feature file usually has pairs of rows, first a gene entry, then usually a CDS
  # (even if a pseudogene), or rRNA or ncRNA type entry, which will have the description under "name".
  # But sometimes the rows are not in order (i.e., two rows for an antisense RNA between
  # the gene and CDS entries, see BSU_32469 in GCF_000009045.1). Or sometimes
  # there's no CDS entry at all.
  # So first, group the entries by locus_tag.
  my %locusTagFeatures = ();
  foreach my $feature (@$features) {
    push @{ $locusTagFeatures{ $feature->{locus_tag} } }, $feature;
  }
  my $nGenes = 0;
  my $nGenesWithProteins = 0;
  foreach my $locus_tag (sort keys %locusTagFeatures) {
    if ($locus_tag eq "") {
      print STDERR "Warning: row(s) with empty locus_tag in $featureFile\n";
      next;
    }
    my $list = $locusTagFeatures{$locus_tag};
    my @genes = grep $_->{"# feature"} eq "gene", @$list;
    my @others = grep $_->{"# feature"} ne "gene", @$list;
    if (@genes == 0) {
      print STDERR "Warning: no gene entry for locus_tag $locus_tag in $featureFile\n";
      next;
    }
    print STDERR "Warning: more than one gene entry for locus_tag $locus_tag in $featureFile\n"
      if @genes > 1;
    print STDERR "Warning: more than one non-gene entry for locus_tag $locus_tag in $featureFile\n"
      if @others > 1;
    my $gene = $genes[0];
    my $other = $others[0]; # or undef

    my $locusTag = $gene->{"locus_tag"};
    my $proteinId = "";
    my $desc = "?";
    if (defined $other) {
      $desc = $other->{name};
      $desc = $desc . " (pseudogene)" if $other->{class} eq "without_protein";
    } else {
      $desc = "pseudogene" if $gene->{attributes} =~ m/pseudo/i;
    }
    # The CDS may have class="with_protein" or empty
    if ($gene->{class} eq "protein_coding" && defined $other
        && $other->{"# feature"} eq "CDS"
        && ($other->{class} eq "with_protein" || $other->{class} eq "")) {
      $proteinId = $other->{"product_accession"};
      if (!exists $protSeen{$proteinId}) {
        print STDERR "Warning: unknown protein id $proteinId in $featureFile for $locusTag\n";
        $proteinId = "";
      }
    }

    # Fields in the Gene table are gid, scaffoldId, begin, end, strand, locusTag, proteinId, desc
    print { $fh{Gene} } join("\t", $gid,
                             $gene->{genomic_accession}, $gene->{start}, $gene->{end}, $gene->{strand},
                             $locusTag, $proteinId, csv_quote($desc)) . "\n";
    $nGenes++;
    $nGenesWithProteins++ if $proteinId ne "";
  }
  print STDERR "Warning: $gid has $nGenes genes but only $nGenesWithProteins with protein sequences !\n"
    if $nGenesWithProteins <= $nGenes/2.0;

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
                             $row->{ncbi_assembly_name}, $row->{ncbi_taxonomy},
                            $nGenes, $nGenesWithProteins)."\n";
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
if (defined $test) {
  print STDERR "Finished with test mode\n";
  exit(0);
}

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

print STDERR "Creating and indexing the sliced mmseqs database\n";
my $cmd = "$RealBin/buildSliced.pl -in $faaOut -slices $nSlices -out $slicedDb";
system($cmd) == 0
  || die "Error running buildSliced.pl: $!\nCommand: $cmd\n";
print STDERR "Finished building:\nsqlite3 $dbFile\nsliceddb $slicedDb\n";

# Remove tab-delimited imports
foreach my $file (values %files) {
  unlink($file);
}

system("gzip --force $faaOut") == 0 || die "gzip failed";

sub csv_quote($) {
  my ($in) = @_;
  return $in unless $in =~ m/"/;
  $in =~ s/"/""/g;
  return '"' . $in . '"';
}
