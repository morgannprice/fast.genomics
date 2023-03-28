#!/usr/bin/perl -w
# Download a table of genomes
use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use lib "../lib";
# neighborWeb.pm relies on various PaperBLAST libraries
use lib "../../PaperBLAST/lib";
use neighborWeb;

# CGI arguments: order (optional)
setOrder(param('order'));

my $genomes = getDbHandle()->selectall_arrayref("SELECT * from Genome"
  . " ORDER BY gtdbDomain, gtdbPhylum, gtdbClass, gtdbOrder, gtdbFamily, gtdbGenus, gtdbSpecies, strain",
  { Slice => {} });
my $fileName = "genomes.tsv";
$fileName = getSubDb() . "_" . $fileName if getOrder() ne "";
print "Content-Type:text/tab-separated-values\n";
print "Content-Disposition: attachment; filename=$fileName\n\n";
my @fields = qw{assemblyId gtdbDomain gtdbPhylum gtdbClass gtdbOrder gtdbFamily gtdbGenus gtdbSpecies strain gtdbAccession assemblyName ncbiTaxonomy nGenes nProteins};
print join("\t", @fields)."\n";
foreach my $genome (@$genomes) {
  $genome->{assemblyId} = $genome->{gid};
  my @out = map $genome->{$_}, @fields;
  print join("\t", @out)."\n";
}
