#!/usr/bin/perl -w
use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use IO::Handle; # for autoflush
use DBI;
use HTML::Entities;
use URI::Escape;
use lib "../lib";
use neighbor;
# neighborWeb.pm relies on various PaperBLAST libraries
use lib "../../PaperBLAST/lib";
use pbweb qw{commify};
use neighborWeb;

# CGI arguments:
# locus (a locus tag in the database) or seqDesc and seq
my $cgi = CGI->new;
my ($gene, $seqDesc, $seq) = getGeneSeqDesc($cgi);

start_page('title' => 'Find homologs');
autoflush STDOUT 1; # show preliminary results

if (defined $gene && !defined $seq) {
  # This shouldn't be reached, but just in case
  print p("Not a protein-coding gene:", encode_entities($seqDesc)),
    p(a({ -href => "gene.cgi?locus=".$gene->{locusTag} }, "see gene"));
  finish_page();
}
#else
die "No sequence input\n" unless defined $seq && $seq =~ m/^[A-Z]+$/;

print p("Finding homologs for", encode_entities($seqDesc)),
  showSequence($seqDesc, $seq),
  "\n";
my $hits = getMMSeqsHits(uc($seq));
my $hitGenes = hitsToGenes($hits);

my $maxScore = 0;
if (@$hitGenes > 0) {
  my $top = $hitGenes->[0];
  my $coverage = ($top->{qEnd} - $top->{qBegin} + 1)/length($seq);
  $maxScore = $top->{bits} / ($top->{identity} * $coverage);
}
my %gidHits = map { $_->{gid} => $_ } @$hitGenes;
my ($nGenomes) = getDbHandle()->selectrow_array("SELECT COUNT(*) FROM Genome");
print p("Found", scalar(@$hitGenes), "hits in",
        scalar(keys %gidHits) . " of", commify($nGenomes), "genomes");
print p("Found",
        scalar(grep $_->{bits} >= $maxScore*0.3, @$hitGenes),
        "hits with &ge;30% of max score (at least " . sprintf("%.1f",0.3 * $maxScore) . " bits)");

my $options = defined $gene ? "locus=".$gene->{locusTag}
  : "seqDesc=" . encode_entities($seqDesc) . "&" . "seq=$seq";
print p("See", a({-href => "neighbors.cgi?${options}"}, "gene neighborhoods"), "of homologs")
  if @$hitGenes > 0;
print p("Or see", a({-href => "gene.cgi?locus=$gene->{locusTag}"}, "gene"))
  if defined $gene;

finish_page();
