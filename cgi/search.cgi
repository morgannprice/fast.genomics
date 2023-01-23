#!/usr/bin/perl -w
use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use IO::Handle; # for autoflush
use DBI;
use HTML::Entities;
use lib "../lib";
use neighbor;
# neighborWeb.pm relies on various PaperBLAST libraries
use lib "../../PaperBLAST/lib";
use neighborWeb;

# CGI arguments:
# query (optional) -- usually an identifier or locus tag,
#   or a genus and word(s) from a protein description.
#   See neighborWeb::parseGeneQuery for details.
my $cgi = CGI->new;
my $query = $cgi->param('query') || "";

# Redirect to locus tag if the query is an exact match to a
# locus tag or proteinId in our database
my $geneRedirect = locusTagToGene($query); # checks upper case as well
if (! $geneRedirect && $query =~ m/^[a-zA-Z0-9_.]+$/) {
  my $genesFromProteinId = getDbHandle()->selectall_arrayref(
    qq{ SELECT * from Gene WHERE proteinId = ? },
    { Slice => {} }, uc($query));
  $geneRedirect = $genesFromProteinId->[0] if @$genesFromProteinId == 1;
}
if ($geneRedirect) {
  print redirect(-url => "gene.cgi?locus=" . $geneRedirect->{locusTag});
  exit(0);
}

start_page('title' => $query eq "" ? "" : "Gene search");
autoflush STDOUT 1; # show preliminary results

my %query = parseGeneQuery($query);
if (!defined $query{gene} && !defined $query{genes} && !defined $query{seq}) {
  print p(b($query{error})) if $query{error};
  print start_form( -name => 'input', -method => 'GET', -action => 'search.cgi' ),
    p("Enter an identifier from UniProt, RefSeq, PDB, or MicrobesOnline,",
      br(),
      "or a protein sequence in FASTA or Uniprot format,",
      br(),
      "or a genus name and a protein description",
      br(),
      textarea( -name  => 'query', -value => '', -cols  => 70, -rows  => 10 ),
      br(),
      br(),
      submit('Search'), reset()),
    end_form;
  finish_page();
}

if (defined $query{genes}) {
  my $genes = $query{genes};
  # show table of genes
  my $gid = $genes->[0]{gid};
  my $genome = getDbHandle()->selectrow_hashref("SELECT * FROM Genome WHERE gid = ?", {}, $gid);
  die $gid unless $genome;
  print p("Found", scalar(@$genes), " matches in $genome->{gtdbSpecies} $genome->{strain} ($genome->{gid})");
  foreach my $gene (@$genes) {
    my $locusTag = $gene->{locusTag};
    print p(a({href => "gene.cgi?locus=$locusTag"}, $locusTag),
            $gene->{proteinId}, $gene->{desc});
  }
  finish_page();
}

if (defined $query{gene}) {
  my $gene = $query{gene};
  my $proteinId = $gene->{proteinId};
  my $genome = getDbHandle()->selectrow_hashref("SELECT * FROM Genome WHERE gid = ?", {}, $gene->{gid});
  print p("Found gene", $gene->{locusTag}, $gene->{proteinId}, encode_entities($gene->{desc}),
          "in",
          i($genome->{gtdbSpecies}), $genome->{strain});
} elsif (defined $query{seqDesc}) {
  print p("Found $query{seqDesc}")."\n";
}

finish_page() unless defined $query{seq};
my $seq = $query{seq};
$seq =~ s/[*]//g;
$seq =~ m/^[a-zA-Z]+$/ || die("Sorry, input sequence has invalid characters");
my $seqDesc = $query{seqDesc}
  || length($seq) . " a.a. beginning with " . substr($seq, 0, 10);
print showSequence($seqDesc, $seq), "\n";

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
        scalar(keys %gidHits) . " of $nGenomes genomes");
print p("Found",
        scalar(grep $_->{bits} >= $maxScore*0.3, @$hitGenes),
        "hits with &ge;30% of max score (at least " . sprintf("%.1f",0.3 * $maxScore) . " bits)");

finish_page();
