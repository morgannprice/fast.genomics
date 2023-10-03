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

# Required CGI arguments:
# locus (a locus tag in the database), or seqDesc and seq
# Optional arguments:
# order (which subdb to use)
# compare1 -- for linking "back" to compare.cgi
#   compare1 is the compare.cgi parameters for the first gene or sequence;
#   this locus/seqDesc/seq will be the second gene or sequence

my $cgi = CGI->new;
setOrder(param('order'));
my ($gene, $seqDesc, $seq) = getGeneSeqDesc($cgi);

my $title = getOrder() eq "" ?
  "Find homologs in diverse bacteria and archaea"
  : "Find homologs in " . getOrder();

start_page('title' => $title);
autoflush STDOUT 1; # show preliminary results

if (defined $gene && !defined $seq) {
  # This shouldn't be reached, but just in case
  print p("Not a protein-coding gene:", encode_entities($seqDesc)),
    p(a({ -href => addOrderToURL("gene.cgi?locus=".$gene->{locusTag}) }, "see gene"));
  finish_page();
}
#else
die "No sequence input\n" unless defined $seq && $seq =~ m/^[A-Z]+$/;

print p("Finding homologs for", encode_entities($seqDesc)),
  showSequence($seqDesc, $seq),
  "\n";
my $hits = getHits(uc($seq));

my $maxScore = 0;
if (@$hits > 0) {
  my $top = $hits->[0];
  my $coverage = ($top->{qEnd} - $top->{qBegin} + 1)/length($seq);
  $maxScore = $top->{bits} / ($top->{identity} * $coverage);
}
print p("Found at least", commify(scalar(@$hits)), "hits");
print p("Found at least",
        commify(scalar(grep $_->{bits} >= $maxScore*0.3, @$hits)),
        "hits with &ge;30% of max score (at least " . sprintf("%.1f",0.3 * $maxScore) . " bits)");

my $options = geneSeqDescSeqOptions($gene,$seqDesc,$seq);
print p("See",
        join(", or ",
             a({-href => addOrderToURL("neighbors.cgi?${options}")}, "gene neighborhoods"),
             a({-href => addOrderToURL("hitTaxa.cgi?${options}")}, "taxonomic distribution")
             . " of homologs",
             a({-href => addOrderToURL("downloadHomologs.cgi?${options}"),
                -title => "tab-delimited table of homologs"}, "download homologs"),
             a({-href => addOrderToURL("compare.cgi?${options}")}, "compare presence/absence")))
  if @$hits > 0;
print p("Or see", a({-href => addOrderToURL("gene.cgi?locus=$gene->{locusTag}")}, "gene"))
  if defined $gene;

if (getOrder() eq "") {
  my $order = getTopHitOrder($hits);
  if (defined $order) {
    my $nSubGenomes = moreGenomesInSubDb("order", $order, $order);
    print p("Or find",
            a({ -href => "findHomologs.cgi?order=$order&${options}" },
              "homologs in", commify($nSubGenomes), $order))
      if $nSubGenomes > 0;
  }
} else {
  print p("Or find",
          a({-href => "findHomologs.cgi?" . geneSeqDescSeqOptionsNoOrder($gene,$seqDesc,$seq) },
            "homologs in diverse bacteria and archaea"));
}

if (param('compare1')) {
  my $options1 = param('compare1');
  my $options2 = defined $gene ? "locus2=".$gene->{locusTag}
    : "seqDesc2=" . uri_escape($seqDesc) . "&seq2=" . $seq;
  print p("Back to",
          a({-href => addOrderToURL("compare.cgi?${options1}&${options2}")},
            "compare gene presence/absence")), "\n";
}

my $hitsFile = hitsFile($seq);
print "<!-- hits are in $hitsFile -->\n"; # aids debugging

print
  h3("Other sequence analysis tools"),
  start_ul,
  map li($_), proteinAnalysisLinks($seqDesc, $seq,
                                   defined $gene ? gidToGenome($gene->{gid}) : undef);
print end_ul;

finish_page();

