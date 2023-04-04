#!/usr/bin/perl -w
use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use HTML::Entities;
use lib "../lib";
use neighbor;
# neighborWeb.pm relies on various PaperBLAST libraries
use lib "../../PaperBLAST/lib";
use neighborWeb;
use pbweb qw{commify};

# required CGI arguments:
# locus (a locus tag in the database) or seqDesc and seq
# optional arguments:
# order (which subdb to use)
#
# Forwards to the gene page if locus is set

my $cgi = CGI->new;
setOrder(param('order'));
my ($gene, $seqDesc, $seq) = getGeneSeqDesc($cgi);

if (defined $gene) {
  # should not be reachable, but redirect to gene page just in case
  print redirect(-url => addOrderToURL("gene.cgi?locus=$gene->{locusTag}"));
  exit(0);
}

my $title = encode_entities($seqDesc);
start_page('title' => $title);

print
h3("Protein sequence"),
    formatFastaHtml($seqDesc, $seq);

my $options = geneSeqDescSeqOptions($gene, $seqDesc, $seq);
if (hasHits($seq)) {
  print p("See",
          join(", or ",
               a({-href => "neighbors.cgi?$options"}, "gene neighborhoods"),
               a({-href => "hitTaxa.cgi?$options"}, "taxonomic distribution")
               . " of its homologs",
               a({-href => "downloadHomologs.cgi?$options",
                  -title => "tab-delimited table of homologs"},
                 "download homologs"),
               a({-href => "compare.cgi?$options",
                  -title => "compare presence/absence of homologs and their proximity"},
                 "compare presence/absence")));
} else {
  print p(a({-href => "findHomologs.cgi?$options"}, "Find homologs with",
            getOrder() eq "" ? "mmseqs2" : "clustered BLAST"),
          "(fast)");
}
if (getOrder() ne "") {
  print p("Or find",
          a({-href => "findHomologs.cgi?seq=$seq&seqDesc=".encode_entities($seqDesc)},
            "homologs in diverse bacteria and archaea"));
} elsif (hasHits($seq)) {
  # Find the order of the top hit and suggest searching there
  my $order = getTopHitOrder($seq);
  if (defined $order) {
    my $nSubGenomes = moreGenomesInSubDb("order", $order, $order);
    print p("Or find",
            a({-href => "findHomologs.cgi?order=$order&$options"},
              "homologs in", commify($nSubGenomes), $order))
      if $nSubGenomes > 0;
  }
}

print
  h3("Other sequence analysis tools"),
  start_ul,
  map li($_), proteinAnalysisLinks($seqDesc, $seq, undef);
print end_ul;
finish_page();
