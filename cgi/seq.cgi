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

# CGI arguments:
# locus (a locus tag in the database) or seqDesc and seq
#
# Forwards to the gene page if locus is set

my $cgi = CGI->new;
my ($gene, $seqDesc, $seq) = getGeneSeqDesc($cgi);

if (defined $gene) {
  # should not be reachable, but redirect to gene page just in case
  print redirect(-url => "gene.cgi?locus=$gene->{locusTag}");
  exit(0);
}

my $title = encode_entities($seqDesc);
start_page('title' => $title);

print
h3("Protein sequence"),
    formatFastaHtml($seqDesc, $seq);

my $options = geneSeqDescSeqOptions($gene, $seqDesc, $seq);
if (hasMMSeqsHits($seq)) {
  print p("See",
          a({-href => "neighbors.cgi?$options"}, "gene neighborhoods"),
          "or",
          a({-href => "hitTaxa.cgi?$options"}, "taxonomic distribution"),
          "of its homologs",
          "or",
          a({-href => "downloadHomologs.cgi?$options",
             -title => "tab-delimited table of homologs"},
            "download"));
} else {
  print p(a({-href => "findHomologs.cgi?$options"}, "Find homologs with mmseqs2"),
          "(fast)");
}
print
  h3("Other sequence analysis tools"),
  start_ul,
  map li($_), proteinAnalysisLinks($seqDesc, $seq, undef);
print end_ul;
print end_html;
