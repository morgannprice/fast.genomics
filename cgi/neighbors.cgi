#!/usr/bin/perl -w
use strict;

use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use IO::Handle; # for autoflush
use DBI;
use HTML::Entities;
use URI::Escape;
use List::Util qw{min max};
use lib "../lib";
use neighbor;
use genesSvg;
# neighborWeb.pm relies on various PaperBLAST libraries
use lib "../../PaperBLAST/lib";
use neighborWeb;

# CGI arguments:
# locus (a locus tag in the database) or seqDesc and seq
# n -- max number of hits to show
my $cgi = CGI->new;
my ($gene, $seqDesc, $seq) = getGeneSeqDesc($cgi);
my $n = $cgi->param('n');
$n = 50 unless defined $n && $n =~ m/^\d+$/;
$n = 200 if $n > 200;
my $ntShown = 6000;

if (defined $gene && ! $seq) {
  # should not be reachable, but redirect to gene page just in case
  print redirect(-url => "gene.cgi?locus=$gene->{locusTag}");
  exit(0);
}
die unless $seq;

my $options = defined $gene ? "locus=".$gene->{locusTag}
  : "seqDesc=" . encode_entities($seqDesc) . "&" . "seq=$seq";
unless (hasMMSeqsHits($seq)) {
  print redirect(-url => "findHomologs.cgi?${options}");
  exit(0);
}

my ($locusTag, $genome);
if (defined $gene) {
  $locusTag = $gene->{locusTag};
  $genome = gidToGenome($gene->{gid}) || die;
}

my $title = 'Gene neighborhoods';
if ($gene) {
  $title .= " for $gene->{locusTag} and homologs";
} else {
  $title .= " for homologs of " . encode_entities($seqDesc);
}
start_page('title' => $title);

autoflush STDOUT 1; # show preliminary results

if (defined $gene) {
  print p(a({-href => "gene.cgi?locus=$locusTag"}, "$locusTag"),
          "from",
          i($genome->{gtdbSpecies}), $genome->{strain}.":",
          encode_entities($gene->{desc}));
} else {
  print showSequence($seqDesc, $seq),
}
print "\n";

my $hits = getMMSeqsHits($seq);
my $geneHits = hitsToTopGenes( $hits, $n);
if (scalar(@$geneHits) == 0) {
  print p("Sorry, no homologs for this sequence");
  finish_page;
}
if (@$geneHits < $n) {
  print p("Showing all", scalar(@$geneHits), "hits");
} else {
  print p("Showing the top $n hits, out of at least", scalar(@$hits));
}
print "\n";


my $nHits = scalar(@$geneHits);
my $yAt = 5;
my @svgLines = ();
my $xMax = 500;

foreach my $hit (@$geneHits) {
  my $hitGenome = gidToGenome($hit->{gid}) || die;
  my $genomeURL = encode_entities("https://www.ncbi.nlm.nih.gov/assembly/$hitGenome->{gid}/");
  $yAt += 15; # make space for genome label
  push @svgLines,
    qq[<a xlink:href="$genomeURL">],
    qq[<text x="0" y="$yAt" font-style="italic">],
    qq[<title>$hitGenome->{gtdbSpecies} strain $hitGenome->{strain} ($hitGenome->{gid})</title>],
    qq[$hitGenome->{gtdbSpecies}</text></a>];
  $yAt += 6;
  my $nearbyGenes = getNearbyGenes($hit);
  my $mid = ($hit->{begin} + $hit->{end})/2;
  my $showBegin = $mid - $ntShown/2;
  my $showEnd = $mid + $ntShown/2;
  my @showGenes = grep $_->{end} >= $showBegin && $_->{begin} <= $showEnd, @$nearbyGenes;
  foreach my $s (@showGenes) {
    $s->{label} = $s->{locusTag};
    $s->{URL} = "gene.cgi?locus=" . $s->{locusTag};
    $s->{color} = $s->{locusTag} eq $hit->{locusTag} ? "lightblue" : "lightgrey";
  }
  my %genesSvg = genesSvg(\@showGenes,
                          'begin' => $showBegin, 'end' => $showEnd,
                          'yTop' => $yAt,
                          # labels only for the top row
                          'showLabel' => $hit->{locusTag} eq $geneHits->[0]{locusTag},
                          'invert' => $hit->{strand} eq "-");
  push @svgLines, $genesSvg{svg};
  $yAt = max($yAt, $genesSvg{yMax}) + 10;
  $xMax = max($xMax, $genesSvg{xMax});
}
my %scaleBarSvg = scaleBarSvg('xLeft' => $xMax * 0.8,
                              'yTop' => $yAt + 5);
my $svgWidth = max($xMax, $scaleBarSvg{xMax});
my $svgHeight = $scaleBarSvg{yMax};

print join("\n",
           qq[<SVG width="$svgWidth" height="$svgHeight" style="position: relative; left: 1em;>],
           qq[<g transform="scale(1.0)">],
           @svgLines,
           $scaleBarSvg{svg},
           "</g>",
           "</svg>") . "\n";

finish_page();
