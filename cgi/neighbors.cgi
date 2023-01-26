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
use pbweb qw{commify};

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
  print p("Sorry, no homologs were found for this sequence");
  finish_page;
}
if (@$geneHits < $n) {
  print p("Showing all", commify(scalar(@$geneHits)), "hits");
} else {
  print p("Showing the top $n hits, out of at least", commify(scalar(@$hits)));
}
print "\n";


my $nHits = scalar(@$geneHits);
my $yAt = 5;
my @svgLines = ();
my $xMax = 500; # it will actually be much higher for the gene track

my %genes = (); # locusTag => gene to show
foreach my $hit (@$geneHits) {
  my $nearbyGenes = getNearbyGenes($hit);
  my $mid = ($hit->{begin} + $hit->{end})/2;
  my $showBegin = $mid - $ntShown/2;
  my $showEnd = $mid + $ntShown/2;
  my @showGenes = grep $_->{end} >= $showBegin && $_->{begin} <= $showEnd, @$nearbyGenes;
  $hit->{showGenes} = \@showGenes;
  foreach my $s (@showGenes) {
    $genes{ $s->{locusTag} } = $s;
  }
}

# Color the genes based on clustering
my $clusters = clusterGenes(values %genes);
my %locusTagToICluster = ();
foreach my $iCluster (0..(scalar(@$clusters)-1)) {
  my $cluster = $clusters->[$iCluster];
  foreach my $gene (@$cluster) {
    $locusTagToICluster{ $gene->{locusTag} } = $iCluster;
  }
}
# From colorbrewer2.org -- qualitative palette, n=12 (the maximum), and save the 1st color
# for the query and its homologs
my $focalColor = '#a6cee3';
my @clusterColorSet = ('#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
              '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928');
# Add hatch patterns for more effective colors
my @defLines = (); # at top of svg
my @hatchFills = ();
foreach my $iColor (0..(scalar(@clusterColorSet)-1)) {
  my $baseColor = $clusterColorSet[$iColor];
  my $angle = $iColor % 2 == 0 ? 45 :-45;
  push @defLines,
    qq[<pattern id="hatch$iColor" patternUnits="userSpaceOnUse" width="4" height="8"],
    qq[  patternTransform="rotate($angle 2 2)" >],
    qq[<path d="M -1,2 l 6,0" stroke="$baseColor" stroke-width="8" />"],
    qq[</pattern>];
  push @hatchFills, "url(#hatch$iColor)";
}
push @clusterColorSet, @hatchFills;

my @iClusterColors = ();
my %locusTagToColor = ();
foreach my $hit (@$geneHits) {
  $locusTagToColor{ $hit->{locusTag} } = $focalColor;
}
# Color the top tracks first
my $iColor = 0;
foreach my $hit (@$geneHits) {
  foreach my $s (@{ $hit->{showGenes} }) {
    my $sTag = $s->{locusTag};
    next if exists $locusTagToColor{ $sTag };
    if (exists $locusTagToICluster{ $sTag }) {
      my $iCluster = $locusTagToICluster{ $sTag };
      if (!defined $iClusterColors[$iCluster]) {
        $iClusterColors[$iCluster] = $clusterColorSet[ $iColor++ ];
        $iColor = 0 if $iColor >= scalar(@clusterColorSet);
      }
      $locusTagToColor{ $s->{locusTag} } = $iClusterColors[$iCluster];
    }
  }
}

foreach my $hit (@$geneHits) {
  my $hitGenome = gidToGenome($hit->{gid}) || die;
  my $genomeURL = encode_entities("https://www.ncbi.nlm.nih.gov/assembly/$hitGenome->{gid}/");
  $yAt += 25; # make space for genome label
  my $identity = int(100 * $hit->{identity} + 0.5);
  my $qShow = $locusTag || "query";
  my $nHitAA = int(0.5 + ($hit->{end} - $hit->{begin} + 1 - 3)/3);
  my $seqLen = length($seq);
  my $eValue = sprintf("%.2g", $hit->{eValue});
  my $hitDetails = "$hit->{qBegin}:$hit->{qEnd}/$seqLen of $qShow"
    . " is ${identity}% identical to $hit->{sBegin}:$hit->{sEnd}/$nHitAA of $hit->{locusTag} (E = $eValue)";
  my $domainChar = $hitGenome->{gtdbDomain} eq "Bacteria" ? "B" : "A";
  my $domainColor = $domainChar eq "B" ? "blue" : "green";
  push @svgLines,
    qq[<text x="0" y="$yAt" font-size="90%"><title>$hitDetails</title>${identity}% id, $hit->{bits} bits</text>],
    qq[<text x="150" y="$yAt" font-family="bold" font-size="90%" fill=$domainColor>],
    qq[<title>$hitGenome->{gtdbDomain}</title>$domainChar</text>],
    qq[<a xlink:href="$genomeURL">],
    qq[<text x="165" y="$yAt">],
    qq[<title>$hitGenome->{gtdbSpecies} strain $hitGenome->{strain} ($hitGenome->{gid})</title>],
    qq[<tspan font-size="75%"> $hitGenome->{gtdbPhylum} $hitGenome->{gtdbClass} $hitGenome->{gtdbOrder} $hitGenome->{gtdbFamily}</tspan>],
    qq[<tspan font-style="italic" font-size="90%">$hitGenome->{gtdbSpecies}</tspan></text></a>];
  $yAt += 6;
  my $mid = ($hit->{begin} + $hit->{end})/2;
  my $showBegin = $mid - $ntShown/2;
  my $showEnd = $mid + $ntShown/2;
  foreach my $s (@{ $hit->{showGenes} }) {
    $s->{label} = $s->{locusTag};
    $s->{URL} = "gene.cgi?locus=" . $s->{locusTag};
    $s->{color} = $locusTagToColor{$s->{locusTag}} || "lightgrey";
    if ($s->{locusTag} eq $hit->{locusTag}
        && ! (defined $gene && $gene->{locusTag} eq $hit->{locusTag})) {
      my $aaLength = ($hit->{end} - $hit->{begin} + 1 - 3)/3;
      $s->{bar} = { beginFraction => ($hit->{sBegin} - 1)/$aaLength,
                    endFraction => min(1, $hit->{sEnd}/$aaLength),
                    title => $hitDetails,
                    color => $focalColor };
    }
  }
  my %genesSvg = genesSvg($hit->{showGenes},
                          'begin' => $showBegin, 'end' => $showEnd,
                          'yTop' => $yAt,
                          # labels only for the top row
                          'showLabel' => $hit->{locusTag} eq $geneHits->[0]{locusTag},
                          'invert' => $hit->{strand} eq "-");
  push @svgLines, $genesSvg{svg};
  $yAt = max($yAt, $genesSvg{yMax}) + 2;
  $xMax = max($xMax, $genesSvg{xMax});
}
my %scaleBarSvg = scaleBarSvg('xLeft' => $xMax * 0.7,
                              'yTop' => $yAt + 5);
my $svgWidth = max($xMax, $scaleBarSvg{xMax});
my $svgHeight = $scaleBarSvg{yMax};

print join("\n",
           qq[<SVG width="$svgWidth" height="$svgHeight" style="position: relative; left: 1em;">],
           "<defs>",
           @defLines,
           "</defs>",
           qq[<g transform="scale(1.0)">],
           @svgLines,
           $scaleBarSvg{svg},
           "</g>",
           "</svg>") . "\n";

finish_page();
