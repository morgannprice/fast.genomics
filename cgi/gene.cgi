#!/usr/bin/perl -w
use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use IO::Handle; # for autoflush
use DBI;
use HTML::Entities;
use List::Util qw{min max};
use lib "../lib";
use genesSvg;
use neighbor;
# neighborWeb.pm relies on various PaperBLAST libraries
use lib "../../PaperBLAST/lib";
use neighborWeb;

# CGI arguments:
# locus -- a locus tag in the database (with correct capitalization)
my $cgi = CGI->new;
my $locusTag = param('locus') || die "locus must be specified";

start_page('title' => encode_entities($locusTag));
autoflush STDOUT 1; # show preliminary results

my $dbh = getDbHandle();
my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE locusTag = ?",
                                   {}, $locusTag);
if (!defined $gene) {
  print p(b("Sorry, cannot find a gene with locus tag ", encode_entities($locusTag)));
  finish_page();
}

my $focalColor = '#a6cee3';; # see neighbors.pm

my $gid = $gene->{gid};
my $genome = gidToGenome($gid) || die "Cannot find genome $gid";

my @lines = ();
push @lines, "Genome: " . i($genome->{gtdbSpecies}) . " " . $genome->{strain}
  . " " . small("(" . a({-title => "$gid at NCBI",
              -href => "https://www.ncbi.nlm.nih.gov/assembly/$gid/" }, $gid) . ")");
push @lines, "Lineage: " . join(" : ",
                                $genome->{gtdbDomain},
                                $genome->{gtdbPhylum},
                                $genome->{gtdbClass},
                                $genome->{gtdbOrder},
                                $genome->{gtdbFamily});
push @lines, "Description: " . encode_entities($gene->{desc}) if $gene->{desc} ne "";
my $seq;
if ($gene->{proteinId} ne "") {
  my $pId = $gene->{proteinId};
  ($seq) = $dbh->selectrow_array("SELECT sequence FROM Protein WHERE proteinId = ?",
                                 {}, $pId);
  push @lines, "Protein: " .
    a({-href => "https://www.ncbi.nlm.nih.gov/protein/$pId",
       -title => "$pId at NCBI"}, $pId)
    . " (" . length($seq) . " amino acids)";
}
my $center = int(($gene->{begin} + $gene->{end})/2);
my $showBegin =  max(1, min($gene->{begin} - 1000, $center - 4000));
my $showEnd = max($gene->{end} + 1000, $center + 4000);
my $showBegin2 =  max(1, min($gene->{begin} - 1000, $center - 8000));
my $showEnd2 = max($gene->{end} + 1000, $center + 8000);
my $ncbiBrowser = "https://www.ncbi.nlm.nih.gov/nuccore/$gene->{scaffoldId}?report=graph"
                 . "&from=$showBegin&to=$showEnd";
push @lines, "Location: " . $gene->{scaffoldId} . " "
  . a({-href => $ncbiBrowser, -title => "NCBI browser"},
      $gene->{begin} . ":" . $gene->{end})
  . " ($gene->{strand})";

print map p({-style => "margin-top: 0.25em; margin-bottom: 0.25em;"}, $_), @lines;
print "\n";

my $nearbyGenes = getNearbyGenes($gene);
my @showGenes = grep $_->{end} >= $showBegin && $_->{begin} <= $showEnd, @$nearbyGenes;
my @showGenes2 = grep $_->{end} >= $showBegin2 && $_->{begin} <= $showEnd2, @$nearbyGenes;
foreach my $s (@showGenes) {
  $s->{label} = $s->{locusTag};
  $s->{URL} = "gene.cgi?locus=" . $s->{locusTag};
  $s->{color} = $s->{locusTag} eq $gene->{locusTag} ? $focalColor : "lightgrey";
}

my %genesSvg = genesSvg(\@showGenes,
                        'begin' => $showBegin, 'end' => $showEnd,
                        'yTop' => 5,
                        'showLabel' => 1);
my %scaleBarSvg = scaleBarSvg('xLeft' => $genesSvg{xMax} * 0.8,
                              'yTop' => $genesSvg{yMax} + 5);
my $svgWidth = max($genesSvg{xMax}, $scaleBarSvg{xMax});
my $svgHeight = $scaleBarSvg{yMax};

my $genesURL = "genes.cgi?" . join("&", map "g=$_->{locusTag}", @showGenes2);
print join("\n",
           qq{<TABLE width=100% border=0><TR><TD align="left" valign="top"><H3>Gene Neighborhood</h3></TD>},
           qq{<TD align="right" valign="top"><A HREF="$genesURL">or see table</A></TD></TR></TABLE>},
           qq[<SVG width="$svgWidth" height="$svgHeight"
                    style="position: relative; left: 1em;>],
           qq[<g transform="scale(1.0)">],
           $genesSvg{svg},
           $scaleBarSvg{svg},
           "</g>",
           </svg>) . "\n";

if (defined $seq) {
  if (hasMMSeqsHits($seq)) {
    print p("See",
            a({-href => "neighbors.cgi?locus=$locusTag"}, "gene neighborhoods"),
            "or",
            a({-href => "hitTaxa.cgi?locus=$locusTag"}, "taxonomic distribution"),
            "of homologs",
           "or",
           a({-href => "downloadHomologs.cgi?locus=$locusTag",
              -title => "tab-delimited table of homologs"}, "download"));
  } else {
    print p(a{ -href => "findHomologs.cgi?locus=$locusTag" }, "Find homologs with mmseqs2",
            "(fast)");
  }
}
print "\n";

# Genes as table -- maybe another page ?
#my $iMax = scalar(@$nearbyGenes) - 1;
#my ($iThis) = grep $nearbyGenes->[$_]{locusTag} eq $gene->{locusTag}, (0..$iMax);
#die "gene is not near itself" unless defined $iThis;
#my $i1 = max(0, $iThis-5);
#my $i2 = min($iThis+5, $iMax);
#print h3("Nearby Genes");
#for my $i ($i1..$i2) {
#  my $nearbyGene = $nearbyGenes->[$i];
#  my $nlt = $nearbyGene->{locusTag};
#  print p(a({-href => "gene.cgi?locus=$nlt"}, $nlt), $nearbyGene->{strand}, encode_entities($nearbyGene->{desc}));
#}

if (defined $seq) {
  print
    h3("Sequence analysis tools"),
    start_ul,
    join("\n", map li($_), proteinAnalysisLinks($locusTag . " " . $gene->{desc},
                                                $seq, $genome));
    end_ul,
    h3("Protein sequence"),
    formatFastaHtml($locusTag . " " . $gene->{desc}, $seq);
}

finish_page();
