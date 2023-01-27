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
use neighborWeb;
use pbweb qw{commify};

# CGI arguments:
# locus (a locus tag in the database) or seqDesc and seq
# taxLevel -- defaults to phylum
# Good -- set if show results for good hits only

my $cgi = CGI->new;
my ($gene, $seqDesc, $seq) = getGeneSeqDesc($cgi);
my $taxLevel = $cgi->param('taxLevel') || 'phylum';
die "Unknown taxLevel"
  unless $taxLevel eq "phylum"
  || $taxLevel eq "class"
  || $taxLevel eq "order"
  || $taxLevel eq "family";
my $showGood = $cgi->param('Good') ? 1 : "";

if (defined $gene && ! $seq) {
  # should not be reachable, but redirect to gene page just in case
  print redirect(-url => "gene.cgi?locus=$gene->{locusTag}");
  exit(0);
}
die unless $seq;

unless (hasMMSeqsHits($seq)) {
  print redirect(-url => "findHomologs.cgi?" . geneSeqDescSeqOptions($gene,$seqDesc,$seq));
  exit(0);
}

my $title = 'Taxonomic prevalence';
if ($gene) {
  $title .= " of $gene->{locusTag} and its homologs";
} else {
  $title .= " for homologs of " . encode_entities($seqDesc);
}
start_page('title' => $title);
autoflush STDOUT 1; # show preliminary results

my $hits = getMMSeqsHits($seq);
if (scalar(@$hits) == 0) {
  print p("Sorry, no homologs were found for this sequence");
  finish_page;
}
print "\n";
my $geneHits = hitsToGenes($hits);
my $maxScore = estimateTopScore($geneHits->[0], $seq);
my $scoreThreshold = sprintf("%.1f", 0.3 * $maxScore);
my @goodHits = grep $_->{bits} >= $scoreThreshold, @$geneHits;
my $nGood = scalar(@goodHits);
print p("Loaded", commify(scalar(@$geneHits)), "hits,",
       "including ", commify($nGood), "good hits (&ge; $scoreThreshold bits)");
$geneHits = \@goodHits if $showGood;

print p("Showing the distribution of",
        $showGood ? "good" : "all",
        "hits at the level of", $taxLevel), "\n";
print
  start_form( -name => 'input', -method => 'GET', -action => 'hitTaxa.cgi' ),
  geneSeqDescSeqHidden($gene, $seqDesc, $seq),
  p("Level:",
    popup_menu(-name => 'taxLevel', -values => [ qw(phylum class order family) ], -default => $taxLevel),
    "&nbsp;",
    checkbox(-name => 'Good', -checked => $showGood),
    "hits only?",
    "&nbsp;",
    submit('Change')),
  end_form;

my @links = ();
my $options = geneSeqDescSeqOptions($gene,$seqDesc,$seq);
if (defined $gene) {
  push @links, a({-href => "gene.cgi?$options"}, "gene");
} else {
  push @links, a({-href => "seq.cgi?$options"}, "sequence");
}
push @links, a({-href => "neighbors.cgi?$options"}, "gene neighborhoods")
  . " of homologs";
push @links, a({-href => "downloadHomologs.cgi?$options",
               -title => "tab-delimited table of homologs"}, "download");
print p("Or see", join(" or ", @links));

my @levelsShow = ("Domain", "Phylum");
unless ($taxLevel eq "phylum") {
  push @levelsShow, "Class";
  unless ($taxLevel eq "class") {
    push @levelsShow, "Order";
    unless ($taxLevel eq "order") {
      push @levelsShow, "Family";
    }
  }
}

# taxString is each taxon in this list of levels, joined by ";;;"

my %taxStringN = ();
my $genomes = getDbHandle()->selectall_hashref("SELECT * from Genome", "gid");
foreach my $genome (values %$genomes) {
  $genome->{taxString} = join(";;;", map $genome->{"gtdb".$_}, @levelsShow);
  $taxStringN{ $genome->{taxString} }++;
}

my %taxHits = (); # taxString to list of hits
foreach my $gh (@$geneHits) {
  my $genome = $genomes->{ $gh->{gid} } || die $gh->{gid};
  push @{ $taxHits{$genome->{taxString}} }, $gh;
}
# Each row includes taxString, taxLevels, nHits, nHitGenomes, nGenomes
my @rows = ();
while (my ($taxString, $tHits) = each %taxHits) {
  my $row = { 'taxString' => $taxString, 'nHits' => scalar(@$tHits) };
  $row->{nGenomes} = $taxStringN{$taxString} || die $taxString;
  # Count distinct genomes
  my %gid = (); # genome id to #hits
  foreach my $hit (@$tHits) {
    $gid{ $hit->{gid} }++;
  }
  $row->{nHitGenomes} = scalar(keys %gid);
  push @rows, $row;
}
@rows = sort { $b->{nHits} <=> $a->{nHits}
                 || $a->{taxString} cmp $b->{taxString} } @rows;
my @header = @levelsShow;
$header[0] = "&nbsp;";
push @header, ("#Genomes", "# with hits", "#Hits");
print qq[<TABLE cellpadding=1 cellspacing=1>], "\n";
print Tr(map th($_), @header), "\n";
my $iRow = 0;
my $maxRows = 200;
foreach my $row (@rows) {
  my @out = split /;;;/, $row->{taxString};
  $out[0] = domainHtml($out[0]);
  my $showNHits = commify($row->{nHits});
  my @tHits = @{ $taxHits{ $row->{taxString} } };
  my $truncate = 0;
  my $maxShow = 250;
  if (@tHits > $maxShow) {
    $truncate = 1;
    splice @tHits, $maxShow;
  }
  my $URL = "genes.cgi?" . join("&", map "g=" . $_->{locusTag}, @tHits);
  my $goodString = $showGood ? "good" : "all";
  my $truncateString = $truncate ? " top $maxShow" : "";
  $showNHits = a({ -href => $URL,
                   -style => "text-decoration: none;",
                   -title => "see$truncateString $goodString hits in " . encode_entities($out[-1]) },
                 $showNHits);
  my $bgColor = $iRow % 2 == 0 ? "lightgrey" : "white";
  print Tr({-style => "background-color: $bgColor;"},
           td(\@out),
           td({-style => "text-align: right;"},
              [ commify($row->{nGenomes}), commify($row->{nHitGenomes}), $showNHits ])) . "\n";
  $iRow++;
  last if $iRow >= $maxRows;
}
print "</TABLE>\n";
print p("Table was limited to $maxRows rows")
  if scalar(@rows) > $maxRows;
finish_page();
