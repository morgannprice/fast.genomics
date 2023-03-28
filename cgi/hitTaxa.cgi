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

# required CGI arguments:
# locus (a locus tag in the database) or seqDesc and seq
# optional CGI arguments:
# order (which subdb to use)
# taxLevel -- defaults to phylum, for top-level db, or family otherwise
# Good -- set if show results for good hits only
# all -- set to all if show results for all taxa, sorted by taxonomy, instead of sorting by frequency

my $cgi = CGI->new;
setOrder(param('order'));
my ($gene, $seqDesc, $seq) = getGeneSeqDesc($cgi);
my $taxLevel = $cgi->param('taxLevel')
  || (getOrder() eq "" ? "phylum" : "family");
die "Unknown taxLevel"
  unless $taxLevel eq "phylum"
  || $taxLevel eq "class"
  || $taxLevel eq "order"
  || $taxLevel eq "family"
  || $taxLevel eq "genus"
  || $taxLevel eq "species";
my @levelsAllowed = taxLevels();
shift @levelsAllowed; # don't allow the top level (domain or order)
if (getOrder eq "") {
  # no genus or species
  pop @levelsAllowed;
  pop @levelsAllowed;
}
my %levelsAllowed = map { $_ => 1 } @levelsAllowed;
die "Invalid taxLevel $taxLevel\n" unless exists $levelsAllowed{$taxLevel};

my $showGood = $cgi->param('Good') ? 1 : "";
my $all = $cgi->param('all') || "freq";
die "Invalid all" unless $all eq "all" || $all eq "freq";

if (defined $gene && ! $seq) {
  # should not be reachable, but redirect to gene page just in case
  print redirect(-url => addOrderToURL("gene.cgi?locus=$gene->{locusTag}"));
  exit(0);
}
die unless $seq;

unless (hasHits($seq)) {
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

my $hits = getHits($seq);
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
    popup_menu(-name => 'taxLevel', -values => \@levelsAllowed, -default => $taxLevel),
    "&nbsp;",
    checkbox(-name => 'Good', -checked => $showGood),
    "homologs only?",
    "&nbsp;",
    "Show:",
    popup_menu(-name => 'all', -values => [ "all", "freq" ],
               -labels => { "all" => "all taxa", "freq" => "taxa with homologs" },
               -default => $all),
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
               -title => "tab-delimited table of homologs"}, "download homologs");
push @links, a({-href => "compare.cgi?$options",
               -title => "compare presence/absence of homologs and their proximity"}, "compare presence/absence");
print p("Or see", join(", or ", @links));

# Which taxon levels to include in the analysis
my @levels = taxLevels(); # all the ones in the db
# Remove levels after $taxLevel
my ($iLevelAt) = grep $levels[$_] eq $taxLevel, (0..(scalar(@levels)-1));
splice @levels, $iLevelAt + 1;
# Remove order if subdb (the order is always the same)
shift @levels if getOrder() ne "" && $levels[0] eq "order";
my @levelsShow = map capitalize($_), @levels;


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
my @rows;
if ($all eq "freq") {
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
} else {
  # all mode
  my $taxa = getTaxa();
  foreach my $tax (values %{ $taxa->{$taxLevel} }) {
    my $levels = taxToParts($tax, $taxa);
    my $taxString = join(";;;", map $levels->{$_}, @levels);
    my $nHits = 0;
    my $nHitGenomes = 0;
    if (exists $taxHits{$taxString}) {
      $nHits = scalar(@{ $taxHits{$taxString} });
      my %gid = map { $_->{gid} => 1 } @{ $taxHits{$taxString} };
      $nHitGenomes = scalar(keys %gid);
    }
    die $taxString unless defined $taxStringN{$taxString};
    push @rows, { 'taxString' => $taxString,
                  'nGenomes' => $taxStringN{$taxString},
                  'nHitGenomes' => $nHitGenomes,
                  'nHits' => $nHits };
  }
  @rows = sort { $a->{taxString} cmp $b->{taxString} } @rows;
}

my @header = @levelsShow;
$header[0] = "&nbsp;" if $levels[0] eq "domain";
push @header, ("#Genomes", "# with hits", "#Hits",
               a({-title => "Bit score ratio for the best homolog from this group of genomes"},
                 "Max ratio"));
print qq[<TABLE cellpadding=1 cellspacing=1>], "\n";
print Tr(map th($_), @header), "\n";
my $iRow = 0;
my $maxRows;
$maxRows = 200 if $all eq "freq";
foreach my $row (@rows) {
  my @out = split /;;;/, $row->{taxString};
  for (my $i = 0; $i < scalar(@levels); $i++) {
    my $levelThis = $levels[$i];
    if ($levelThis eq "domain") {
      $out[$i] = domainHtml($out[$i]);
    } else {
      $out[$i] = a({ -style => "text-decoration: none;",
                     -href => addOrderToURL("taxon.cgi?level=$levelThis&taxon=".uri_escape($out[$i])) },
                   encode_entities($out[$i]));
    }
  }
  my $showNHits = commify($row->{nHits});
  my $showRatio = "&nbsp;";
  if ($row->{nHits} > 0) {
    my @tHits = @{ $taxHits{ $row->{taxString} } };
    my $maxBits = 0;
    foreach my $hit (@tHits) {
      $maxBits = $hit->{bits} if $hit->{bits} > $maxBits;
    }
    $showRatio = sprintf("%.2f", $maxBits / $maxScore);
    my $truncate = 0;
    my $maxShow = 250;
    if (@tHits > $maxShow) {
      $truncate = 1;
      splice @tHits, $maxShow;
    }
    my $URL = addOrderToURL("genes.cgi?" . join("&", map "g=" . $_->{locusTag}, @tHits));
    my $goodString = $showGood ? "good" : "all";
    my $truncateString = $truncate ? " top $maxShow" : "";
    $showNHits = a({ -href => $URL,
                   -style => "text-decoration: none;",
                   -title => "see$truncateString $goodString hits in " . encode_entities($out[-1]) },
                   $showNHits);
  } else {
    $showNHits = "&nbsp;";
  }
  my $bgColor = $iRow % 2 == 0 ? "lightgrey" : "white";
  print Tr({-style => "background-color: $bgColor;"},
           td(\@out),
           td({-style => "text-align: right;"},
              [ commify($row->{nGenomes}),
                $row->{nHitGenomes} > 0 ? commify($row->{nHitGenomes}) : "&nbsp;",
                $showNHits, $showRatio ])) . "\n";
  $iRow++;
  last if defined $maxRows && $iRow >= $maxRows;
}
print "</TABLE>\n";
print p("Table was limited to $maxRows rows")
  if defined $maxRows && scalar(@rows) >= $maxRows;
finish_page();
