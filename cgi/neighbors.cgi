#!/usr/bin/perl -w
use strict;

use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use IO::Handle; # for autoflush
use DBI;
use HTML::Entities;
use URI::Escape;
use List::Util qw{min max sum};
use lib "../lib";
use neighbor;
use genesSvg;
# neighborWeb.pm relies on various PaperBLAST libraries
use lib "../../PaperBLAST/lib";
use neighborWeb;
use pbweb qw{commify};
use MOTree; # from PaperBLAST/lib/

# Required CGI arguments:
# locus (a locus tag in the database) or seqDesc and seq
#
# Optional arguments:
# order (which subdb to use)
# n -- max number of hits to show
# kb -- how many kilobases to show (default 6 except for long genes)
# hitType -- top, topCollapse, random, or randomAny;
#   top means show top hits;
#   topCollapse means show top hits but collapse cluster x species (subdbs only);
#   random means select that many random good hits;
#   randomAny means random N hits
#   The default is either top or topCollapse
# compact -- use far less vertical space
# format -- empty (default -- shows an svg) or fasta (proteins sequences of homologs shown)
#     or tsv (tab-delimited table of all genes shown)
#     or newick (the tree) or svg
# tree -- set to build a tree (if needed) and render it

my $cgi = CGI->new;
setOrder(param('order'));
my ($gene, $seqDesc, $seq) = getGeneSeqDesc($cgi);
my $n = $cgi->param('n');
$n = 50 unless defined $n && $n =~ m/^\d+$/;
$n = 200 if $n > 200;
my $compact = defined $cgi->param('compact') && $cgi->param('compact') eq "0" ? 0 : 1;
my $hitType = $cgi->param('hitType') || (getOrder ne "" ? 'topCollapse' : 'top');
my $showTree = $cgi->param('tree') || 0;

if (defined $gene && ! $seq) {
  # should not be reachable, but redirect to gene page just in case
  print redirect(-url => addOrderToURL("gene.cgi?locus=$gene->{locusTag}"));
  exit(0);
}
die unless $seq;

my @kbValues = (6, 9, 12, 18, 25, 40);
my $kbShown = $cgi->param('kb');
if (!defined $kbShown) {
  # Default kbShown is 6, unless the gene is long, in
  # which case it is twice its length
  $kbShown = max(6, (length($seq)+1) * 3 * 2/1000);
}
die "Invalid kb" unless $kbShown =~ m/^\d+/;
# Set kbShown to the next highest value in @kbValues
for (my $i = 0; $i < scalar(@kbValues); $i++) {
  if ($kbShown <= $kbValues[$i]) {
    $kbShown = $kbValues[$i];
    last;
  }
}
$kbShown = max(@kbValues) if $kbShown > max(@kbValues);
my $ntShown = $kbShown * 1000;
my $kbWidth = int(0.5 + 150 / ($kbShown/6));

my ($locusTag, $genome);
if (defined $gene) {
  $locusTag = $gene->{locusTag};
  $genome = gidToGenome($gene->{gid}) || die;
}
$kbWidth *= 0.75 if $compact;

unless (hasHits($seq)) {
  print redirect(-url => "findHomologs.cgi?" . geneSeqDescSeqOptions($gene, $seqDesc, $seq));
  exit(0);
}

my $format = $cgi->param('format') || "";
$format = "" unless $format eq "fasta" || $format eq "tsv" || $format eq "newick" || $format eq "svg";
$format = "" unless hasHits($seq);
neighborWeb::setQuietMode() unless $format eq "";

if ($format eq "") {
  my $title = getOrder() eq "" ? "Homologs of "
    : "Homologs from " . getOrder() . " of ";
  if ($gene) {
    $title .= $gene->{locusTag};
  } else {
    $title .= encode_entities($seqDesc);
  }
  start_page('title' => $title);
  autoflush STDOUT 1; # show preliminary results
} else {
  my $fileName = "browse_";
  if (defined $gene) {
    $fileName .= $gene->{locusTag};
  } else {
    my $firstWord = $seqDesc;
    if ($firstWord =~ m/^(sp|tr)[|]/) { # UniProt entries
      $firstWord =~ s/^[a-zA-Z]+[|]//;
    }
    $firstWord =~ s/[| ].*//;
    if ($firstWord =~ m/^[a-zA-Z0-9._-]+$/) {
      $fileName .= $firstWord;
    } else {
      $fileName .= length($seq);
    }
  }
  $fileName .= "." . $format;
  if ($format eq "fasta") {
    print "Content-Type:text\n";
  } elsif ($format eq "tsv") {
    print "Content-Type:text/tab-separated-values\n";
  } elsif ($format eq "newick") {
    $showTree = 1;
    print "Content-Type:text\n";
  } elsif ($format eq "svg") {
    print "Content-Type:image/svg+xml\n";
  }
  print "Content-Disposition: attachment; filename=$fileName\n\n";
}

# Put the xml tag at the end so that any HTML/XML comments, such as from neighborWeb::hitsToGopGenes(),
# are acceptable.
if ($format eq "svg") {
  print qq{<?xml version="1.0"?>}."\n";
}

if (defined $gene) {
  print p(a({-href => addOrderToURL("gene.cgi?locus=$locusTag")}, "$locusTag"),
          "from",
          a({ -href => addOrderToURL("genome.cgi?gid=$genome->{gid}") },
              i(encode_entities($genome->{gtdbSpecies})),
              encode_entities($genome->{strain}))
          .":",
          encode_entities($gene->{desc}))
    if $format eq "";;
}
print "\n" if $format eq "";

my $hits = getHits($seq);
if (scalar(@$hits) == 0) {
  if ($format eq "") {
    print p("Sorry, no homologs were found for this sequence");
    finish_page;
  } else {
    exit(0); # nothing to report
  }
}

my $geneHits; # reference to a list of gene hits
my $listLabel = "hits";
my $nTot = scalar(@$hits);
if ($hitType eq "top") {
  $geneHits = hitsToTopGenes($hits, $n);
} elsif ($hitType eq "topCollapse") {
  die "Cannot collapse for top-level db" if getOrder() eq "";
  $geneHits = hitsToTopGenesClustered($hits, $n,
                                      defined $gene ? $gene->{locusTag} : "");
} elsif ($hitType eq "randomAny") {
  $geneHits = hitsToRandomGenes($hits, $n, 0); # threshold score of 0
} elsif ($hitType eq "random") {
  my $maxScore = estimateTopScore($hits->[0], $seq);
  my $scoreThreshold = sprintf("%.1f", 0.3 * $maxScore);
  $geneHits = hitsToRandomGenes($hits, $n, $scoreThreshold);
  $nTot = scalar(grep { $_->{bits} >= $scoreThreshold } @$hits);
  $listLabel = "good hits (&ge; $scoreThreshold bits)";
} else {
  die "Unknown hitType $hitType";
}

# Force self to be the top hit
if (defined $gene && $geneHits->[0]{locusTag} ne $gene->{locusTag}) {
  my @self = grep $_->{locusTag} eq $gene->{locusTag}, @$geneHits;
  my @other = grep $_->{locusTag} ne $gene->{locusTag}, @$geneHits;
  my @gh = @self;
  push @gh, @other;
  $geneHits = \@gh;
}

my $options = geneSeqDescSeqOptions($gene,$seqDesc,$seq);
if ($format eq "") {
  if ($hitType eq "topCollapse") {
    print p(a({-title => "Within each species, proteins are clustered at 70% identity and 90% overlap, and only one representative is shown from each cluster"},
              "Showing $kbShown kb around",
              @$geneHits < $n ? "all" : "the top",
              commify(scalar(@$geneHits)),
              "clusters x species"));
  } elsif (@$geneHits < $n) {
    print p("Showing $kbShown kb around all",
            commify(scalar(@$geneHits)),
            $listLabel);
  } elsif ($hitType eq "top") {
    print p("Showing $kbShown kb around the top",
            commify(scalar(@$geneHits)),
            "hits, out of at least",
            commify($nTot));
  } else {
    print p("Showing $kbShown kb around",
            commify(scalar(@$geneHits)),
            "random $listLabel, out of at least",
            commify($nTot));
  }
  print "\n";

  # Show options form
  my $randomOption = "";
  my $treeChecked = $showTree ? "CHECKED" : "";
  my $topValues = getOrder() eq "" ? [qw{top randomAny random}]
    : [qw{top topCollapse randomAny random}];
  print
    start_form( -name => 'input', -method => 'GET', -action => 'neighbors.cgi' ),
      geneSeqDescSeqHidden($gene, $seqDesc, $seq),
      p('Hits to show:',
        popup_menu(-name => 'n', -values => [25, 50, 100, 150, 200], -default => $n),
        popup_menu(-name => 'hitType', -values => $topValues,
                   -labels => { 'top' => 'top hits',
                                'topCollapse' => 'top clusters x species',
                                'randomAny' => 'random hits',
                                'random' => 'random good hits' },
                   -default => $hitType),
        "&nbsp;",
        'Kilobases:', popup_menu(-name => 'kb', -values => \@kbValues, -default => $kbShown),
        "&nbsp;",
        qq{<INPUT name='tree' TYPE='checkbox' $treeChecked><label for='tree'>Show tree?</label>},
        "&nbsp;",
        'Taxonomy:',
        popup_menu(-name => 'compact', -values => [1,0], -default => ! $compact,
                   -labels => {'0' => 'detailed', '1' => 'compact'}),
        "&nbsp;",
        submit('Change')),
    end_form;

  my @links = ();
  if (defined $gene) {
    push @links, a({-href => addOrderToURL("gene.cgi?$options")}, "gene");
  } else {
    push @links, a({-href => "seq.cgi?$options"}, "sequence");
  }
  push @links, a({-href => "hitTaxa.cgi?$options"},
                 "taxonomic distribution"). " of its homologs";
  my @downloads = ();
  push @downloads, a({ -href => "neighbors.cgi?${options}&n=${n}&hitType=${hitType}&tree=${showTree}&format=svg&compact=${compact}",
                       -title => "figure in SVG format"}, "SVG");
  push @downloads, a({ -href => "neighbors.cgi?${options}&n=${n}&hitType=${hitType}&format=fasta",
                       -title => "homologs shown (fasta format)"}, "protein sequences");
  push @downloads, a({ -href => "neighbors.cgi?${options}&n=${n}&hitType=${hitType}&kb=${kbShown}&format=tsv",
                       -title => "all genes shown (tab-delimited)"}, "table of genes");
  push @downloads, a({ -href => "neighbors.cgi?${options}&n=${n}&hitType=${hitType}&kb=${kbShown}&tree=1&format=newick",
                       -title => "tree (newick format, rooted)"}, "phylogenetic tree")
    if $showTree;
  print p("Or see", join(" or ", @links).".",
          small("Downloads:", join(", ", @downloads))), "\n";
}

if ($format eq "fasta") {
  foreach my $hit (@$geneHits) {
    my $genome = gidToGenome($hit->{gid}) || die;
    print ">" . join(" ",
                     $hit->{locusTag},
                     $hit->{proteinId},
                     $hit->{bits}, "bits",
                     int(100 * $hit->{identity} + 0.5) . "%", "identity",
                     "from",
                     $genome->{gtdbSpecies}, $genome->{strain},
                     "assemblyId",
                     $hit->{gid})."\n";
    die "No sequence for $hit->{locusTag}" unless $hit->{proteinId} ne "";
    my ($seq) = getDbHandle()->selectrow_array("SELECT sequence FROM Protein WHERE proteinId = ?",
                                               {}, $hit->{proteinId});
    my @seqPieces = $seq =~ /.{1,60}/g;
    print join("\n", @seqPieces)."\n";
  }
  exit(0);
}

my $nHits = scalar(@$geneHits);
# In compact show-tree mode, need extra space at top for the tree's scale bar
my $yAt = $compact && $showTree ? 25 : 1;
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
# for the alignment bar
my $focalColor = "white";
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

if ($format eq "tsv") {
  my @fields = qw{track gtdbDomain gtdbPhylum gtdbClass gtdbOrder gtdbFamily gtdbGenus gtdbSpecies strain assemblyId
                  locusTag proteinId scaffoldId begin end strand desc clusterId};
  print join("\t", @fields)."\n";
  my $iTrack = 0;
  foreach my $hit (@$geneHits) {
    $iTrack++;
    my $hitGenome = gidToGenome($hit->{gid}) || die;
    # additional fields for output
    $hitGenome->{assemblyId} = $hit->{gid};
    $hitGenome->{track} = $iTrack;
    my @show = @{ $hit->{showGenes} };
    @show = reverse @show if $hit->{strand} eq "-";
    foreach  my $s (@show) {
      $s->{clusterId} = "";
      if ($s->{locusTag} eq $hit->{locusTag}) {
        $s->{clusterId} = "query";
      } elsif (exists $locusTagToICluster{ $s->{locusTag} }) {
        $s->{clusterId} = $locusTagToICluster{ $s->{locusTag} };
      }
      my @out = map { exists $hitGenome->{$_} ? $hitGenome->{$_} : $s->{$_} } @fields;
      print join("\t", @out) . "\n";
    }
  }
  exit(0);
}

my $tree = undef;
my $genesLeft = 0;
my @nodesSorted;
my %leafToHit;

if ($showTree) {
  $tree = geneHitsToTree($geneHits);
  $tree->rerootMidpoint();
  if ($format eq "newick") {
    print $tree->toNewick();
    print "\n";
    exit(0);
  }
  $genesLeft = 175;

  # Add node to each gene hit
  my %locusTagToGene = map { $_->{locusTag} => $_ } @$geneHits;
  foreach my $node ($tree->get_leaf_nodes) {
    my $id = $tree->id($node);
    die "Unknown id $id" unless exists $locusTagToGene{$id};
    $locusTagToGene{$id}{node} = $node;
    $leafToHit{$node} = $locusTagToGene{ $tree->id($node) };
  }

  # Order the tree so that the highest-scoring hits are at the top
  # First, compute the max bits for each internal node
  my %nodeToMaxBits = ();
  foreach my $hit (@$geneHits) {
    $nodeToMaxBits{ $hit->{node} } = $hit->{bits};
  }
  # Ensure that the gene itself comes first
  if (defined $gene && $geneHits->[0]{locusTag} eq $gene->{locusTag}) {
    my $selfHit = $geneHits->[0];
    $nodeToMaxBits{ $selfHit->{node} } += 1e9;
  }

  # Reverse of depth first traversal is bottom-up
  my $dfs = $tree->depthfirst();
  foreach my $node (reverse @$dfs) {
    if (!exists $nodeToMaxBits{$node}) {
      my @children = $tree->children($node);
      die unless @children > 0;
      my @values = ();
      foreach my $child (@children) {
        die $child unless exists $nodeToMaxBits{$child};
        push @values, $nodeToMaxBits{$child};
      }
      $nodeToMaxBits{$node} = max(@values);
    }
  }
  # Depth first traversal, sorted by max bits
  $dfs = $tree->depthfirst(\%nodeToMaxBits);
  @nodesSorted = reverse @$dfs;
  my @leavesSorted = grep $tree->is_Leaf($_), @nodesSorted;
  my @hitsSorted = ();
  foreach my $leaf (@leavesSorted) {
    push @hitsSorted, $leafToHit{$leaf};
  }
  # Reassign sorting of hits
  $geneHits = \@hitsSorted;
}

die unless $format eq "" || $format eq "svg";

foreach my $hit (@$geneHits) {
  my $trackHeight = $compact ? 10 : 18;
  my $hitGenome = gidToGenome($hit->{gid}) || die;
  my $genomeURL = encode_entities(addOrderToURL("genome.cgi?gid=$hit->{gid}"));

  my $identity = int(100 * $hit->{identity} + 0.5);
  my $qShow = $locusTag || "query";
  my $nHitAA = int(0.5 + ($hit->{end} - $hit->{begin} + 1 - 3)/3);
  my $seqLen = length($seq);
  my $eValue = sprintf("%.2g", $hit->{eValue});
  my $hitDetails = "$hit->{qBegin}:$hit->{qEnd}/$seqLen of $qShow"
    . " is ${identity}% identical to $hit->{sBegin}:$hit->{sEnd}/$nHitAA of $hit->{locusTag} (E = $eValue)";
  my $hitDetailsShort = "${identity}% identical to the query over " . ($hit->{sEnd}-$hit->{sBegin}+1) . " residues";

  if ($compact) {
    # Minimal genome label at right, with genus or species
    my @levels = ();
    if (getOrder() eq "") {
      @levels = qw{domain phylum class order family species};
    } else {
      @levels = qw{family species};
    }
    my $lineage = join(" : ",
                       map $hitGenome->{"gtdb" . capitalize($_)}, @levels);
    $lineage .= " " . $hitGenome->{strain}
      unless exists $hit->{genes} && scalar(@{ $hit->{genes}}) > 1;
    my $showTax;
    if (getOrder() eq "") {
      $showTax = $hitGenome->{"gtdbGenus"};
    } else {
      $showTax = $hitGenome->{"gtdbSpecies"};
    }
    $showTax = substr($showTax, 0, 45)."." if length($showTax) > 45;
    $showTax = encode_entities($showTax);
    my $xLabel = $genesLeft + $kbWidth * $kbShown + 5;
    my $yCenter = $yAt + $trackHeight/2;
    my $genomeURL = encode_entities(addOrderToURL("genome.cgi?gid=".$hitGenome->{gid}));
    if (exists $hit->{genes} && exists $hit->{nGenomes}
        && scalar(@{ $hit->{genes} }) > 1) {
      # Link to genus not genome if >1 gene
      $genomeURL = encode_entities(addOrderToURL("taxon.cgi?level=species&taxon=".encode_entities($hitGenome->{gtdbSpecies})));
    }
    push @svgLines,
      qq[<text x="$xLabel" y="$yCenter" font-size="70%" dominant-baseline="middle">],
      qq[<a xlink:href="$genomeURL"><title>$lineage</title><tspan font-style="italic">$showTax</tspan></a>];
    if (exists $hit->{genes} && exists $hit->{nGenomes}
        && scalar(@{ $hit->{genes} }) > 1) {
      # information about the collapsed group of genes
      my $nGenes = scalar(@{ $hit->{genes} });
      my $genesURL = encode_entities(addOrderToURL("genes.cgi?" . join("&", map "g=$_", @{ $hit->{genes} })));
      my ($nSpeciesGenomes) = getDbHandle()->selectrow_array(
          qq{SELECT nGenomes FROM Taxon WHERE taxon = ? AND level = "species"}, {}, $hitGenome->{gtdbSpecies});
      my $strainsString = $hit->{nGenomes} > 1 ? "$hit->{nGenomes}/$nSpeciesGenomes strains" : "1/$nSpeciesGenomes strain";
      my $nGenesLabel = "$hit->{nGenomes}/$nSpeciesGenomes "
        . ($nSpeciesGenomes == 1 ? "strain" : "strains");
      $nGenesLabel = "$nGenes in $nGenesLabel" if $nGenes != $hit->{nGenomes};
      push @svgLines, qq[<a xlink:href="$genesURL"><title>$nGenes similar genes (over 70% identity and 90% overlap) in $strainsString of $hitGenome->{gtdbSpecies}</title>($nGenesLabel)</a>];
    }
    push @svgLines, "</text>";
  } else {
    # regular (not compact) view has a line at the top with %identity, bits, and lineage
    $yAt += 25; # make space for genome label
    my @lineage = ();
    my @levels = ();
    if (getOrder() eq "") {
      @levels = qw{domain phylum class order family};
    } else {
      @levels = qw{family}; # skip order, is always the same
    }
    foreach my $level (@levels) {
      my $taxon = $hitGenome->{"gtdb" . capitalize($level)};
      my $taxonURL = encode_entities(addOrderToURL("taxon.cgi?level=$level&taxon=".uri_escape($taxon)));
      my $taxShow;
      my $taxTitle;
      if ($level eq "domain") {
        my $domainChar = $taxon eq "Bacteria" ? "B" : "A";
        my $domainColor = $domainChar eq "B" ? "blue" : "green";
        $taxShow = qq[<tspan font-style="bold" font-size="85%" fill="$domainColor">$domainChar</tspan>];
        $taxTitle = $taxon;
      } else {
        my $size = getOrder() eq "" ? "70%" : "85%";
        $taxShow = qq{<tspan font-size="$size">$taxon</tspan>};
        $taxTitle = $level;
      }
      push @lineage, qq[<a xlink:href="$taxonURL"><title>$taxTitle</title>$taxShow</a>];
    }

    # Add species and sometimes strain information to the lineage (gtdbSpecies includes the genus name)
    # In top-level db, always show the species (but not strain) and link to the genome
    # In sub-db, link species to the taxon, and show the strain, except if this is a collapsed group, show the
    # #genes/#genomes instead.
    my $species = $hitGenome->{gtdbSpecies};
    my $speciesShort = $species;
    $speciesShort = substr($speciesShort, 0, 45) . "." if length($species) > 45;
    my $speciesE = qq{<tspan font-style="italic">}
      . encode_entities($speciesShort) . qq[</tspan>];
    my $strainE = encode_entities($hitGenome->{strain});
    if (getOrder() eq "") {
      push @lineage, qq{<a xlink:href="$genomeURL"><title>strain $strainE</title><tspan font-size="85%">$speciesE</tspan></a>};
    } else {
      my $speciesURL = addOrderToURL("taxon.cgi?level=species&taxon=".uri_escape($species));
      push @lineage, qq[<a xlink:href="$speciesURL"><title>genus/species</title>$speciesE</a>];
      if (exists $hit->{genes} && exists $hit->{nGenomes}
          && scalar(@{ $hit->{genes} }) > 1) {
        # information about the collapsed group of genes
        my $nGenes =  scalar(@{ $hit->{genes} });
        my $genesURL = encode_entities(addOrderToURL("genes.cgi?" . join("&", map "g=$_", @{ $hit->{genes} })));
        my ($nSpeciesGenomes) = getDbHandle()->selectrow_array(
          qq{SELECT nGenomes FROM Taxon WHERE taxon = ? AND level = "species"}, {}, $species);
        die "Unknown species $species" unless $nSpeciesGenomes;
        my $strains = $hit->{nGenomes} > 1 ? "$hit->{nGenomes}/$nSpeciesGenomes strains" : "1/$nSpeciesGenomes strain";
        push @lineage, qq{<a xlink:href="$genesURL"><title>$nGenes similar genes (over 70% identity and 90% overlap) in $strains of $speciesE</title><tspan font-size="90%">($nGenes genes in $strains)</tspan></a>};
      } else {
        # link the strain to the genome
        push @lineage, qq{<a xlink:href="$genomeURL"><title>strain</title>$strainE</a>};
      }
    } # end else has order
    my $xLineage = $genesLeft + 150;
    my $bitsShow = sprintf("%d bits", $hit->{bits});
    # subdb e-values (from LAST) may be of the form 1.23e+03 instead of 1230
    # (also, the bit values are missing the last digit of precision)
    # So I tried changing them to dot 10^3, but it just takes up too much space
    #if ($bitsShow =~ m/e[+]0(\d+)/) {
    #  # &#xb7; is &middot; (but entities are not supported in SVG)
    #  #$bitsShow =~ s!e[+]0(\d+)!&#xb7;10<tspan dy="-5" font-size="80%">$1</tspan> <tspan dy="5">bits</tspan>!;
    #} else {
    #  $bitsShow = "$bitsShow bits";
    #}
    push @svgLines,
      qq[<text x="$genesLeft" y="$yAt" font-size="90%"><title>$hitDetails</title>${identity}% id, $bitsShow</text>],
        qq[<text x="$xLineage" y="$yAt">],
          join(" ",@lineage),
            qq[</text>];
    $yAt += 6;
  } # end else not compact

  $hit->{yTrack} = $yAt + $trackHeight/2; # remember gene location for the future. This is the middle of the track
  my $mid = ($hit->{begin} + $hit->{end})/2;
  my $showBegin = $mid - $ntShown/2;
  my $showEnd = $mid + $ntShown/2;
  foreach my $s (@{ $hit->{showGenes} }) {
    $s->{label} = $s->{locusTag};
    $s->{URL} = addOrderToURL("gene.cgi?locus=" . $s->{locusTag});
    $s->{color} = $locusTagToColor{$s->{locusTag}} || "lightgrey";
    if ($s->{locusTag} eq $hit->{locusTag}) {
      if (! (defined $gene && $gene->{locusTag} eq $hit->{locusTag})) {
        # no bar (or %identity hover text) for the query itself
        my $aaLength = ($hit->{end} - $hit->{begin} + 1 - 3)/3;
        $s->{bar} = { beginFraction => ($hit->{sBegin} - 1)/$aaLength,
                      endFraction => min(1, $hit->{sEnd}/$aaLength),
                      URL => "alignPair.cgi?${options}&locus2=" . $s->{locusTag},
                      # in compact mode, difficult to click on one vs. other so have the same
                      # information in the hover text for the bar as for the ORF
                      title => $compact ? $s->{desc}. " ($hitDetailsShort, see alignment)": $hitDetails,
                      # The original color '#a6cee3' was too close to white, use a cyan instead
                      color => '#44DDDD'  };
        $s->{desc} .= " ($hitDetailsShort)";
      }
    }
  }
  my ($scaffoldLen) = getDbHandle->selectrow_array(
    "SELECT length FROM Scaffold WHERE gid = ? AND scaffoldId = ?",
    {}, $hit->{gid}, $hit->{scaffoldId});
  my %genesSvg = genesSvg($hit->{showGenes},
                          'begin' => $showBegin, 'end' => $showEnd,
                          'kbWidth' => $kbWidth,
                          'xLeft' => $genesLeft,
                          'geneHeight' => $trackHeight,
                          'yTop' => $yAt,
                          # labels only for the top row, and only if not compact
                          'showLabel' => $hit->{locusTag} eq $geneHits->[0]{locusTag} && ! $compact,
                          'invert' => $hit->{strand} eq "-",
                          'scaffoldLength' => $scaffoldLen);
  push @svgLines, $genesSvg{svg};
  $yAt = max($yAt, $genesSvg{yMax}) + ($compact ? 0 : 2);
  $xMax = max($xMax, $genesSvg{xMax});
}
my %scaleBarSvg = scaleBarSvg('xLeft' => $genesLeft + ($xMax - $genesLeft) * 0.7,
                              'yTop' => $yAt + 5);
push @svgLines, $scaleBarSvg{svg};
# Would be ~900; increase to make sure species name shows up in top level view on some browsers
my $svgWidth = max($xMax, $scaleBarSvg{xMax}, 1000);
if ($compact) { # do not truncate the taxon labels at the right
  $svgWidth = max($xMax, $genesLeft + $kbWidth * $kbShown + 400);
}
my $svgHeight = $scaleBarSvg{yMax};

# Now that the gene tracks are laid out, we have a position for each gene, so
# lay out the tree
if ($tree) {
  my %nodeY;
  # Leaves are at yTrack (top of the gene diagram)
  foreach my $hit (@$geneHits) {
    my $node = $hit->{node};
    $nodeY{ $hit->{node} } = $hit->{yTrack};
  }
  # For internal nodes, y is average of childrens' y
  my $dfs = $tree->depthfirst;
  foreach my $node (reverse @$dfs) {
    if (!exists $nodeY{$node}) {
      my @children = $tree->children($node);
      die unless @children > 0;
      my @values = ();
      foreach my $child (@children) {
        die unless defined $nodeY{$child};
        push @values, $nodeY{$child};
      }
      $nodeY{$node} = sum(@values) / scalar(@values);
    }
  }

  my %nodeX; # proportionate to distance from root
  my $treeLeft = 1;
  my $treeRight = $genesLeft - 5;
  my $nodeDepth = $tree->nodeDepth();   # Distance from root to every node
  my $maxDepth = max(values %$nodeDepth);
  $maxDepth = 0.001 if $maxDepth < 0.001;
  while (my ($node, $depth) = each %$nodeDepth) {
    $nodeX{$node} = $treeLeft + ($treeRight - $treeLeft) * $depth/$maxDepth;
  }

  my %nodeTitle = ();
  my %nodeURL = ();
  my %nodeColor = (); # defaults to black
  my %nodeRadius = ();
  foreach my $node (@$dfs) {
    if ($tree->is_Leaf($node)) {
      my $hit = $leafToHit{$node} || die;
      my $hitGenome = gidToGenome($hit->{gid}) || die;
      $nodeTitle{$node} = encode_entities(join(" ", $hit->{locusTag}, "from",
                                               $hitGenome->{gtdbSpecies}, $hitGenome->{strain}));
      $nodeURL{$node} = addOrderToURL("gene.cgi?locus=" . $hit->{locusTag});
      $nodeRadius{$node} = $compact ? 3 : 4;
    } else {
      my $id = $tree->id($node);
      if (defined $id && $id ne "") {
        $nodeRadius{$node} = $compact ? 3 : 4;
        $nodeTitle{$node} = sprintf("Support: %.2f", $tree->id($node));
        $nodeColor{$node} = "darkgrey" if $id < 0.8;
      }
    }
  }
  push @svgLines, $tree->drawSvg('nodeX' => \%nodeX,
                                 'nodeY' => \%nodeY,
                                 'nodeURL' => \%nodeURL,
                                 'nodeTitle' => \%nodeTitle,
                                 'nodeColor' => \%nodeColor,
                                 'nodeRadius' => \%nodeRadius);
  push @svgLines, $tree->drawSvgScaleBar('nodeX' => \%nodeX, 'y' => 10);
}

if ($format eq "svg") {
  # I tried using version="1.1" and setting xml:base to allow the links to work after the SVG is downloaded,
  # but that didn't work
  # Instead, rewrite relative links to point at http://fast.genomics.lbl.gov/cgi/
  foreach (@svgLines) { s! xlink:href *= *"([^/][^"]+)"! xlink:href="https://fast.genomics.lbl.gov/cgi/$1"!g; }
  # xmlns:xlink= allows the SVG to contain xlink:href attributes for <A> tags
  print join("\n",
             qq{<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="$svgWidth" height="$svgHeight">},
             "<defs>", @defLines, "</defs>",
             qq{<g transform="scale(1)">},
             @svgLines,
             "</g>",
             "</svg>")."\n";
  exit(0);
}

print join("\n",
           qq[<SVG width="$svgWidth" height="$svgHeight" style="position: relative; left: 1em;">],
           "<defs>", @defLines, "</defs>",
           qq[<g transform="scale(1.0)">],
           @svgLines,
           "</g>",
           "</svg>") . "\n";

if ($tree) {
  print p({-style => "font-size: 90%;"},
          "The alignment was computed with",
          a({-href => "https://drive5.com/muscle/downloads_v3.htm"}, "MUSCLE 3"),
          "and fast settings (-maxiters 2 -maxmb 1000). Only the",
          "part of each sequence that is similar to the query was used.",
          "The phylogenetic tree was computed with",
          a({-href => "http://www.microbesonline.org/fasttree/"}, "FastTree 2"),
          "using default settings",
          "and was then rooted at the midpoint.");
}

finish_page();
