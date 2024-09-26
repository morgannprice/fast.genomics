#!/usr/bin/perl -w
use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use IO::Handle; # for autoflush
use DBI;
use HTML::Entities;
use LWP::Simple;
use List::Util qw{min max};
use lib "../lib";
use genesSvg;
use neighbor;
# neighborWeb.pm relies on various PaperBLAST libraries
use lib "../../PaperBLAST/lib";
use neighborWeb;
use pbweb qw{commify};

# required CGI arguments:
# locus -- a locus tag in the database (with correct capitalization)
# optional CGI arguments:
# order (which subdb to use)
# get=ntseq (report the nucleotide sequence as text)
my $cgi = CGI->new;
my $locusTag = param('locus') || die "locus must be specified";
setOrder(param('order'));

my $get = param('get') || "";
die unless $get eq "" || $get eq "ntseq";
if ($get eq "") {
  start_page('title' => encode_entities($locusTag));
  autoflush STDOUT 1; # show preliminary results
} else {
  print "Content-Type:text\n\n";
}
  

my $dbh = getDbHandle();
my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE locusTag = ?",
                                   {}, $locusTag);
if (!defined $gene) {
  die "No such locustag" if $get ne "";
  print p(b("Sorry, cannot find a gene with locus tag ", encode_entities($locusTag)));
  finish_page();
}

if ($get eq "ntseq") {
  my $fetch = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.cgi?"
    . join("&", "db=nuccore",
           "rettype=fasta",
           "id=" . $gene->{scaffoldId},
           "seq_start=" . $gene->{begin},
           "seq_stop=" . $gene->{end},
           "strand=" . ($gene->{strand} eq "+" ? 1 : 2));
  my $fasta  = get($fetch);
  die "Cannot fetch fasta from NCBI using\n$fetch\n" unless $fasta;
  my $header = $locusTag . " " . $gene->{scaffoldId} . ":" . $gene->{begin} . "-" . $gene->{end} . ":" . $gene->{strand};
  $fasta =~ s/^>[^\n]+/>$header/;
  print $fasta;
  exit(0);
}

my $focalColor = '#a6cee3'; # see neighbors.pm

my $gid = $gene->{gid};
my $genome = gidToGenome($gid) || die "Cannot find genome $gid";

my @lines = ();
push @lines, "Genome: "
  . a({-href => addOrderToURL("genome.cgi?gid=$gid")},
      i($genome->{gtdbSpecies}) . " " . $genome->{strain})
  . " " . small("(" . a({-title => "$gid at NCBI",
              -href => "https://www.ncbi.nlm.nih.gov/assembly/$gid/" }, $gid) . ")");
push @lines, "Lineage: "
  . join(" : ",
         map a({ -href => addOrderToURL("taxon.cgi?level=" . lc($_) . "&taxon=" . $genome->{"gtdb".$_}),
                 -title => lc($_) },
               $genome->{"gtdb".$_}),
         qw{Domain Phylum Class Order Family});
push @lines, "Description: " . encode_entities($gene->{desc}) if $gene->{desc} ne "";
my $seq;
if ($gene->{proteinId} ne "") {
  my $pId = $gene->{proteinId};
  ($seq) = $dbh->selectrow_array("SELECT sequence FROM Protein WHERE proteinId = ?",
                                 {}, $pId);
  push @lines, "Protein: " .
    a({-href => "https://www.ncbi.nlm.nih.gov/protein/$pId",
       -title => "$pId at NCBI"}, $pId)
    . " (" . commify(length($seq)) . " amino acids)";
}
my $center = int(($gene->{begin} + $gene->{end})/2);
# Nearby genes shown in the graphic or in the NCBI link
my $showBegin =  max(1, min($gene->{begin} - 2000, $center - 4000));
my $showEnd = max($gene->{end} + 2000, $center + 4000);
# Nearby genes shown in the linked-to table
my $showBegin2 =  max(1, min($gene->{begin} - 2000, $center - 8000));
my $showEnd2 = max($gene->{end} + 2000, $center + 8000);
my $ncbiBrowser = "https://www.ncbi.nlm.nih.gov/nuccore/$gene->{scaffoldId}?report=graph"
                 . "&from=$showBegin&to=$showEnd";
my $scaffold = getDbHandle()->selectrow_hashref("SELECT * from Scaffold WHERE scaffoldId = ?",
                                                {}, $gene->{scaffoldId});

my $ntURL = addOrderToURL("gene.cgi?locus=$locusTag&get=ntseq");

push @lines, "Location: "
  . a({-title => $scaffold->{scaffoldDesc}, -style => "text-decoration: none;"},
      $gene->{scaffoldId})
  . " "
  . a({-href => $ncbiBrowser, -title => "NCBI browser"},
      $gene->{begin} . ":" . $gene->{end})
  . " ($gene->{strand})"
  . " " . small(a({-href => $ntURL, -title => "fetch the nucleotide sequence from NCBI"}, "nt. sequence"));

print map p({-style => "margin-top: 0.25em; margin-bottom: 0.25em;"}, $_), @lines;
print "\n";

my $nearbyGenes = getNearbyGenes($gene);
my @showGenes = grep $_->{end} >= $showBegin && $_->{begin} <= $showEnd, @$nearbyGenes;
my @showGenes2 = grep $_->{end} >= $showBegin2 && $_->{begin} <= $showEnd2, @$nearbyGenes;
foreach my $s (@showGenes) {
  $s->{label} = $s->{locusTag};
  $s->{URL} = addOrderToURL("gene.cgi?locus=" . $s->{locusTag});
  $s->{color} = $s->{locusTag} eq $gene->{locusTag} ? $focalColor : "lightgrey";
}

my $showWidth = $showEnd - $showBegin + 1;
my $kbWidth = $showWidth < 9*1000  ?  150 : 150/($showWidth/9000);
my %genesSvg = genesSvg(\@showGenes,
                        'kbWidth' => $kbWidth,
                        'begin' => $showBegin, 'end' => $showEnd,
                        'yTop' => 5,
                        'scaffoldLength' => $scaffold->{length},
                        'showLabel' => 1);
my %scaleBarSvg = scaleBarSvg('xLeft' => $genesSvg{xMax} * 0.8,
                              'yTop' => $genesSvg{yMax} + 5);
my $svgWidth = max($genesSvg{xMax}, $scaleBarSvg{xMax});
my $svgHeight = $scaleBarSvg{yMax};

my $genesURL = addOrderToURL("genes.cgi?" . join("&", map "g=$_->{locusTag}", @showGenes2));
print join("\n",
           qq{<TABLE width=100% border=0><TR><TD align="left" valign="top"><H3>Gene Neighborhood</h3></TD>},
           qq{<TD align="right" valign="top"><A HREF="$genesURL">or see table</A></TD></TR></TABLE>},
           qq[<SVG width="$svgWidth" height="$svgHeight"
                    style="position: relative; left: 1em;>],
           qq[<g transform="scale(1.0)">],
           $genesSvg{svg},
           $scaleBarSvg{svg},
           "</g>",
           "</svg>") . "\n";

if (defined $seq) {
  if (hasHits($seq)) {
    my $hitsFile = hitsFile($seq);
    print "<!-- hits are in $hitsFile -->\n"; # aids debugging
    print p("See",
            join(", or ",
                 a({-href => addOrderToURL("neighbors.cgi?locus=$locusTag")}, "gene neighborhoods"),
                 a({-href => addOrderToURL("hitTaxa.cgi?locus=$locusTag")}, "taxonomic distribution")
                 . " of homologs",
                 a({-href => addOrderToURL("downloadHomologs.cgi?locus=$locusTag"),
                    -title => "tab-delimited table of homologs"}, "download homologs"),
                 a({-href => addOrderToURL("compare.cgi?locus=$locusTag"),
                    -title => "compare presence/absence of homologs and their proximity"},
                   "compare presence/absence")));
  } else {
    print p(a({ -href => addOrderToURL("findHomologs.cgi?locus=$locusTag") },
              getOrder() eq "" ? "Find homologs with mmseqs2" : "Find homologs with clustered BLAST"),
            "(fast)");
  }

  if (getOrder() eq "") {
    # Link to subdb homologs, if subdb has more genomes for this order
    my $order = $genome->{gtdbOrder};
    my ($nSubGenomes) = moreGenomesInSubDb("order", $order, $order);
    if ($nSubGenomes > 0) {
      print p("Or find",
              a({-href => "findHomologs.cgi?locus=$locusTag&order=$order"},
                "homologs in", commify($nSubGenomes), $order));
    }
  } else {
    # Link to main db homologs
    # Is this gene in the main db?
    my $main = getTopDbHandle()->selectrow_arrayref(
       "SELECT * FROM Gene WHERE locusTag = ? AND proteinId = ?",
       {}, $locusTag, $gene->{proteinId});
    my $spec;
    if ($main) {
      $spec = "locus=$locusTag";
    } else {
      $spec = "seq=$seq&seqDesc=" . encode_entities("$locusTag $gene->{proteinId} $gene->{desc}");
    }
    print p("Or find",
              a({-href => "findHomologs.cgi?${spec}"},
                "homologs in diverse bacteria and archaea"));
  }

  print showSequence($locusTag . " " . $gene->{desc}, $seq);
}

print "\n";

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
