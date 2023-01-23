#!/usr/bin/perl -w
use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use IO::Handle; # for autoflush
use DBI;
use HTML::Entities;
use List::Util qw{min max};
use lib "../lib";
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

my $gid = $gene->{gid};
my $genome = $dbh->selectrow_hashref("SELECT * FROM Genome WHERE gid = ?",
                                     {}, $gid)
  || die "Cannot find genome $gid";

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
my $leftShow = max(1, min($gene->{begin}, $center - 5000));
my $rightShow = max($gene->{end}, $center + 5000);
my $ncbiBrowser = "https://www.ncbi.nlm.nih.gov/nuccore/$gene->{scaffoldId}?report=graph"
                 . "&from=$leftShow&to=$rightShow";
push @lines, "Location: " . $gene->{scaffoldId} . " "
  . a({-href => $ncbiBrowser, -title => "NCBI browser"},
      $gene->{begin} . ":" . $gene->{end})
  . " ($gene->{strand})";

print map p({-style => "margin-top: 0.25em; margin-bottom: 0.25em;"}, $_), @lines;

if (defined $seq) {
  my $newline = "%0A";

  my ($psortType, $psortShow) = ("negative", "Gram-negative bacteria");
  ($psortType, $psortShow) = ("positive", "Gram-positive bacteria")
    if $genome->{gtdbPhylum} =~ m/^Firmicutes|Actinobacteriota/;
  ($psortType, $psortShow) = ("archaea", "archaea")
    if $genome->{gtdbDomain} eq "Archaea";

  print h3("Sequence analysis tools");
  my @tools =
    ( a({-href => "http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?query=>${locusTag}$newline$seq"},
        "PaperBLAST") .
      " (search for papers about homologs of this protein)",
      a({-href => "http://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi?seqinput=>${locusTag}$newline$seq"},
        "Search the Conserved Domains Database"),
      a({ -href => "http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/FindSequence.pl?pasted=$seq",
          -title => "Find similar proteins with known structures (PDBsum)"},
      "Search structures"),
      "Predict protein localization: " .
      a({-href => "https://papers.genomics.lbl.gov/cgi-bin/psortb.cgi?name=${locusTag}&type=${psortType}&seq=${seq}",
         -title => "PSORTb v3.0 for $psortShow"},
        "PSORTb") . " ($psortShow)",
      "Find homologs in the " .
      a({-href => "https://iseq.lbl.gov/genomes/seqsearch?sequence=>${locusTag}%0A$seq"},
      "ENIGMA genome browser")
      . " (slow)");
  print map p({-style => "margin-top: 0.25em; margin-bottom: 0.25em;"}, $_), @tools;
  # I could try to look up the locus tag and the protein id in InterPro
  # and link to that...

  print h3("Protein sequence"),
    formatFastaHtml($locusTag . " " . $gene->{desc}, $seq);
}

finish_page();
