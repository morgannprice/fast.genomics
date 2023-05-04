#!/usr/bin/perl -w
use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use IO::Handle; # for autoflush
use HTML::Entities;
use URI::Escape;
use lib "../lib";
use neighbor;
# neighborWeb.pm relies on various PaperBLAST libraries
use lib "../../PaperBLAST/lib";
use neighborWeb;
use pbweb qw{commify};

# required CGI arguments:
# taxon and level -- which taxon to show information about
#   (level must be one of domain phylum class order family genus species)
# optional arguments:
# order (which subdb to use)
# format = tsv or faa

my $cgi = CGI->new;
setOrder(param('order'));
my $taxon = $cgi->param('taxon') || die "Must specify taxon";
my $level = $cgi->param('level') || die "Must specify level";
my $taxObj = getDbHandle()->selectrow_hashref("SELECT * FROM Taxon WHERE taxon = ? AND level = ?",
                                              {}, $taxon, $level);
if (!defined $taxObj && getOrder ne "") {
  # Look for this taxon in the main database
  $taxObj = getTopDbHandle()->selectrow_hashref("SELECT * FROM Taxon WHERE taxon = ? AND level = ?",
                                              {}, $taxon, $level);
  if (defined $taxObj) {
    # forward to the main database page
    print redirect(-url => "taxon.cgi?level=$level&taxon=$taxon");
    exit(0);
  }
}

start_page('title' => capitalize($level) . " " . encode_entities($taxon));
autoflush STDOUT 1; # show preliminary results

if (!defined $taxObj) {
  # show links to NCBI and GTDB taxonomy search
  print p("Sorry,", encode_entities($level), encode_entities($taxon),
          "is not in this database. You can try searching for this taxon at",
          a({ -href => "https://gtdb.ecogenomic.org/searches?s=al&q=".uri_escape($taxon),
             -title => "Genome Taxonomy Database"}, "GTDB"),
          "or",
        a({ -href => "https://www.ncbi.nlm.nih.gov/taxonomy/?term=".uri_escape($taxon) },
          "NCBI's taxonomy"));
  print p(start_form(-name => 'input', -method => 'GET', -action => 'findTaxon.cgi'),
          orderToHidden(),
          "Find another taxon:",
          textfield(-name => 'query', -size => 20, -maxlength => 1000, -value => ''),
          submit('Search'),
          end_form);
  finish_page();
}

die "Unknown taxon" unless defined $taxObj;

my $taxa = getTaxa();
my $levelToParent = taxToParts($taxObj, $taxa);
my @taxLevels = taxLevels();
my @lineage = ();
foreach my $parentLevel (@taxLevels) {
  my $parentTaxon = $levelToParent->{$parentLevel};
  push @lineage, a({ -href => addOrderToURL("taxon.cgi?level=$parentLevel&taxon=" . uri_escape($parentTaxon)),
                     -title => $parentLevel },
                   encode_entities($parentTaxon));
  last if $parentLevel eq $level;
}

print p("Lineage:", join(", ", @lineage)), "\n";

my $order = getOrder();
if (getOrder() ne "") {
  my $topDbh = getTopDbHandle();
  my ($class) = $topDbh->selectrow_array(qq{SELECT parent FROM Taxon WHERE level = "order" AND taxon = ?},
                                         {}, $order);
  my ($phylum) = $topDbh->selectrow_array(qq{SELECT parent FROM Taxon WHERE level = "class" AND taxon = ?},
                                         {}, $class);
  my ($domain) = $topDbh->selectrow_array(qq{SELECT parent FROM Taxon WHERE level = "phylum" AND taxon = ?},
                                         {}, $phylum);
  my ($n) = $topDbh->selectrow_array("SELECT COUNT(*) from Genome WHERE gtdbOrder = ?",
                                     {}, $order);
  print p({-style => "font-size: 80%; margin-left: 2em;"},
          "$order is in $domain : $phylum : $class"), "\n";
}

if ($level eq "species" || $taxObj->{nGenomes} == 1) {
  # show the genome(s)
  my $levelCol = "gtdb" . capitalize($level);
  my $genomes = getDbHandle->selectall_arrayref("SELECT * FROM Genome WHERE $levelCol = ?",
                                                { Slice => {} }, $taxon);
  if (@$genomes == 1) {
    my $genome = $genomes->[0];
    print h3("Genome"),
      p(a({-href => addOrderToURL("genome.cgi?gid=" . $genome->{gid})},
          i(encode_entities($genome->{gtdbSpecies})), encode_entities($genome->{strain})));
  } else {
    print h3("Genomes"), start_ul();
    foreach my $genome (@$genomes) {
      print li(a({-href => addOrderToURL("genome.cgi?gid=" . $genome->{gid})},
                 i(encode_entities($genome->{gtdbSpecies})), encode_entities($genome->{strain})));
    }
    print end_ul;
  }
} else {
  # count the genomes
  print p("Genomes:", commify($taxObj->{nGenomes}));
}

my $otherDbStyle = "margin-left: 4em; font-size: 90%;";
if (getOrder() eq "") {
  if ($level eq "order" || $level eq "family" || $level eq "genus" || $level eq "species") {
    # link to subdb if it has more genomes
    my $order = $levelToParent->{order};
    my $nSubGenomes = moreGenomesInSubDb($level, $taxon, $order);
    print p({-style => $otherDbStyle},
            "Or see",
            a({-href => "taxon.cgi?level=$level&taxon=$taxon&order=$order"},
              $nSubGenomes, "genomes"),
            "in the database for", $order, "only")
      if $nSubGenomes > 0;
  }
} else {
  # show a link to the main database, if it has this taxon
  my $topRow = getTopDbHandle()->selectrow_hashref("SELECT * FROM Taxon WHERE level = ? AND taxon = ?",
                                                   {}, $level, $taxon);
  print p({-style => $otherDbStyle},
          "Or see",
          a({-href => "taxon.cgi?level=$level&taxon=$taxon" }, $topRow->{nGenomes}, "genomes"),
          "in the database of diverse bacteria and archaea")
    if $topRow;
}
print "\n";

print p(start_form(-name => 'input', -method => 'GET', -action => 'findTaxon.cgi'),
        orderToHidden(),
        "Find another taxon:",
        textfield(-name => 'query', -size => 20, -maxlength => 1000, -value => ''),
        submit('Search'),
        end_form);

my $childLevel = taxLevelToChild($level);
my $children = [];
$children = getDbHandle->selectall_arrayref("SELECT * FROM Taxon WHERE parent = ? AND level = ?",
                                            { Slice => {} },
                                            $taxon, $childLevel)
  if defined $childLevel;

if (@$children > 0) {
  #my %levelToName = qw{domain Phyla phylum Classes class Orders order Families family Genera genus Species};
  #print h3($levelToName{$level} || die),
  print 
    qq{<TABLE cellpadding=2 cellspacing=2>},
    Tr(th([ "Child $childLevel", "#Genomes" ]));
  my $iRow = 0;
  foreach my $child (@$children) {
    $iRow++;
    my $nG = $child->{nGenomes};
    print Tr({-bgcolor => $iRow % 2 == 0 ? "white" : "lightgrey"},
             td(a({ -style => 'text-decoration: none;',
                    -href => addOrderToURL("taxon.cgi?level=" . $child->{level}
                                           . "&taxon=" . uri_escape($child->{taxon})) },
                  $child->{taxon})),
             td({ -align => "right" }, $nG )), "\n";
  }
  print "</TABLE>";
  print p({-style => "font-size:smaller;"},
          "Only taxa that have high-quality genomes are included in <i>fast.genomics</i>.",
          getOrder() ne "" ? "A limit of 10 genomes per species are included." : "");
}
finish_page();
