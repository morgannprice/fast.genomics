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

# CGI arguments:
# taxon and level -- which taxon to show information about
#   (level must be one of domain phylum class order family genus species)
# optional arguments: format = tsv or faa

my $cgi = CGI->new;
my $taxon = $cgi->param('taxon') || die "Must specify taxon";
my $level = $cgi->param('level') || die "Must specify level";
my $taxObj = getDbHandle()->selectrow_hashref("SELECT * FROM Taxon WHERE taxon = ? AND level = ?",
                                              {}, $taxon, $level);
die "Unknown taxon" unless defined $taxObj;

start_page('title' => capitalize($level) . " " . encode_entities($taxon));
autoflush STDOUT 1; # show preliminary results
my $taxa = getTaxa();
my $levelToParent = taxToParts($taxObj, $taxa);
my @taxLevels = taxLevels();
my @lineage = ();
foreach my $parentLevel (@taxLevels) {
  my $parentTaxon = $levelToParent->{$parentLevel};
  push @lineage, a({ -href => "taxon.cgi?level=$parentLevel&taxon=" . uri_escape($parentTaxon),
                     -title => $parentLevel },
                   encode_entities($parentTaxon));
  last if $parentLevel eq $level;
}
print p("Lineage:", join(", ", @lineage)), "\n";

if ($level eq "species" || $taxObj->{nGenomes} == 1) {
  # show the genome(s)
  my $levelCol = "gtdb" . capitalize($level);
  my $genomes = getDbHandle->selectall_arrayref("SELECT * FROM Genome WHERE $levelCol = ?",
                                                { Slice => {} }, $taxon);
  if (@$genomes == 1) {
    my $genome = $genomes->[0];
    print h3("Genome"),
      p(a({-href => "genome.cgi?gid=" . $genome->{gid}},
        encode_entities($genome->{gtdbSpecies} . " " . $genome->{strain})));
  } else {
    print h3("Genomes"), start_ul();
    foreach my $genome (@$genomes) {
      print li(a({-href => "genome.cgi?gid=" . $genome->{gid}},
                 encode_entities($genome->{gtdbSpecies} . " " . $genome->{strain})));
    }
    print end_ul;
  }
} else {
  # count the genomes
  print p("Genomes:", commify($taxObj->{nGenomes}));
}
print "\n";

print p(start_form(-name => 'input', -method => 'GET', -action => 'findTaxon.cgi'),
        "Or find taxon",
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
                    -href => "taxon.cgi?level=" . $child->{level} . "&taxon=" . uri_escape($child->{taxon}) },
                  $child->{taxon})),
             td({ -align => "right" }, $nG )), "\n";
  }
  print "</TABLE>";
  print p({-style => "font-size:smaller;"},
          "Only taxa that have high-quality genomes are included in <i>fast.genomics</i>");
}
finish_page();
