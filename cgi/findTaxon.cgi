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

# optional CGI arguments:
# query -- which taxon or strain to search for
# order -- which subdb to use

my $cgi = CGI->new;
setOrder(param('order'));
my $query = $cgi->param('query');
$query = "" if !defined $query;
$query =~ s/^\s+//;
$query =~ s/\s+$//;

if ($query =~ m/^GC[A-Z]_\d+[.]?\d*$/i) {
  # search for this assembly
  $query = uc($query);
  $query =~ s/[.]\d*$//; # remove version number, if any
  my $query2 = $query;
  if ($query =~ m/^GCA/) {
    $query2 =~ s/^GCA/GCF/;
  } else {
    $query2 =~ s/^GC[A-Z]/GCA/;
  }
  my $genomes = getTopDbHandle()->selectall_arrayref(
    "SELECT * from AllGenome WHERE gid LIKE ? OR gid LIKE ?",
    { Slice => {} }, $query . "%", $query2 . "%");
  if (@$genomes == 1) {
    my $genome = $genomes->[0];
    my $subdbSpec = "";
    $subdbSpec = "&order=" . $genome->{prefix} unless $genome->{inTop};
    print redirect(-url => "genome.cgi?gid=". $genome->{gid} . $subdbSpec);
    exit(0);
  }
}

my $taxList = getDbHandle()->selectall_arrayref("SELECT * from Taxon WHERE taxon LIKE ? LIMIT 200",
                                                {Slice => {} }, $query);

if (@$taxList == 0 && getOrder() ne "" && $query ne "") {
  # Try searching in top-level database?
  $taxList = getTopDbHandle()->selectall_arrayref("SELECT * from Taxon WHERE taxon LIKE ? LIMIT 200",
                                                  {Slice => {} }, $query);
  if (@$taxList > 0) {
    # forward to the top-level database
    print redirect(-url => "findTaxon.cgi?query=".uri_escape($query));
    exit(0);
  }
}

if (@$taxList == 0 && $query ne "") {
  # Search in AllGenome for this species? Put in-top genomes first
  my $genomes = getTopDbHandle()->selectall_arrayref(
    "SELECT * FROM AllGenome WHERE gtdbSpecies LIKE ? ORDER BY inTop DESC",
   { Slice => {} }, $query);
  if (@$genomes > 0) {
    my $genome = $genomes->[0];
    my $subdbSpec = "";
    $subdbSpec = "&order=" . $genome->{prefix} unless $genome->{inTop};
    print redirect(-url => "taxon.cgi?level=species&taxon="
                   . uri_escape($genome->{gtdbSpecies}) . $subdbSpec);
    exit(0);
  }
}

if (@$taxList == 0 && $query ne "" && lc($query) ne "none") {
  # Look for this strain
  # First look in the current database
  my $genomes = getDbHandle()->selectall_arrayref("SELECT * from Genome WHERE strain LIKE ? OR strain LIKE ?",
                                                  { Slice => {} }, $query, $query.";%");
  if (@$genomes > 0) {
    print redirect(-url => addOrderToURL("genome.cgi?gid=" . $genomes->[0]{gid}));
    exit(0);
  }

  # Else look in the top-level AllGenome table
  $genomes = getTopDbHandle()->selectall_arrayref("SELECT * from AllGenome WHERE strain LIKE ? OR strain LIKE ?",
                                                  { Slice => {} }, $query, $query.";%");
  if (@$genomes > 0) {
    my $genome = $genomes->[0];
    my $subdbSpec = "";
    $subdbSpec = "&order=" . $genome->{prefix} unless $genome->{inTop};
    print redirect(-url => "genome.cgi?gid=" . $genomes->[0]{gid} . $subdbSpec);
    exit(0);
  }
}

if (@$taxList == 1) {
  my $taxObj = $taxList->[0];
  print redirect(-url => addOrderToURL("taxon.cgi?level=$taxObj->{level}&taxon=$taxObj->{taxon}"));
  exit(0);
}

start_page('title' => "Taxon search");
if (@$taxList == 0 && $query ne "") {
  my $message = "Sorry, no matching taxa or strains were found for "
    . q{"} . encode_entities($query) . q{".};
  $message .= join(" ",
                   " Or search",
                   a({ -href => "https://gtdb.ecogenomic.org/searches?s=al&q=".uri_escape($query),
                       -title => "Genome Taxonomy Database"}, "GTDB"),
                   "or",
                   a({ -href => "https://www.ncbi.nlm.nih.gov/taxonomy/?term=".uri_escape($query) },
                     "NCBI's taxonomy") . ".");
  print p($message), "\n";
}

if (@$taxList > 0) {
  print p("Found", scalar(@$taxList), "matching taxa in the database of",
          getOrder() || "representative genomes"),
    "<TABLE cellpadding=2 cellspacing=2>",
    Tr(th([ "Level", "Taxon", "#Genomes" ]));
  my $iRow = 0;
  foreach my $taxObj (@$taxList) {
    $iRow++;
    print Tr({-bgcolor => $iRow % 2 == 0 ? "white" : "lightgrey" },
             td($taxObj->{level}),
             td(a({-href => addOrderToURL("taxon.cgi?level=$taxObj->{level}&taxon=$taxObj->{taxon}") },
                  encode_entities($taxObj->{taxon}))),
             td({-align => 'right'}, $taxObj->{nGenomes})), "\n";
  }
  print "</TABLE>", "\n";
}

print p({-style => "font-size:smaller;"},
        "Only taxa that have high-quality genomes are included in <i>fast.genomics</i>.")
  if $query ne "";

print p(start_form(-name => 'input', -method => 'GET', -action => 'findTaxon.cgi'),
        orderToHidden(),
        "Find taxon:",
        qq{<INPUT name='query' type='text' size=20 maxlength=1000 />},
        submit('Search'),
        br(),
        span({-style => "font-size:smaller;"}, "use % for wild cards"),
        end_form);

finish_page();
