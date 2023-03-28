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
# query -- which taxon to search for
# order -- which subdb to use

my $cgi = CGI->new;
setOrder(param('order'));
my $query = $cgi->param('query');
$query = "" if !defined $query;
$query =~ s/^\s+//;
$query =~ s/\s+$//;

my $taxList = getDbHandle()->selectall_arrayref("SELECT * from Taxon WHERE taxon LIKE ?",
                                                {Slice => {} }, $query);

if (@$taxList == 1) {
  my $taxObj = $taxList->[0];
  print redirect(-url => addOrderToURL("taxon.cgi?level=$taxObj->{level}&taxon=$taxObj->{taxon}"));
  exit(0);
}

start_page('title' => "Taxon search");
if (@$taxList == 0 && $query ne "") {
  my $message = "Sorry, no matching taxa were found for "
    . q{"} . encode_entities($query) . q{".};
  $message .= " Search the "
    . a({-href => "findTaxon.cgi?query=" . uri_escape($query)}, "main database")
    . " instead."
      if getOrder() ne "";
  print p($message), "\n";
}

if (@$taxList > 0) {
  print p(scalar(@$taxList), "matching taxa"),
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
        "Only taxa that have high-quality genomes are included in <i>fast.genomics</i>")
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
