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

# required CGI arguments:
# gid -- which genome to search in
# query -- terms to search, either match the locus tag or protein id, or search in the description
# optional CGI arguments:
# order (which subdb to use)
my $cgi = CGI->new;
setOrder(param('order'));
my $gid = $cgi->param('gid') || "Must specify gid";
my $query = $cgi->param('query') || "";
$query =~ s/^\s+//;
$query =~ s/\s+$//;

my $geneRedirect = locusTagToGene($query); # checks upper case as well
if (! $geneRedirect && $query =~ m/^[a-zA-Z0-9_.]+$/) {
  my $genesFromProteinId = getDbHandle()->selectall_arrayref(
    qq{ SELECT * from Gene WHERE proteinId = ? },
    { Slice => {} }, uc($query));
  $geneRedirect = $genesFromProteinId->[0] if @$genesFromProteinId == 1;
}
if ($geneRedirect) {
  print redirect(-url => addOrderTourl("gene.cgi?locus=" . $geneRedirect->{locusTag}));
  exit(0);
}

my $genome = getDbHandle()->selectrow_hashref("SELECT * from Genome WHERE gid = ?", {}, $gid);
die "Unknown gid $gid\n" unless defined $genome;
start_page('title' => "Search in " . encode_entities($genome->{gtdbSpecies} . " " . $genome->{strain}));
autoflush STDOUT 1; # show preliminary results

if ($query ne "") {
  my $genes = getDbHandle()->selectall_arrayref(qq{ SELECT * FROM Gene
                                           WHERE gid=?
                                           AND (desc LIKE ? OR desc LIKE ? OR desc LIKE ? OR desc LIKE ?)
                                           LIMIT 200 },
                                       { Slice => {} },
                                       $gid,
                                       "${query}%", "%-${query}%",
                                       "% ${query}%", "% (${query}%");
  my $genomeLink = a({-href => addOrderToURL("genome.cgi?gid=$gid")},
                     i(encode_entities($genome->{gtdbSpecies})), 
                     encode_entities($genome->{strain}));
  if (@$genes == 0) {
    print p("Sorry, no matches for", encode_entities($query), "in", $genomeLink);
  } else {
    print p("Found", scalar(@$genes), "genes matching", encode_entities($query), "in", $genomeLink);
    foreach my $gene (@$genes) {
      my $locusTag = $gene->{locusTag};
      print p(a({href => addOrderToURL("gene.cgi?locus=$locusTag")}, $locusTag),
              $gene->{proteinId}, $gene->{desc});
    }
  }
}

# Search form
print
  start_form(-name => 'search', -method => 'GET', -action => 'genomeSearch.cgi'),
  orderToHidden(),
  qq{<INPUT type="hidden" name="gid" value="$gid" />},
  p("Try another search:", br(),
    qq{<INPUT type="text" name="query" size=40 maxLength=1000 />},
    submit('Go')),
  end_form,
  start_form(-name => 'curated', -method => 'GET', -action => "https://papers.genomics.lbl.gov/cgi-bin/genomeSearch.cgi"),
  qq{<INPUT type="hidden" name="gdb" value="NCBI" />},
  qq{<INPUT type="hidden" name="gid" value="$gid" />},
  p("Or try",
    a({-title => "Examples: perchlorate or 1.2.1.88"}, "Curated BLAST:"), br(),
    qq{<INPUT type="text" name="query" size=40 maxLength=1000 />},
    checkbox(-name => 'word', -label => 'Match whole words only?', -checked => 0),
    submit('Go')),
  end_form;

finish_page();
