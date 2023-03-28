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
# g -- list of locus tags
# optional CGI arguments:
# order (which subdb to use)
my $cgi = CGI->new;
setOrder(param('order'));
my @g = $cgi->param('g');
splice @g, 250 if @g > 250;
my @genes = map { locusTagToGene($_) || die "Unknown locusTag $_" } @g;

# If they all in the same genome, show the genome info at the top,
# and give more details

my %gid = map { $_->{gid} => 1 } @genes;
my $constGenome = gidToGenome( (keys %gid)[0] )
  if scalar(keys %gid) == 1;

my $title = scalar(@g) . " genes";
$title .= " in $constGenome->{gtdbSpecies} $constGenome->{strain}"
  if $constGenome;

start_page('title' => $title);
autoflush STDOUT 1; # show preliminary results

if ($constGenome) {
  my @lineage = ();
  foreach my $level (qw{phylum class order family}) {
    my $tax = $constGenome->{"gtdb" . capitalize($level)};
    my $URL = "taxon.cgi?level=$level&taxon=$tax";
    if (getOrder() ne "" && ($level eq "phylum" || $level eq "class")) {
      push @lineage, a({-href => $URL,
                        -title => "See $level $tax in the main database",
                        -style => "text-decoration: none; color: black;" },
                       $tax);
    } else {
      push @lineage, a({-href => addOrderToURL($URL),
                        -title => $level,
                        -style => "text-decoration: none;"},
                       $tax);
    }
  }
  print p("Lineage:",
          domainHtml($constGenome),
          @lineage,
          a({-href => addOrderToURL("genome.cgi?gid=$constGenome->{gid}"),
             -style => "text-decoration:none;"},
            i($constGenome->{gtdbSpecies})),
          $constGenome->{strain});
  print "<TABLE cellpadding=1 cellspacing=1>\n";
  my @header = qw{Locus Protein Description Begin End Strand};
  print Tr(th([ @header ])), "\n";
  my $iRow = 0;
  foreach my $gene (@genes) {
    my $bgColor = $iRow % 2 == 0 ? "lightgrey" : "white";
    print Tr({-style => "background-color: $bgColor;"},
             td([ a({-href => addOrderToURL("gene.cgi?locus=$gene->{locusTag}") }, $gene->{locusTag}),
                  small($gene->{proteinId}),
                  encode_entities($gene->{desc}),
                  commify($gene->{begin}),
                  commify($gene->{end}) ])
             . td({ -align => "center"}, $gene->{strand})), "\n";
    $iRow++;
  }
  print "</TABLE>\n";
} else {
  foreach my $gene (@genes) {
    my $genome = gidToGenome($gene->{gid}) || die $gene->{gid};
    my @lineage = ();
    foreach my $level (qw{phylum class order family}) {
      my $field = "gtdb" . capitalize($level);
      my $tax = $genome->{$field};
      my $URL = "taxon.cgi?level=$level&taxon=".uri_escape($tax);
      if (getOrder() ne "" && ($level eq "phylum" || $level eq "class")) {
        # Cannot link to this taxon in the subdb.
        push @lineage, a( {-href => $URL,
                           -title => "See $level $tax in the main database",
                           -style => "text-decoration: none; color: black;" },
                          $tax);
      } else {
        push @lineage, a( {-href => addOrderToURL($URL),
                           -title => $level,
                           -style => "text-decoration: none;" },
                          $tax);
      }
    }
    print p(domainHtml($genome),
            ":",
            join(" : ", @lineage),
            ":",
            a({-href => addOrderToURL("genome.cgi?gid=$genome->{gid}"),
               -style => "text-decoration:none;"},
              i($genome->{gtdbSpecies})),
            $genome->{strain},
            br(),
            a({-href => addOrderToURL("gene.cgi?locus=$gene->{locusTag}")}, $gene->{locusTag}),
            small($gene->{proteinId}),
            $gene->{desc});
  }
}
finish_page;
