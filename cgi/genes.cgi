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

# CGI arguments:
# g -- list of locus tags
my $cgi = CGI->new;
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
  print p("Lineage:",
          domainHtml($constGenome),
          map(a({-href => "taxon.cgi?level=".lc($_)."&taxon=" . uri_escape($constGenome->{"gtdb$_"}),
                 -title => lc($_),
                 -style => "text-decoration: none;"}, $constGenome->{"gtdb$_"}),
                qw{Phylum Class Order Family}),
          a({-href => "genome.cgi?gid=$constGenome->{gid}",
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
             td([ a({-href => "gene.cgi?locus=$gene->{locusTag}" }, $gene->{locusTag}),
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
    print p(domainHtml($genome),
            map(a({ -href => "taxon.cgi?level=".lc($_)."&taxon=" . uri_escape($genome->{"gtdb$_"}),
                    -title => lc($_),
                    -style => "text-decoration: none;"}, $genome->{"gtdb$_"}),
                qw{Phylum Class Order Family}),
            a({-href => "genome.cgi?gid=$genome->{gid}",
               -style => "text-decoration:none;"},
              i($genome->{gtdbSpecies})),
            $genome->{strain},
            br(),
            a({-href => "gene.cgi?locus=$gene->{locusTag}"}, $gene->{locusTag}),
            small($gene->{proteinId}),
            $gene->{desc});
  }
}
finish_page;
