#!/usr/bin/perl -w
use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
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
# gid -- which genome to show information about
# optional arguments: format = tsv or faa

my $cgi = CGI->new;
my $gid = $cgi->param('gid') || die "Must specify gid";

my $genome = gidToGenome($gid);
die "Unknown gid" unless defined $genome;

my $format = $cgi->param('format') || "";

if ($format eq "tsv") {
  print "Content-Type:text/tab-separated-values\n";
  print "Content-Disposition: attachment; filename=genes_$gid.tsv\n\n";
  my $genes = getDbHandle()->selectall_arrayref(qq{ SELECT * FROM Gene WHERE gid = ? },
                                                { Slice => {} }, $gid);
  my @col = qw{scaffoldId begin end strand locusTag proteinId desc};
  print join("\t", @col)."\n";
  foreach my $gene (@$genes) {
    print join("\t", map $gene->{$_}, @col) . "\n";
  }
  exit(0);
} elsif ($format eq "faa") {
  print "Content-Type:text\n";
  print "Content-Disposition: attachment; filename=proteins_$gid.faa\n\n";
  my $prot = getDbHandle()->selectall_arrayref(qq{ SELECT * FROM Gene JOIN Protein USING (proteinId)
                                                   WHERE gid = ? AND proteinId <> "" },
                                               { Slice => {} }, $gid);
  foreach my $p (@$prot) {
    print ">" . $p->{locusTag} . " " . $p->{proteinId} . " " . $p->{desc} . "\n";
    my @seqPieces = $p->{sequence} =~ /.{1,60}/g;
    print join("\n", @seqPieces)."\n";
  }
  exit(0);
}

# else
my $title = encode_entities($genome->{gtdbSpecies} . " " . $genome->{strain});
start_page('title' => $title);

print p("Genome identifier:", a({-href => "https://www.ncbi.nlm.nih.gov/data-hub/genome/$gid" }, $gid));

print p("NCBI assembly name:", encode_entities($genome->{assemblyName}));

print p("GTDB taxonomy:",
        join(", ",
           map a({ -title => lc($_),
                   -style => "text-decoration: none;",
                   -href => "taxon.cgi?level=".lc($_)."&taxon="
                   . uri_escape($genome->{"gtdb" . $_}) },
                 encode_entities($genome->{"gtdb" . $_})),
             qw[Domain Phylum Class Order Family Genus Species]));

my @ncbi = map { s/^[a-z]__//; $_ } split /;/, $genome->{ncbiTaxonomy};
@ncbi = grep $_ ne "", @ncbi;
print p("NCBI taxonomy:",
        join(", ",
             map a({-href => "https://www.ncbi.nlm.nih.gov/taxonomy/?term=" . uri_escape($_),
                    -style => "text-decoration: none;"},
                   encode_entities($_)), @ncbi));

my ($nScaffolds) = getDbHandle()->selectrow_array("SELECT COUNT(DISTINCT scaffoldId) FROM Gene WHERE gid = ?",
                                                  {}, $gid);
print p("Scaffolds:", commify($nScaffolds),
        "Genes:",
        a({-href => "genome.cgi?gid=$gid&format=tsv", -title => "tab-delimited table",
           -style => "text-decoration: none;"},
          commify($genome->{nGenes})),
        "Protein-coding:",
        a({-href => "genome.cgi?gid=$gid&format=faa", -title => "fasta protein sequences",
           -style => "text-decoration: none;"},
          commify($genome->{nProteins})));

finish_page();
