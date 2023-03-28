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
# optional arguments:
# order (which subdb to use)
# format = tsv or faa

my $cgi = CGI->new;
my $gid = $cgi->param('gid') || die "Must specify gid";
setOrder(param('order'));

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

my @lineage = ();
foreach my $level (qw[Domain Phylum Class Order Family Genus Species]) {
  my $title = lc($level);
  my $href = "taxon.cgi?level=".lc($level)."&taxon=" . uri_escape($genome->{"gtdb" . $level});
  my $style = "text-decoration: none;";
  if (getOrder() eq ""
      || ($level ne "Domain" && $level ne "Phylum" && $level ne "Class")) {
    $href = addOrderToURL($href);
  } else {
    # link out
    $title = "View this " . lc($level) . " in the main database";
    $style .= "color: black;";
  }
  push @lineage, a({ -title => $title, -style => $style, -href => $href },
                   encode_entities($genome->{"gtdb" . $level}));
}
print p("GTDB taxonomy:", join(", ", @lineage));

my @ncbi = map { s/^[a-z]__//; $_ } split /;/, $genome->{ncbiTaxonomy};
@ncbi = grep $_ ne "", @ncbi;
print p("NCBI taxonomy:",
        join(", ",
             map a({-href => "https://www.ncbi.nlm.nih.gov/taxonomy/?term=" . uri_escape($_),
                    -title => "See taxon at NCBI",
                    -style => "text-decoration: none;"},
                   encode_entities($_)), @ncbi));

my ($nScaffolds) = getDbHandle()->selectrow_array("SELECT COUNT(DISTINCT scaffoldId) FROM Gene WHERE gid = ?",
                                                  {}, $gid);
print p("Scaffolds:", commify($nScaffolds),
        "Genes:",
        a({-href => addOrderToURL("genome.cgi?gid=$gid&format=tsv"),
           -title => "tab-delimited table",
           -style => "text-decoration: none;"},
          commify($genome->{nGenes})),
        "Protein-coding:",
        a({-href => addOrderToURL("genome.cgi?gid=$gid&format=faa"),
           -title => "fasta protein sequences",
           -style => "text-decoration: none;"},
          commify($genome->{nProteins})));

print
  h3("Tools"),
  start_form(-name => 'search', -method => 'GET', -action => 'genomeSearch.cgi'),
  qq{<INPUT type="hidden" name="gid" value="$gid" />},
  orderToHidden(),
  p("Search annotations:", br(),
    qq{<INPUT type="text" name="query" size=40 maxLength=1000 />},
    submit('Go')),
  end_form,
  start_form(-name => 'curated', -method => 'GET', -action => "https://papers.genomics.lbl.gov/cgi-bin/genomeSearch.cgi"),
  qq{<INPUT type="hidden" name="gdb" value="NCBI" />},
  qq{<INPUT type="hidden" name="gid" value="$gid" />},
  p(a({-href => "https://papers.genomics.lbl.gov/cgi-bin/genomeSearch.cgi",
       -title => "Find characterized proteins whose descriptions match the query, and then search the genome for homologs of those proteins",
       -style => "text-decoration: none;"},
      "Curated BLAST:"),
    br(),
    qq{<INPUT type="text" name="query" size=40 maxLength=1000 />},
    checkbox(-name => 'word', -label => 'Match whole words only?', -checked => 0),
    submit('Go'),
    br(),
    span({-style => "font-size: 80%;"}, "Examples: perchlorate or 1.2.1.88")),
  end_form,
  p(a({-href => "https://papers.genomics.lbl.gov/cgi-bin/gapView.cgi",
       -title => "Automated annotation of metabolic pathways",
       -style => "text-decoration: none;"},
      "GapMind")
    ,"for",
    a({-href => "https://papers.genomics.lbl.gov/cgi-bin/gapView.cgi?gdb=NCBI&gid=$gid&set=aa"},
      "amino acid biosynthesis"),
    "or",
    a({-href => "https://papers.genomics.lbl.gov/cgi-bin/gapView.cgi?gdb=NCBI&gid=$gid&set=carbon"},
      "carbon catabolism"),
    span({-style => "font-size: 80%;"}, "(takes 10-40 seconds)")
   );

finish_page();
