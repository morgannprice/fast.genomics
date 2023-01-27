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

# CGI arguments:
# query (optional) -- usually an identifier or locus tag,
#   or a genus and word(s) from a protein description.
#   See neighborWeb::parseGeneQuery for details.
my $cgi = CGI->new;
my $query = $cgi->param('query') || "";

# Redirect to locus tag if the query is an exact match to a
# locus tag or proteinId in our database
my $geneRedirect = locusTagToGene($query); # checks upper case as well
if (! $geneRedirect && $query =~ m/^[a-zA-Z0-9_.]+$/) {
  my $genesFromProteinId = getDbHandle()->selectall_arrayref(
    qq{ SELECT * from Gene WHERE proteinId = ? },
    { Slice => {} }, uc($query));
  $geneRedirect = $genesFromProteinId->[0] if @$genesFromProteinId == 1;
}
if ($geneRedirect) {
  print redirect(-url => "gene.cgi?locus=" . $geneRedirect->{locusTag});
  exit(0);
}

start_page('title' => $query eq "" ? "" : "Gene search");
autoflush STDOUT 1; # show preliminary results

my %query = parseGeneQuery($query);
if (!defined $query{genes} && !defined $query{seq}) {
  print p(b($query{error})) if $query{error};
  print start_form( -name => 'input', -method => 'GET', -action => 'search.cgi' ),
    p("Enter an identifier from UniProt, RefSeq, PDB, or MicrobesOnline,",
      br(),
      "or a protein sequence in FASTA or Uniprot format,",
      br(),
      "or a genus name and a protein description",
      br(),
      textarea( -name  => 'query', -value => '', -cols  => 70, -rows  => 10 ),
      br(),
      br(),
      submit('Search'), reset()),
    end_form;
  finish_page();
}

if (defined $query{genes}) {
  my $genes = $query{genes};
  # show table of genes
  my $gid = $genes->[0]{gid};
  my $genome = gidToGenome($gid) || die $gid;
  die $gid unless $genome;
  print p("Found", scalar(@$genes), " matches in $genome->{gtdbSpecies} $genome->{strain} ($genome->{gid})");
  foreach my $gene (@$genes) {
    my $locusTag = $gene->{locusTag};
    print p(a({href => "gene.cgi?locus=$locusTag"}, $locusTag),
            $gene->{proteinId}, $gene->{desc});
  }
  my $query2 = $query; $query2 =~ s/^\S+\s+//;
  my $curatedURL = "https://papers.genomics.lbl.gov/cgi-bin/genomeSearch.cgi?gdb=NCBI"
    . "&gid=${gid}&query=" . uri_escape($query2);
  print p("Or search for", encode_entities($query2), "in $genome->{gtdbSpecies} $genome->{strain} using", a({-href =>  $curatedURL }, "Curated BLAST"));
  finish_page();
}

finish_page() unless defined $query{seq};

my $seq = $query{seq};
$seq =~ s/[*]//g;
$seq =~ m/^[a-zA-Z]+$/ || die("Sorry, input sequence has invalid characters");
my $seqDesc = $query{seqDesc}
  || length($seq) . " a.a. beginning with " . substr($seq, 0, 10);

print p("Found", encode_entities($seqDesc));
print showSequence($seqDesc, $seq);

my $seqDescE = uri_escape($seqDesc);
if (hasMMSeqsHits($seq)) {
  print p("See",
          a({-href => "neighbors.cgi?seqDesc=${seqDescE}&seq=${seq}"},
            "gene neighborhoods"),
          "or",
          a({-href => "hitTaxa.cgi?seqDesc=${seqDescE}&seq=${seq}"},
            "taxonomic distribution"),
          "of its homologs",
         "or",
          a({-href => "downloadHomologs.cgi?seqDesc=${seqDescE}&seq=${seq}",
             -title => "tab-delimited table of homologs"},
            "download"));
} else {
  print p(a({-href => "findHomologs.cgi?seqDesc=$seqDescE&seq=${seq}"}, "Find homologs with mmseqs2"),
          "(fast)");
}
print
  h3("Other sequence analysis tools"),
  start_ul,
  map li($_), proteinAnalysisLinks($seqDesc, $seq, undef);
print end_ul;

finish_page();

