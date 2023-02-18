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
      textarea( -name  => 'query', -value => '', -cols  => 80, -rows  => 4 ),
      br(),
      br(),
      submit('Search')),
    end_form;
  my @examples = ("ING2E5A_RS06865",
                  "3osdA",
                  "Escherichia thymidylate synthase");
  my @exampleLabels = ("by locus tag", "by Protein Data Bank entry", "by genus and description");
  my @exampleLinks = ();
  for (my $i = 0; $i < scalar(@examples); $i++) {
    push @exampleLinks, a({ -href => "search.cgi?query=" . uri_escape($examples[$i]),
                            -style => "text-decoration: none;",
                            -title => $exampleLabels[$i] },
                          $examples[$i]);
  }
  print p({-style => "font-size: 90%;"},
          "Example searches:", join(", ", @exampleLinks));
  my ($nGenomes) = getDbHandle()->selectrow_array("SELECT COUNT(*) FROM Genome");
  $nGenomes = commify($nGenomes);
  print <<END
<H3>About <i>fast.genomics</i></H3>

<P>Fast.genomics includes representative genomes for $nGenomes genera
of Bacteria and Archaea. These were classified by using
the <A HREF="https://gtdb.ecogenomic.org/">Genome Tree Database</A>.
Only high-quality genomes are included. Potential chimeras were excluded using
<A HREF="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02393-0">GUNC</A>.
Where possible, genomes were
taken from NCBI's <A
HREF="https://www.ncbi.nlm.nih.gov/refseq/">RefSeq</A>.

<P>Fast.genomics uses <A
HREF="https://github.com/soedinglab/MMseqs2">mmseqs2</A> to find
homologs for a protein sequence of interest. This usually takes a few
seconds. To speed up the search, fast.genomics splits the protein database into pieces,
which allows parallel
analysis of a single query.
Fast.genomics also keeps the indexes in memory.
The protein of interest need not be in fast.genomics' database.

<P>Once the homologs are identified, fast.genomics can quickly show:

<UL>
<LI><A HREF="neighbors.cgi?locus=ING2E5A_RS06865">Gene neighborhoods</A>
<LI><A HREF="hitTaxa.cgi?locus=ING2E5A_RS06865">Which taxa contain homologs</A>
<LI>Compare the <A HREF="compare.cgi?locus=ING2E5A_RS06865&query2=ING2E5A_RS06860">presence/absence</A> of two proteins
</UL>

<P style="font-size:90%;">(These examples are for a putative 3-ketoglycoside hydrolase, ING2E5A_RS06865. This family of proteins was formerly known as DUF1080.)</P>
END
    ;
  my ($nPhyla) = getDbHandle()->selectrow_array("SELECT COUNT(DISTINCT gtdbPhylum) FROM Genome");
  $nPhyla = commify($nPhyla);
  my ($nClasses) = getDbHandle()->selectrow_array("SELECT COUNT(DISTINCT gtdbClass) FROM Genome");
  $nClasses = commify($nClasses);
  my ($nOrders) = getDbHandle()->selectrow_array("SELECT COUNT(DISTINCT gtdbOrder) FROM Genome");
  $nOrders = commify($nOrders);
  my ($nFamilies) = getDbHandle()->selectrow_array("SELECT COUNT(DISTINCT gtdbFamily) FROM Genome");
  $nFamilies = commify($nFamilies);
  my $dbDate = `date -r ../data/neighbor.db '+%b %-d, %Y'`;
  chomp $dbDate;
  print <<END
<H3>Statistics</H3>

<TABLE cellpadding=2 cellspacing=2>
<TR style="background-color: lightgrey;"><TD>Phyla</TD><TD align="right">$nPhyla</TD>
<TR style="background-color: lightyellow;"><TD>Classes</TD><TD align="right">$nClasses</TD>
<TR style="background-color: lightgrey;"><TD>Orders</TD><TD align="right">$nOrders</TD>
<TR style="background-color: lightyellow;"><TD>Families</TD><TD align="right">$nFamilies</TD>
<TR style="background-color: lightgrey;"><TD>Genomes</TD><TD align="right">$nGenomes</TD>
END
    ;
  my $buildLog = "../data/build.log";
  if (-e $buildLog) {
    my $nProteins = `egrep '^nProteins' $buildLog | cut -f 2`;
    chomp $nProteins;
    $nProteins = commify($nProteins);
    my $nGenes = `egrep '^nGenes' $buildLog | cut -f 2`;
    $nGenes = commify($nGenes);
    print <<END
<TR style="background-color: lightyellow;"><TD>Protein sequences</TD><TD align="right">$nProteins</TD>
<TR style="background-color: lightgrey;"><TD>Genes</TD><TD align="right">$nGenes</TD>
END
      ;
  }
  print <<END
<TR style="background-color: lightyellow;"><TD>Database date</TD><TD align="right">$dbDate</TD></TR>
</TABLE>


<H3>Downloads</H3>

<UL>
<LI><A HREF="downloadGenomes.cgi">Genomes</A> (tab-delimited)
<LI><A HREF="../data/neighbor.faa.gz">Protein sequences</A> (fasta format, gzipped)
<LI><A HREF="../data/neighbor.db">SQLite3 database</A> (see <A HREF="../lib/neighbor.sql">schema</A>)
<LI><A HREF="https://github.com/morgannprice/fast.genomics">Source code</A> (github repository)
</UL>
END
    ;
  finish_page();
}

if (defined $query{genes}) {
  my $genes = $query{genes};
  # show table of genes
  my $gid = $genes->[0]{gid};
  my @gids = map $_->{gid}, @$genes;
  my %gid = map { $_ => 1 } @gids;
  if (scalar(keys %gid) == 1) {
    my $genome = gidToGenome($gid) || die $gid;
    die $gid unless $genome;
    print p("Found", scalar(@$genes), " matches in $genome->{gtdbSpecies} $genome->{strain} ($genome->{gid})");
    foreach my $gene (@$genes) {
      my $locusTag = $gene->{locusTag};
      print p(a({href => "gene.cgi?locus=$locusTag"}, $locusTag),
              $gene->{proteinId}, $gene->{desc});
    }
    my $query2 = $query; $query2 =~ s/^\S+\s+//; $query2 =~ s/\s+$//;
    my $curatedURL = "https://papers.genomics.lbl.gov/cgi-bin/genomeSearch.cgi?gdb=NCBI"
      . "&gid=${gid}&query=" . uri_escape($query2);
    print p("Or search for", encode_entities($query2),
            "in $genome->{gtdbSpecies} $genome->{strain} using",
            a({-href =>  $curatedURL }, "Curated BLAST"));
  } else {
    print p("Found", scalar(@$genes), " matches");
    foreach my $gene (@$genes) {
      my $locusTag = $gene->{locusTag};
      my $genome = gidToGenome($gene->{gid}) || die $gid;
      print p(a({href => "gene.cgi?locus=$locusTag"}, $locusTag),
              $gene->{proteinId}, $gene->{desc},
              "from", i($genome->{gtdbSpecies}), $genome->{strain}, small("(".$genome->{gid}.")"));
    }
  }
  finish_page();
}

finish_page() unless defined $query{seq};

my $seq = $query{seq};
$seq =~ s/[*]//g;
$seq =~ m/^[a-zA-Z]+$/ || die("Sorry, input sequence has invalid characters");
my $seqDesc = $query{seqDesc}
  || length($seq) . " a.a. beginning with " . substr($seq, 0, 10);

print p("Sequence name:", encode_entities($seqDesc));
print showSequence($seqDesc, $seq);

my $seqDescE = uri_escape($seqDesc);
if (hasMMSeqsHits($seq)) {
  print p("See",
          join(", or ",
               a({-href => "neighbors.cgi?seqDesc=${seqDescE}&seq=${seq}"},
                 "gene neighborhoods"),
               a({-href => "hitTaxa.cgi?seqDesc=${seqDescE}&seq=${seq}"},
                 "taxonomic distribution"),
               a({-href => "downloadHomologs.cgi?seqDesc=${seqDescE}&seq=${seq}",
                  -title => "tab-delimited table of homologs"},
                 "download homologs"),
               a({-href => "compare.cgi?seqDesc=${seqDescE}&seq=${seq}",
                  -title => "compare presence/absence of homologs and their proximity"},
                 "compare presence/absence")));
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

