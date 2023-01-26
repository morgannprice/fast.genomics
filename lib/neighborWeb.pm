# Utilities for the CGI code for the gene neighbor viewer
# Many of these routines assume that .. is the main directory
package neighborWeb;
require Exporter;
use strict;
use CGI;
use DBI;
use Time::HiRes qw{gettimeofday tv_interval};
use Digest::MD5 qw{md5_hex};
use HTML::Entities;

# from the PaperBLAST code base
use pbweb qw{GetMotd commify
             VIMSSToFasta RefSeqToFasta UniProtToFasta FBrowseToFasta pdbToFasta};
use pbutils qw{NewerThan};
use neighbor;
use clusterProteins;

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw{getDbHandle parseGeneQuery hasMMSeqsHits getMMSeqsHits parseGeneQuery
             hitsToGenes hitsToTopGenes
             start_page finish_page
             locusTagToGene getNearbyGenes gidToGenome
             showSequence formatFastaHtml proteinAnalysisLinks
             getGeneSeqDesc
             clusterGenes};

my $dbh = undef;
sub getDbHandle() {
  if (!defined $dbh) {
    my $sqldb = "../data/neighbor.db";
    $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;
  }
  return $dbh;
}

sub getHitsFile($) {
  my ($seq)= @_;
  $seq = uc($seq);
  die "Invalid sequence for getHitsFile: $seq" unless $seq =~ m/^[A-Z]+$/;
  my $md5 = md5_hex($seq);
  return "../tmp/hits/$md5.hits";
}

sub getMMSeqsDb {
  return "../data/neighbor.mmseqs";
}

sub hasMMSeqsHits($) {
  my ($seq) = @_;
  return NewerThan(getHitsFile($seq), getMMSeqsDb());
}

# From a protein sequence to a reference list of hits (each a hash),
# as parsed by readMMSeqsHits.
# Caches the results in ../tmp/hits/md5.hits
# Uses a temporary subdirectory of ../data/tmp as the temporary directory for mmseqs2
sub getMMSeqsHits($) {
  my ($seq) = @_;
  my $mmseqsDb = getMMSeqsDb();
  my $hitsFile = getHitsFile($seq);
  unless (NewerThan($hitsFile, $mmseqsDb)) {
    die "No such file: $mmseqsDb" unless -e $mmseqsDb;
    my $mmseqs = "../bin/mmseqs";
    die "No such executable: $mmseqs" unless -x $mmseqs;

    my $md5 = md5_hex($seq);
    my $faaFile = "../tmp/$md5.$$.faa";
    open(my $fhFaa, ">", $faaFile) || die "Cannot write to $faaFile";
    print $fhFaa ">query\n$seq\n";
    close($fhFaa) || die "Error writing to $faaFile";
    my $tmpOut = "$hitsFile.$$";
    die "No such directory: ../data/tmp" unless -d "../data/tmp";
    my $procId = $$;
    my $startTime = [gettimeofday()];
    my $timestamp = join(".",@$startTime);
    my $mmseqsTmpDir = "../data/tmp/mm.$procId.$timestamp";
    mkdir($mmseqsTmpDir) || die "Cannot mkdir $mmseqsTmpDir";
    my ($nGenomes) = getDbHandle()->selectrow_array("SELECT COUNT(*) FROM Genome");
    my @cmd = ($mmseqs, "easy-search", $faaFile, $mmseqsDb, $tmpOut, $mmseqsTmpDir,
               "--max-seqs", $nGenomes,
               "--threads", 1,
               "--db-load-mode", 2,
               ">", "$tmpOut.log");
    print CGI::p("Searching for similar proteins with mmseqs2..."), "\n";
    system(join(" ",@cmd)) == 0
      || die "Error running @cmd -- $!";
    rename($tmpOut, $hitsFile) || die "Renaming $tmpOut to $hitsFile failed";
    unlink($faaFile);
    unlink("$tmpOut.log");
    system("rm -Rf $mmseqsTmpDir");
    my $elapsed = tv_interval($startTime);
    print CGI::p(sprintf("mmseqs2 finished in %.1f seconds", $elapsed))."\n";
  }
  return readMMSeqsHits($hitsFile);
}

# The query may be a locus tag in the database,
#   or, an identifier from UniProt, RefSeq, MicrobesOnline, Fitness Browser, PDB,
#	i.e. METE_CORGL Q8NRB3 WP_012018426.1 F605_RS0102535 VIMSS14484 b1110 3osd 3osdA 3osd_1
#   or, a query of the form genus search_words,
#	where the words must be in the protein description in the specified order,
#	and % is the wild card character
#   or a sequence in fasta format (header line optional) or uniprot format
# Returns a hash with the fields
#   genes -- if gene(s) matched in the database
#   seq -- if a single gene or known identifier was found
#   seqDesc -- the description
#   error -- set if there was an error (and with entities encoded)
# On an empty query, it will return an empty hash.
# Otherwise, either genes or seq will be set. If seq is set, gene and
# seqDesc may be set as well. If a non-protein-coding gene was
# searched for, then only gene will be set.

sub parseGeneQuery($) {
  my ($query) = @_;
  # Trim leading and trailing whitespace
  $query =~ s/^\s+//;
  $query =~ s/\s+$//;
  return () if $query eq "";
  my $dbh = getDbHandle();

  if ($query =~ m/^[a-zA-Z0-9_.-]+$/ && $query =~ m/[0-9_]/) {
    # try to find the identifier
    # first, look in the SQL database
    my $gene = locusTagToGene($query);
    return ('genes' => [$gene]) if defined $gene;
    #else
    # Check MicrobesOnline and fitness browser first because they are fast
    # Check pdb next because it is very specific about which identifiers to look for
    # Check uniprot before refseq  because it is faster
    my $fasta = VIMSSToFasta($query) # MicrobesOnline
      || FBrowseToFasta("../fbrowse_data", $query)
      || pdbToFasta($query)
      || UniProtToFasta($query)
      || RefSeqToFasta($query);
    return("error" => "Sorry, could not find protein " .  encode_entities($query))
      unless $fasta;
    $query = $fasta; # handled below
  }

  if ($query =~ m/^>/) {
    # fasta format
    my @lines = split /\n/, $query;
    fail("Sorry, query has a FASTA header but no sequence") unless @lines > 0;
    my $header = shift @lines;
    my $seqDesc = $header;
    $seqDesc =~ s/^>//;
    my $seq = join("", @lines);
    return ('seq' => $seq, 'seqDesc' => $seqDesc);
  }

  if ($query =~ m/\n/ || ($query =~ m/^[A-Z*]+$/ && length($query) >= 10)) {
    # multi-line protein sequence, possibly in uniprot format (leading digits)
    # or a sequence on a line by itself
    my $seq = $query;
    $seq =~ s/\s//g;
    $seq =~ s/\d//g;
    $seq =~ s/[*]//g; # stop codons are sometimes rendered as *
    return ("error" => "Sorry, could not interpret the query")
      unless $seq =~ m/^[A-Za-z]+$/;
    return ('seq' => $seq);
  }

  if ($query =~ m/ /) {
    # Should be a query of the form genus-name description_words
    my @parts = split / /, $query;
    my $genus = shift @parts;
    # use LIKE so matching is not case sensitive
    my $genomes = $dbh->selectall_arrayref("SELECT * FROM Genome WHERE gtdbGenus LIKE ?",
                                           { Slice => {} }, $genus);
    return ("error" => "Text queries must start with a genus, but genus "
            . encode_entities($genus) . " is not in the database")
      if @$genomes == 0;
    return ("error" => "More than one genome matches " . encode_entities($genus))
      if @$genomes > 1;
    my $genome = $genomes->[0];
    my $wordQuery = join(" ", @parts);
    my $genes = $dbh->selectall_arrayref(qq{ SELECT * FROM Gene
                                           WHERE gid=?
                                           AND (desc LIKE ? OR desc LIKE ? OR desc LIKE ? OR desc LIKE ?)
                                           LIMIT 200 },
                                         { Slice => {} },
                                         $genome->{gid},
                                         "${wordQuery}%", "%-${wordQuery}%",
                                         "% ${wordQuery}%", "% (${wordQuery}%");
    return ("error" => "Sorry, no genes in $genome->{gtdbSpecies} $genome->{strain} ($genome->{gid}) match "
            . encode_entities($wordQuery))
      if @$genes == 0;
    return ('genes' => $genes);
  }
  return ('error'=> "Sorry, did not understand the query.");
}

# Case insensitive if the locus tag is all upper case
sub locusTagToGene($) {
  my ($locusTag) = @_;
  return getDbHandle()->selectrow_hashref("SELECT * FROM Gene WHERE locusTag = ? OR locusTag = ?",
                                          {}, $locusTag, uc($locusTag));
}

sub gidToGenome($) {
  my ($gid) = @_;
  return getDbHandle()->selectrow_hashref("SELECT * FROM Genome WHERE gid = ?",
                                          {}, $gid);
}

# Convert hits to proteinIds (from readMMSeqsHits) to a list of genes
# Both input and output should be a list of references; the output rows
# hves all of the fields from the Gene table as well.
sub hitsToGenes($) {
  my ($hits) = @_;
  my $dbh = getDbHandle();
  my @hitGenes = ();
  foreach my $hit (@$hits) {
    my $proteinId = $hit->{subject};
    my $hg = $dbh->selectall_arrayref("SELECT * FROM Gene WHERE proteinId=?",
                                      { Slice => {} }, $proteinId);
    die "No genes for protein $proteinId" unless @$hg > 0;
    foreach my $hitGene (@$hg) {
      while (my ($key,$value) = each %$hit) {
        $hitGene->{$key} = $value;
      }
      push @hitGenes, $hitGene;
    }
  }
  return \@hitGenes;
}

# Like hitsToGenes, but returns a maximum of $n entries
sub hitsToTopGenes($$) {
  my ($hits, $n) = @_;
  return [] unless $n >= 1;
  my @top = @$hits;
  $#top = $n-1 if scalar(@top) > $n;
  my $hitGenes = hitsToGenes(\@top);
  splice $hitGenes, $n if scalar(@$hitGenes) > $n;
  return $hitGenes;
}

sub TopDivHtml($$) {
  my ($banner, $URL) = @_;
  return <<END
<div style="background-color: #40C0CB; display: block; position: absolute; top: 0px; left: -1px;
  width: 100%; padding: 0.25em; z-index: 400;">
<H2 style="margin: 0em;">
<A HREF="$URL" style="color: grey90; font-family: 'Montserrat', sans-serif; text-decoration: none;">
$banner
</A></H2></div>
<P style="margin: 0em;">&nbsp;</P>
END
;
}

sub start_page {
  my (%param) = @_;
  my $title = $param{title} || "";
  my ($nGenomes) = getDbHandle()->selectrow_array("SELECT COUNT(*) FROM Genome");
  my $banner = $param{banner} || "<i>FastComper</i> &ndash; "
    . qq{<SPAN style="font-size:smaller;"> browse }
    . commify($nGenomes) . " representative genomes</SPAN>";
  my $bannerURL = $param{bannerURL} || 'search.cgi';
  print
    CGI::header(-charset => 'utf-8'),
    CGI::start_html(-head => CGI::Link({-rel => "shortcut icon", -href => "../static/favicon.ico"}),
               -title => $title eq "" ? "FastComper" : "FastComper $title"),
    TopDivHtml($banner, $bannerURL),
    $title ? CGI::h2($title) : "",
    "\n";
}

sub finish_page() {
  print <<END
<P>
<small>
<center>by <A HREF="http://morgannprice.org/">Morgan Price</A>,
<A HREF="http://genomics.lbl.gov/">Arkin group</A><BR>
Lawrence Berkeley National Laboratory
</center>
</small>
</P>
END
;
  print CGI::end_html();
  exit(0);
}

# Returns the HTML, including javascript, to show a "Show sequence" link that
# replaces itself with the fasta sequence.
# Because the ids are hardcoded in javascript, This can only be used once per page
sub showSequence($$) {
  my ($seqDesc, $seq) = @_;

  my $out = <<END
<SCRIPT>
showQuery = function() {
 document.getElementById("showSequenceLink").style.display = "none";
 document.getElementById("querySequence").style.display = "block";
 return false;
}
</SCRIPT>
END
;

  $out .= CGI::p({ -id => "showSequenceLink", -style => "font-size:90%" },
           CGI::a({ -href => "#", -onclick => "return showQuery()" }, "Show sequence"));
  $out .=  CGI::p({ -id => "querySequence",
                    -style => "font-family: monospace; display:none; font-size:90%; padding-left: 2em;" },
                  formatFastaHtml($seqDesc, $seq));
  return $out;
}

sub formatFastaHtml($$) {
  my ($seqDesc, $seq) = @_;
  my @seqPieces = $seq =~ /.{1,60}/g;
  return CGI::span({ -style => "font-family: monospace;" },
                   join(CGI::br(), ">".encode_entities($seqDesc), @seqPieces));
}

# Assumes that genes are not too enormous; fetches ~40 kb on each side
# Returns a reference to a list of rows from the Gene table
sub getNearbyGenes($) {
  my ($gene) = @_;
  my $mid = int( ($gene->{begin} + $gene->{end})/2 );
  return getDbHandle()->selectall_arrayref(
    qq[ SELECT * FROM Gene WHERE gid = ? AND scaffoldId = ? AND begin >= ? AND end <= ? ORDER BY begin ],
    { Slice => {} },
    $gene->{gid}, $gene->{scaffoldId}, $mid - 40*1000, $mid + 40*1000);
}

# returns a list of HTML strings
# If $genome is undef, just uses gram negative models for psortb
sub proteinAnalysisLinks($$$) {
  my ($header, $seq, $genome) = @_;
  $header = encode_entities($header);

  my ($psortType, $psortShow) = ("negative", "Gram-negative bacteria");
  if (defined $genome) {
    ($psortType, $psortShow) = ("positive", "Gram-positive bacteria")
      if $genome->{gtdbPhylum} =~ m/^Firmicutes|Actinobacteriota/;
    ($psortType, $psortShow) = ("archaea", "archaea")
      if $genome->{gtdbDomain} eq "Archaea";
  }

  my $newline = "%0A";
  my $fitnessBlastJs = <<END
<SCRIPT src="https://fit.genomics.lbl.gov/d3js/d3.min.js"></SCRIPT>
<SCRIPT src="https://fit.genomics.lbl.gov/images/fitblast.js"></SCRIPT>
<SCRIPT>
var server_root = "https://fit.genomics.lbl.gov/";
var seq = "$seq";
fitblast_load_short("fitblast_short", server_root, seq);
</SCRIPT>
END
    ;
  return
    ( CGI::a({-href => "http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?query=>${header}$newline$seq"},
             "PaperBLAST") .
      " (search for papers about homologs of this protein)",
      qq[<A TITLE="Fitness BLAST compares your sequence to bacterial proteins that have mutant phenotypes"
            NAME="#fitness">Fitness BLAST:</A>
         <SPAN ID="fitblast_short"><SMALL>loading...</SMALL></SPAN>]
      . $fitnessBlastJs,
      CGI::a({-href => "http://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi?seqinput=>${header}$newline$seq"},
             "Search the Conserved Domains Database"),
      CGI::a({ -href => "http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/FindSequence.pl?pasted=$seq",
               -title => "Find similar proteins with known structures (PDBsum)"},
             "Search protein structures"),
      "Predict protein localization: " .
      CGI::a({-href => "https://papers.genomics.lbl.gov/cgi-bin/psortb.cgi?name=${header}&type=${psortType}&seq=${seq}",
              -title => "PSORTb v3.0 for $psortShow"},
        "PSORTb") . " ($psortShow)",
      "Find homologs in the " .
      CGI::a({-href => "https://iseq.lbl.gov/genomes/seqsearch?sequence=>${header}$newline$seq"},
             "ENIGMA genome browser")
    );
}

# Given a CGI object, which may specify locus or seq and seqDesc,
# extract gene, seq, and seqDesc
# (If a locus tag is specified, seq and seqDesc are from the database, and seq will be undef
#  if it is not a protein-coding gene)
#
sub getGeneSeqDesc($) {
  my ($cgi) = @_;
  my $locus = $cgi->param('locus');
  if (defined $locus && $locus ne "") {
    my $gene = locusTagToGene($locus) || die "Unknown locus tag in locus parameter";
    my ($seq) = $dbh->selectrow_array("SELECT sequence FROM Protein WHERE proteinId = ?",
                                      { Slice => {} }, $gene->{proteinId})
      if $gene->{proteinId} ne "";
    return ($gene, $gene->{locusTag} . " " . $gene->{desc}, $seq);
  }
  #else
  my $seqDesc = $cgi->param('seqDesc');
  die "Must specify seqDesc if no locus" unless defined $seqDesc && $seqDesc ne "";
  my $seq = $cgi->param('seq') || die "Must specify seq if no locus";
  $seq =~ m/^[A-Z]+$/ || die "Invalid sequence";
  return (undef, $seqDesc, $seq);
}

# Given a list of objects, fetch their protein sequences and cluster them using lastal.
# Returns a list of clusters, each of which is a list of gene objects.
# Non-coding genes are never put in a cluster.
# Caches the clustering
sub clusterGenes {
  my (@genes) = @_;
  @genes = grep $_->{proteinId} ne "", @genes;
  my %proteinToGenes = ();
  foreach my $gene (@genes) {
    push @{ $proteinToGenes{$gene->{proteinId}} }, $gene;
  }
  my @proteinIds = sort keys %proteinToGenes;

  my $md5 = md5_hex(join(",", @proteinIds));
  my $n = scalar(@proteinIds);
  my $clusterFile = "../tmp/hits/${n}_${md5}.cluster";
  my @proteinClusters;
  if (-e $clusterFile) {
    open (my $fh, "<", $clusterFile) || die "Cannot read $clusterFile";
    while (my $line = <$fh>) {
      chomp $line;
      my @F = split /\t/, $line;
      foreach my $proteinId (@F) {
        die "Unknown protein $proteinId in $clusterFile"
          unless exists $proteinToGenes{$proteinId};
      }
      push @proteinClusters, \@F;
    }
    close($fh) || die "Error reading $clusterFile";
  } else {
    my %proteinSeq = ();
    foreach my $proteinId (keys %proteinToGenes) {
      my ($seq) = $dbh->selectrow_array("SELECT sequence FROM Protein WHERE proteinId = ?",
                                        {}, $proteinId);
      die "Unknown proteinId $proteinId" unless $seq;
      $proteinSeq{$proteinId} = $seq;
    }
    print CGI::p("Running LAST to cluster and color the proteins...")."\n";
    @proteinClusters = @{ clusterProteins(\%proteinSeq, 0.5) };
    foreach my $proteinCluster (@proteinClusters) {
      foreach my $proteinId (@$proteinCluster) {
        die "Unknown protein $proteinId from clusterProteins()"
          unless exists $proteinToGenes{$proteinId};
      }
    }
    my $tmpFile = "$clusterFile.$$.tmp";
    open (my $fh, ">", $tmpFile) || die "Cannot write to $tmpFile";
    foreach my $proteinCluster (@proteinClusters) {
      print $fh join("\t", @$proteinCluster)."\n";
    }
    close($fh) || die "Error writing to $tmpFile";
    rename($tmpFile, $clusterFile) || die "rename $tmpFile to $clusterFile failed";
  }

  # convert protein clusters to gene clusters
  my @geneClusters = ();
  foreach my $cluster (@proteinClusters) {
    my @geneCluster = ();
    foreach my $proteinId (@$cluster) {
      die unless defined $proteinToGenes{$proteinId};
      push @geneCluster, @{ $proteinToGenes{$proteinId} };
    }
    push @geneClusters, \@geneCluster;
  }
  return \@geneClusters;
}

1;
