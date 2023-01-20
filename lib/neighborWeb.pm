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

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw{getDbHandle parseGeneQuery getMMSeqsHits parseGeneQuery hitsToGenes
             start_page finish_page
             showSequence};

my $dbh = undef;
sub getDbHandle() {
  if (!defined $dbh) {
    my $sqldb = "../data/neighbor.db";
    $dbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;
  }
  return $dbh;
}

# From a protein sequence to a reference list of hits (each a hash),
# as parsed by readMMSeqsHits.
# Caches the results in ../tmp/hits/md5.hits
# Uses ../data/tmp as the 
sub getMMSeqsHits($) {
  my ($seq) = @_;
  $seq = uc($seq);
  die "Invalid sequence for getMMSeqsHits: $seq" unless $seq =~ m/^[A-Z]+$/;
  my $md5 = md5_hex($seq);
  my $hitsFile = "../tmp/hits/$md5.hits";
  my $mmseqsDb = "../data/neighbor.mmseqs";
  unless (NewerThan($hitsFile, $mmseqsDb)) {
    die "No such file: $mmseqsDb" unless -e $mmseqsDb;
    my $mmseqs = "../bin/mmseqs";
    die "No such executable: $mmseqs" unless -x $mmseqs;

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
    print CGI::p("Searching for similar proteins with mmseqs2")."\n";
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
#   gene -- if a single gene in the database was found
#   genes -- if multiple genes matched
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
    my $gene = $dbh->selectrow_hashref("SELECT * FROM Gene WHERE locusTag = ? OR locusTag = ?",
                                       {}, $query, uc($query));
    if (defined $gene) {
      # seq will be undef if not protein-coding
      my $seq;
      if ($gene->{proteinId} ne "") {
        ($seq) = $dbh->selectrow_array("SELECT sequence FROM Protein WHERE proteinId = ?",
                                       {}, $gene->{proteinId});
        die "Could not find sequence for protein $gene->{proteinId}"
          unless $seq;
      }
      my $seqDesc = $gene->{locusTag} . " " . $gene->{desc};
      return ('gene' => $gene, 'seq' => $seq,
              'seqDesc' => $gene->{locusTag} . " " . $gene->{desc} );
    }
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
                                           AND (desc LIKE ? OR desc LIKE ? OR desc LIKE ?)
                                           LIMIT 200 },
                                         { Slice => {} },
                                         $genome->{gid},
                                         "${wordQuery}%", "%-${wordQuery}%", "% ${wordQuery}%");
    return ("error" => "Sorry, no genes in $genome->{gtdbSpecies} $genome->{strain} ($genome->{gid}) match "
            . encode_entities($wordQuery))
      if @$genes == 0;
    if (@$genes == 1) {
      # return the single gene
      my $gene = $genes->[0];
      my $seq;
      if ($gene->{proteinId} ne "") {
        ($seq) = $dbh->selectrow_array("SELECT sequence FROM Protein WHERE proteinId = ?",
                                       {}, $gene->{proteinId});
        die "Could not find sequence for protein $gene->{proteinId}"
          unless $seq;
      }
      my $seqDesc = $gene->{locusTag} . " " . $gene->{desc};
      return ('gene' => $gene, 'seq' => $seq,
              'seqDesc' => $gene->{locusTag} . " " . $gene->{desc} );
    }
    return ('genes' => $genes);
  }
  return ('error'=> "Sorry, did not understand the query.");
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

sub TopDivHtml($$) {
  my ($banner, $URL) = @_;
  return <<END
<div style="background-color: #40C0CB; display: block; position: absolute; top: 0px; left: -1px;
  width: 100%; padding: 0.25em; z-index: 400;">
<H2 style="margin: 0em;">
<A HREF="$URL" style="color: gold; font-family: 'Montserrat', sans-serif; font-style:italic;
  text-shadow: 1px 1px 1px #000000; text-decoration: none;">
$banner
</A></H2></div>
<P style="margin: 0em;">&nbsp;</P>
END
;
}

sub start_page {
  my (%param) = @_;
  my $title = $param{title} || "FastComper";
  my ($nGenomes) = getDbHandle()->selectrow_array("SELECT COUNT(*) FROM Genome");
  my $banner = $param{banner} || "FastComper &mdash;"
    . qq{<SPAN style="font-size:smaller;"> browse }
    . commify($nGenomes) . " genera of bacteria and archaea</SPAN>";
  my $bannerURL = $param{bannerURL} || 'search.cgi';
  print
    CGI::header(-charset => 'utf-8'),
    CGI::start_html(-head => CGI::Link({-rel => "shortcut icon", -href => "../static/favicon.ico"}),
               -title => $title),
    TopDivHtml($banner, $bannerURL),
    CGI::h2($title),
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
  my @seqPieces = $seq =~ /.{1,60}/g;
  $out .=  CGI::p({ -id => "querySequence",
            -style => "font-family: monospace; display:none; font-size:90%; padding-left: 2em;" },
          join(CGI::br(), ">".encode_entities($seqDesc), @seqPieces));
  return $out;
}

1;
