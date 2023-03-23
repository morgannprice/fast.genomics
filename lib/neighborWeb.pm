# Utilities for the CGI code for the gene neighbor viewer
# Many of these routines assume that .. is the main directory
package neighborWeb;
require Exporter;
use strict;
use CGI;
use DBI;
use Time::HiRes qw{gettimeofday tv_interval};
use Digest::MD5 qw{md5_hex};
use URI::Escape;
use HTML::Entities;

# from the PaperBLAST code base
use pbweb qw{GetMotd commify
             VIMSSToFasta RefSeqToFasta UniProtToFasta FBrowseToFasta pdbToFasta
             runTimerHTML runWhileCommenting};
use pbutils qw{NewerThan ReadFastaEntry};
use neighbor;
use clusterProteins;
use clusterBLASTp;

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw{getDbHandle getSubDbHandle
             parseGeneQuery hasMMSeqsHits getMMSeqsHits
             parseGeneQuery
             getGeneSeqDesc geneSeqDescSeqOptions geneSeqDescSeqHidden
             hitsToGenes hitsToTopGenes hitsToRandomGenes
             start_page finish_page
             locusTagToGene getNearbyGenes gidToGenome
             showSequence formatFastaHtml proteinAnalysisLinks
             clusterGenes
             domainHtml
             getTaxa taxLevels taxToParts taxLevelToParent taxLevelToChild
             capitalize
             geneHitsToProteinAlignment geneHitsToTree
             getSubDbHomologs};

my $memDbh = undef;
sub getDbHandle() {
  if (!defined $memDbh) {
    my $sqldb = "../data/neighbor.db";
    $memDbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;
  }
  return $memDbh;
}

my $memSubDbh = undef;
my $memSubDb = undef;
sub getSubDbHandle($) {
  my ($subDb) = @_;
  die "Invalid subdb: $subDb\n" unless $subDb =~ m/^[a-zA-Z0-9_.-]+$/;
  if (defined $memSubDbh) {
    die "New subDb $subDb vs. $memSubDb" unless $memSubDb eq $subDb;
  } else {
    my $sqldb = "../data/$subDb/sub.db";
    $memSubDbh = DBI->connect("dbi:SQLite:dbname=$sqldb","","",{ RaiseError => 1 }) || die $DBI::errstr;
    $memSubDb = $subDb;
  }
  return $memSubDbh;
}

sub getHitsFile($) {
  my ($seq)= @_;
  $seq = uc($seq);
  die "Invalid sequence for getHitsFile: $seq" unless $seq =~ m/^[A-Z]+$/;
  my $md5 = md5_hex($seq);
  return "../tmp/hits/$md5.hits";
}

sub getMMSeqsDb {
  return "../data/neighbor.sliced";
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

    my $md5 = md5_hex($seq);
    my $faaFile = "../tmp/$md5.$$.faa";
    open(my $fhFaa, ">", $faaFile) || die "Cannot write to $faaFile";
    print $fhFaa ">query\n$seq\n";
    close($fhFaa) || die "Error writing to $faaFile";
    my $tmpOut = "$hitsFile.$$";
    my $procId = $$;
    my $startTime = [gettimeofday()];
    my $timestamp = join(".",@$startTime);
    # higher sensitivity setting for shorter queries
    my $mmseqsSens = length($seq) <= 150 ? 7 : 6;
    my ($nGenomes) = getDbHandle()->selectrow_array("SELECT COUNT(*) FROM Genome");
    my $maxSeqs = int(1.5 * $nGenomes + 0.5);
    my $cmd = "../bin/searchSliced.pl -in $faaFile -sliced $mmseqsDb -out $tmpOut"
      . " -max-seq $maxSeqs -limit $nGenomes -db-load-mode 2 -s $mmseqsSens";
    print CGI::p("Searching for similar proteins with mmseqs2",
                 CGI::small(runTimerHTML())), "\n";
    runWhileCommenting($cmd) == 0
      || die "Error running $cmd -- $!";
    rename($tmpOut, $hitsFile) || die "Renaming $tmpOut to $hitsFile failed";
    unlink($faaFile);
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

    # or in the protein table
    my $genes = $dbh->selectall_arrayref("SELECT * from Gene WHERE proteinId = ?",
                                         { Slice => {} }, $query);
    return ('genes' => $genes) if @$genes > 0;

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
    $seq =~ s/\r//g; # sometimes present due to cut/paste
    $seq =~ s/\s*$//;
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
    my $curatedURL = "https://papers.genomics.lbl.gov/cgi-bin/genomeSearch.cgi?gdb=NCBI"
      . "&gid=$genome->{gid}&query=" . uri_escape($wordQuery);
    return ("genome" => $genome,
            "error" => "Sorry, no genes in $genome->{gtdbSpecies} $genome->{strain} ($genome->{gid}) match "
            . encode_entities($wordQuery)
            . ". Try " . CGI::a({-href => $curatedURL}, "Curated BLAST"))
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

sub hitsToRandomGenes($$$) {
  my ($hits, $n, $minBits) = @_;
  return [] unless $n >= 1;
  my @hits = @$hits;
  # Always keep the top hit
  my $top = shift @hits;
  @hits = grep $_->{bits} >= $minBits, @hits;
  if (scalar(@hits) > $n-1) {
    srand(01262023);
    my @indices = 0..(scalar(@hits) - 1);
    my @rand = map rand, @indices;
    @indices = sort { $rand[$a] <=> $rand[$b] } @indices;
    splice @indices, $n-1;
    @hits = map $hits[$_], @indices;
  }
  unshift @hits, $top;
  @hits = sort { $b->{bits} <=> $a->{bits} } @hits;
  return hitsToTopGenes(\@hits, $n);
}

sub TopDivHtml($$) {
  my ($banner, $URL) = @_;
  return <<END
<div style="background-color: lightgreen; display: block; position: absolute; left:-1px; top:0px;
  width: 100%; padding: 0em; z-index: 400;">
<P style="margin: 0.1em; margin-left: 0.4em; font-size: 160%;">
<A HREF="$URL" style="color: #660066; font-family: 'Montserrat', sans-serif; text-decoration: none;">
$banner
</A></P></div>
<P style="margin: 0em;">&nbsp;</P>
END
;
}

sub start_page {
  my (%param) = @_;
  my $title = $param{title} || "";
  my ($nGenomes) = getDbHandle()->selectrow_array("SELECT COUNT(*) FROM Genome");
  $nGenomes = commify($nGenomes);
  my $banner = $param{banner} || "<i>fast.genomics</i> &ndash; "
    . qq{<SPAN style="font-size:smaller;"> compare $nGenomes genera of bacteria and archaea</SPAN>};
  my $bannerURL = $param{bannerURL} || 'search.cgi';
  print
    CGI::header(-charset => 'utf-8'),
    CGI::start_html(-head => CGI::Link({-rel => "shortcut icon", -href => "../static/favicon.ico"}),
               -title => $title eq "" ? "fast.genomics" : $title),
    TopDivHtml($banner, $bannerURL),
    pbweb::GetMotd(),
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

# Assumes that genes are not too enormous; fetches ~50 kb on each side
# Returns a reference to a list of rows from the Gene table
sub getNearbyGenes($) {
  my ($gene) = @_;
  my $mid = int( ($gene->{begin} + $gene->{end})/2 );
  return getDbHandle()->selectall_arrayref(
    qq[ SELECT * FROM Gene WHERE gid = ? AND scaffoldId = ? AND begin >= ? AND end <= ? ORDER BY begin ],
    { Slice => {} },
    $gene->{gid}, $gene->{scaffoldId}, $mid - 50*1000, $mid + 50*1000);
}

# returns a list of HTML strings
# If $genome is undef, just uses gram negative models for psortb
sub proteinAnalysisLinks($$$) {
  my ($header, $seq, $genome) = @_;
  my $psortType = "negative";
  if (defined $genome) {
    $psortType = "archaea" if $genome->{gtdbDomain} eq "Archaea";
    $psortType = "positive" if $genome->{gtdbPhylum} =~ m/^Firmicutes|Actinobacteriota/;
  }
  return pbweb::analysisLinks('desc' => $header,
                              'seq' => $seq,
                              'skip' => {'fast.genomics' => 1},
                              'psortType' => $psortType,
                              'fbLoad' => 1);
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
    my ($seq) = getDbHandle()->selectrow_array(
        "SELECT sequence FROM Protein WHERE proteinId = ?",
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

sub geneSeqDescSeqOptions($$$) {
  my ($gene, $seqDesc, $seq) = @_;
  return "locus=$gene->{locusTag}" if defined $gene;
  my $seqDescE = encode_entities($seqDesc);
  return "seqDesc=$seqDescE&seq=$seq";
}

sub geneSeqDescSeqHidden($$$) {
  my ($gene, $seqDesc, $seq) = @_;
  return qq[<INPUT type="hidden" name="locus" value="$gene->{locusTag}">]
    if defined $gene;
  my $seqDescE = encode_entities($seqDesc);
  return qq[<INPUT type="hidden" name="seqDesc" value="$seqDescE">]
    . qq[<INPUT type="hidden" name="seq" value="$seq">];
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
  # RECOMPUTE environment variable is for testing
  if (-e $clusterFile && ! $ENV{RECOMPUTE}) {
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
    print qq{<P id="ClusterInfo">Clustering the proteins for coloring...</P>}, "\n";
    foreach my $proteinId (keys %proteinToGenes) {
      my ($seq) = getDbHandle()->selectrow_array("SELECT sequence FROM Protein WHERE proteinId = ?",
                                                 {}, $proteinId);
      die "Unknown proteinId $proteinId" unless $seq;
      $proteinSeq{$proteinId} = $seq;
    }
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
    print qq{<SCRIPT>document.getElementById("ClusterInfo").innerHTML = "";</SCRIPT>}, "\n";
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

# Given a domain value, or a reference to a hash that has the gtdbDomain field,
# return html for a blue B or green A
sub domainHtml($) {
  my ($domain) = @_;
  $domain = $domain->{gtdbDomain} || die
    if ref $domain;
  my $domainChar = $domain eq "Bacteria" ? "B" : "A";
  my $domainColor = $domainChar eq "B" ? "blue" : "green";
  my $URL = "taxon.cgi?level=domain&taxon=$domain";
  return qq[<a title="$domain" style="color: ${domainColor}; text-decoration: none;" href="$URL">$domainChar</a>];
}

# Returns a hash of level => taxon => row from Taxon table
sub getTaxa {
  my $taxa = getDbHandle()->selectall_arrayref("SELECT * from Taxon", { Slice => {} });
  # assume level x taxon is unique
  my %taxa = (); # level => taxon => row
  foreach my $tax (@$taxa) {
    $taxa{ $tax->{level} }{ $tax->{taxon} } = $tax;
  }
  return \%taxa;
}

my @levelsOrdered = qw{domain phylum class order family genus species};
my %levelToParent = map { $levelsOrdered[$_] => $levelsOrdered[$_-1] } (1..6);
my %levelToChild = map { $levelToParent{$_} => $_ } (keys %levelToParent);

sub taxLevels() { return @levelsOrdered };

sub taxLevelToParent($) {
  my ($level) = @_;
  return undef if $level eq "domain";
  die "Unknown level $level" unless exists $levelToParent{$level};
  return $levelToParent{$level};
}

sub taxLevelToChild($) {
  my ($level) = @_;
  return undef if $level eq "species";
  die "Unknown level $level" unless exists $levelToChild{$level};
  return $levelToChild{$level};
}

# Given a row from the Taxon table and the taxa hash (from getTaxa),
# returns a reference to hash of level => taxon for this taxon and all its parents
sub taxToParts($$) {
  my ($tax, $taxa) = @_;
  my %out = ();
  for(;;) {
    $out{ $tax->{level} } = $tax->{taxon};
    if ($tax->{parent} eq "") {
      last;
    } else {
      my $level = $levelToParent{ $tax->{level} } || die $tax->{level};
      my $parent = $tax->{parent};
      die "No taxon for $parent at level $level"
        unless exists $taxa->{$level}{$parent};
      $tax = $taxa->{$level}{$parent};
    }
  }
  return \%out;
}

sub capitalize($) {
  my ($string) = @_;
  return $string if !defined $string || $string eq "";
  return uc(substr($string, 0, 1)) . lc(substr($string, 1));
}

# Returns a hash of protSpec to gene, where each prot spec is proteinId___sBegin___sEnd
sub geneHitsToProtSpec {
  my ($genes) = @_;
  my %protSpecToGenes;
  foreach my $gene (@$genes) {
    die "Not protein-coding: " . $gene->{locusTag}
      unless $gene->{proteinId} ne "";
    die "Empty alignment range for " . $gene->{locusTag}
      unless $gene->{sEnd} > $gene->{sBegin};
    my $protSpec = join("___", $gene->{proteinId}, $gene->{sBegin}, $gene->{sEnd});
    push @{ $protSpecToGenes{$protSpec} }, $gene;
  }
  return %protSpecToGenes;
}

# Given a list of gene hits, which contain locusTag, proteinId, sBegin, and sEnd,
# adds an alignedSeq entry to each gene hit.
# No return value.
sub geneHitsToProteinAlignment {
  my ($genes) = @_;
  my %protSpecToGenes = geneHitsToProtSpec($genes);
  my @protSpec = sort keys %protSpecToGenes;
  my $md5 = md5_hex(join(",", @protSpec));
  my $n = scalar(@protSpec);
  my $alnFile = "../tmp/hits/${md5}_${n}.aln";
  # RECOMPUTE environment variable is for testing
  if (-e $alnFile && ! $ENV{RECOMPUTE}) {
    # reuse the alignment
  } else {
    my $muscle = "../bin/muscle3";
    die "No such executable: $muscle\n" unless -x $muscle;
    my $tmpFile = "$alnFile.$$.tmp";
    my $faaFile = "/tmp/neighborWeb.$$.faa";
    open(my $fhFaa, ">", $faaFile) || die "Cannot write to $faaFile";
    foreach my $protSpec (@protSpec) {
      my ($proteinId, $sBegin, $sEnd) = split /___/, $protSpec;
      die unless $proteinId ne "";
      my ($seq) = getDbHandle()->selectrow_array("SELECT sequence FROM Protein WHERE proteinId = ?",
                                                 {}, $proteinId);
      die "No sequence for $proteinId" unless defined $seq;
      die "Invalid begin:end $sBegin : $sEnd for $proteinId"
        unless $sBegin >=1 && $sEnd <= length($seq);
      my $subseq = substr($seq, $sBegin-1, $sEnd-$sBegin+1);
      print $fhFaa ">$protSpec\n$subseq\n";
    }
    close($fhFaa) || die "Error writing to $faaFile\n";
    my $muscleOptions = "-maxiters 2 -maxmb 1000";
    print CGI::p(CGI::small("Running muscle")), "\n";
    my $cmd = "$muscle -quiet $muscleOptions < $faaFile > $tmpFile";
    system($cmd) == 0 || die "$cmd\nfailed: $!";
    rename($tmpFile, $alnFile) || die "Cannot rename to $tmpFile to $alnFile";
    unlink($faaFile);
  }

  # read the alignment
  open(my $fh, "<", $alnFile) || die "Cannot read $alnFile";
  my $state = {};
  while (my ($protSpec, $seq) = ReadFastaEntry($fh, $state)) {
    die "Unexpected alignment id: $protSpec" unless exists $protSpecToGenes{$protSpec};
    foreach my $gene (@{ $protSpecToGenes{$protSpec} }) {
      $gene->{alignedSeq} = $seq;
    }
  }
  close($fh) || die "Error reading $alnFile";
  foreach my $gene (@$genes) {
    die "No aligned sequence for " . $gene->{locusTag}
      unless $gene->{alignedSeq};
  }
  return 1;
}

# Given gene hits, infer a tree. The returned tree has locus tags as the node identifiers,
# so each locus tag must be present just once.
sub geneHitsToTree {
  my ($genes) = @_;
  my %protSpecToGenes = geneHitsToProtSpec($genes);
  my @protSpec = sort keys %protSpecToGenes;
  my $md5 = md5_hex(join(",", @protSpec));
  my $n = scalar(@protSpec);
  my $treeFile = "../tmp/hits/${md5}_${n}.tree";
  # RECOMPUTE environment variable is for testing
  if (-e $treeFile && ! $ENV{RECOMPUTE}) {
    # reuse the tree
  } else {
    # ensure that each gene shows up just once
    my %seen = ();
    foreach my $gene (@$genes) {
      die "Duplicate locus tag: " . $gene->{locusTag} . "\n"
        if exists $seen{ $gene->{locusTag} };
      $seen{ $gene->{locusTag} } = 1;
      die "Invalid locus tag " . $gene->{locusTag}
        unless $gene->{locusTag} =~ m/^[0-9A-Za-z_.-]+$/;
    }
    my $fastTree = "../bin/FastTree";
    die "No such executable: $fastTree\n" unless -x $fastTree;
    geneHitsToProteinAlignment($genes);
    my $alnFile = "/tmp/neighborWeb.$$.aln";
    open(my $fhAln, ">", $alnFile) || die "Cannot write to $alnFile\n";
    foreach my $gene (@$genes) {
      print $fhAln ">" . $gene->{locusTag} . "\n" . $gene->{alignedSeq} . "\n";
    }
    close($fhAln) || die "Error writing to $alnFile\n";
    my $tmpFile = "$treeFile.tmp";
    my $cmd = "$fastTree -quiet < $alnFile > $tmpFile";
    print CGI::p(CGI::small("Running FastTree"));
    system($cmd) == 0 || die "$cmd\nfailed: $!\n";
    rename($tmpFile, $treeFile) || die "Cannot rename to $tmpFile to $treeFile";
    unlink($alnFile);
  }
  return MOTree::new( file => $treeFile );
}

sub computeSubDbHomologs($$) {
  my ($subDb, $seq) = @_;
  die "Invalid sequence" unless $seq =~ m/^[A-Z]+$/;
  my $subDir = "../data/$subDb";
  my $clusterDb = "$subDir/cluster.faa";
  die "No such file: $clusterDb.pin" unless -e "$clusterDb.pin";

  my $blastall = "../bin/blast/blastall";
  die "No such executable: $blastall" unless -x $blastall;
  my $formatdb = "../bin/blast/formatdb";
  die "No such executable: $formatdb" unless -x $formatdb;

  my $tmpPre = "/tmp/neighborWeb.$$";
  print CGI::p("Running BLASTp");
  my $hits1 = runBLASTp('blastall' => $blastall, 'query' => $seq, 'db' => $clusterDb,
                        'eValue' => 1e-3, 'nCPUs' => 10);
  return [] if @$hits1 == 0;

  print CGI::p("Found", scalar(@$hits1), " clustered hits"), "\n";
  print CGI::p("Processing the hits and re-running BLASTp"), "\n";
  my %clusterIds = map { $_->{subject} => 1 } @$hits1;
  my @clusterIds = sort keys %clusterIds;
  my $subDbh = getSubDbHandle($subDb);
  my $proteinIds = expandByClusters(\@clusterIds, $subDbh);
  print CGI::p("Expanded to", scalar(@$proteinIds), "proteins"),"\n";
  proteinsToFaa($proteinIds, "$tmpPre.faa", $subDbh);
  formatBLASTp($formatdb, "$tmpPre.faa");
  my ($dbSize) = $subDbh->selectrow_array("SELECT nClusteredAA FROM ClusteringInfo");
  my $hits = runBLASTp('blastall' => $blastall, 'query' => $seq, 'db' => "$tmpPre.faa",
                       'eValue' => 1e-3, 'nCPUs' => 2, 'dbSize' => $dbSize);
  unlink("$tmpPre.faa");
  foreach my $suffix (qw{phr pin psq}) {
    unlink("$tmpPre.faa.$suffix");
  }
  print CGI::p("Returning",scalar(@$hits),"hits"),"\n";
  return $hits;
}

# Given a subdb and a sequence, compute the homologs. Returns a reference to a list of hashes,
# each containg subject, eValue, etc.
sub getSubDbHomologs($$) {
  my ($subDb, $seq) = @_;
  $seq = uc($seq);
  my $subDir = "../data/$subDb";
  my $clusterDb = "$subDir/cluster.faa";
  die "No such file: $clusterDb.pin" unless -e "$clusterDb.pin";
  my $md5 = md5_hex($seq);
  my $md5 = md5_hex($seq);
  my $hitsFile = "../tmp/hits/${subDb}_$md5.hits";

  if (NewerThan($hitsFile, $clusterDb.".pin")) {
    return parseBLASTpHits($hitsFile);
  } else {
    my $hits = computeSubDbHomologs($subDb, $seq);
    saveBLASTpHits($hits, $hitsFile);
    return $hits;
  }
}
1;
