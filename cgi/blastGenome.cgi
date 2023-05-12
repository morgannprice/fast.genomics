#!/usr/bin/perl -w
use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use IO::Handle; # for autoflush
use DBI;
use HTML::Entities;
use URI::Escape;
use List::Util qw{min max};
use lib "../lib";
use neighbor;
# neighborWeb.pm relies on various PaperBLAST libraries
use lib "../../PaperBLAST/lib";
use neighborWeb;
use pbweb qw{commify};

# required CGI arguments:
# gid -- which genome to show information about
# optional arguments:
# order (which subdb to use)
# query -- usually an identifier or locus tag,
#   or a genus and word(s) from a protein description.
#   See neighborWeb::parseGeneQuery for details.
#
# Optional arguments:
# format = tblastn -- show the raw tblastn output

my $cgi = CGI->new;
setOrder(param('order'));
my $query = $cgi->param('query') || "";
my $gid = $cgi->param('gid') || die "Must specify gid\n";
my $genome = gidToGenome($gid) || die "Unknown gid\n";
my $format = $cgi->param('format') || "";
die "Invalid format\n" unless $format eq "" || $format eq "tblastn";

start_page('title' => "BLAST against a genome");
autoflush STDOUT 1; # show preliminary results
print p("Compare to",
        a({-href => addOrderToURL("genome.cgi?gid=$gid")},
          $genome->{gtdbGenus}, $genome->{strain}));

my %query = parseGeneQuery($query);
my $queryGene;

if (defined $query{genes}) {
  my $genes = $query{genes};
  if (scalar(@$genes) > 1) {
    my $URL = addOrderToURL("genes.cgi?" . join("&", map "g=" . $_->{locusTag}, @$genes));
    print
      p("Sorry, more than one gene matched", encode_entities($query)),
      p("See", a{-href => $URL}, "table");
    finish_page();
  } else {
    $queryGene = $genes->[0];
  }
} elsif (!defined $query{seq}) {
  print p(b($query{error})) if $query{error};
  print
    start_form(-name => 'input', -method => 'GET', -action => 'blastGenome.cgi'),
    orderToHidden(),
    qq{<INPUT type="hidden" name="gid" value="$gid">},
    p("Enter an identifier from UniProt, RefSeq, PDB, or MicrobesOnline,",
      br(),
      "or a protein sequence in FASTA or Uniprot format,",
      br(),
      textarea( -name  => 'query', -value => '', -cols  => 70, -rows  => 10 ),
      br(),
      br(),
      submit('Search'), reset()),
    end_form;
  finish_page();
}

my ($seq, $seqDesc);
if (defined $queryGene) {
  if ($queryGene->{proteinId} eq "") {
    print p("Sorry,",
            a({-href => addOrderToURL("gene.cgi?locus=".$queryGene->{locusTag})},
              $queryGene->{locusTag}),
            "is not a protein-coding gene");
    finish_page();
  }
  $seqDesc = $queryGene->{locusTag} . " " . $queryGene->{desc};
  ($seq) = getDbHandle->selectrow_array("SELECT sequence FROM Protein WHERE proteinid = ?",
                                        {}, $queryGene->{proteinId});
  die "No sequence for $queryGene->{locusTag}" unless $seq;
} else {
  $seq = $query{seq};
  $seq =~ s/[*]//g;
  $seq =~ m/^[a-zA-Z]+$/ || die("Sorry, input sequence has invalid characters");
  $seqDesc = $query{seqDesc}
    || length($seq) . " a.a. beginning with " . substr($seq, 0, 10);
}

print
  p("Searching for homologs of", encode_entities($seqDesc),
    br(),
    "in", a({-href => addOrderToURL("genome.cgi?gid=$gid")}, $genome->{gtdbGenus}, $genome->{strain})),
  showSequence($seqDesc, $seq);

my $blastall = "../bin/blast/blastall";
my $formatdb = "../bin/blast/formatdb";
foreach my $x ($blastall, $formatdb) {
  die "No such executable: $x\n" unless -x $x;
}

my $tmpPre = "/tmp/blastGenome.$$";
open(my $fhDb, ">", "$tmpPre.db") || die "Cannot write to $tmpPre.db";
my $protGenes = getDbHandle()->selectall_arrayref(
    "SELECT locusTag, sequence FROM Gene JOIN Protein USING (proteinId)"
    . " WHERE gid = ? AND proteinId <> ''",
   {}, $gid);
foreach my $row (@$protGenes) {
  my ($locusTag, $seq) = @$row;
  print $fhDb ">$locusTag\n$seq\n";
}
close($fhDb) || die "Error writing to $tmpPre.db";

open(my $fhFaa, ">", "$tmpPre.faa") || die "Cannot write too $tmpPre.faa";
print $fhFaa ">query\n$seq\n";
close($fhFaa) || die "Error writing to $tmpPre.faa";

my @aaHits = ();
if ($format eq "") {
  system("$formatdb -p T -i $tmpPre.db") == 0 || die "formatdb on $tmpPre.db failed: $!";
  system(qq{$blastall -p blastp -F "m S" -i $tmpPre.faa -d $tmpPre.db -e 1e-3 -m 8 -o $tmpPre.hits}) == 0
    || die "blastall failed: $!";
  open(my $fhHits, "<", "$tmpPre.hits") || die "Cannot read $tmpPre.hits\n";
  foreach my $line (<$fhHits>) {
    chomp $line;
    my @F = split /\t/, $line;
    push @aaHits, \@F;
  }
  close($fhHits) || die "Error reading $tmpPre.hits";

  print p("Found", scalar(@aaHits), "hits to annotated proteins. Or try",
          a({-href => addOrderToURL("blastGenome.cgi?gid=$gid")},
            "another search")."."), "\n";
}

my $qLen = length($seq);
my @geneHits = (); # the hit genes, and the alignment extents (saved in sBegin/sEnd/qBegin/qEnd)
if (@aaHits > 0) {
  print qq{<TABLE cellpadding=2 cellspacing=2>},
    Tr(th([ "locus",  "description",
            a({-title => '%identity'}, "id."),
            a({-title => '%coverage of query'}, "cov.") ])),
    "\n";
  my $iRow = 0;
  foreach my $hit (@aaHits) {
    my (undef, $locusTag, $identity, $alen, $mm, $gap, $qBeg, $qEnd, $sBeg, $sEnd, $evalue, $bits) = @$hit;
    my $genes = getDbHandle()->selectall_arrayref("SELECT * FROM Gene JOIN Protein USING (proteinId)"
                                                  . " WHERE locusTag = ?",
                                                  { Slice => {} }, $locusTag);
    die "Cannot find protien for $locusTag" unless @$genes == 1;
    my ($geneHit) = @$genes;
    $geneHit->{sBegin} = $sBeg;
    $geneHit->{sEnd} = $sEnd;
    $geneHit->{qBegin} = $qBeg;
    $geneHit->{qEnd} = $qEnd;
    push @geneHits, $geneHit;
    my $sLen = length($geneHit->{sequence});
    print Tr({-bgcolor => $iRow % 2 == 1 ? "white" : "lightgrey"},
             td({-valign => "top", -align => "left" },
                a({-href => addOrderToURL("gene.cgi?locus=".$locusTag)}, $locusTag)),
             td({-valign => "top", -align => "left" }, encode_entities($geneHit->{desc})),
             td({-valign => "top", -align => "right" },
                a({-title => "$bits bits, E = $evalue"}, int(0.5 + $identity)."%")),
             td({-valign => "top", -align => "right" },
                a({-title => "$qBeg:$qEnd/$qLen of query is similar to $sBeg:$sEnd/$sLen of $locusTag",
                   -href => addOrderToURL("alignPair.cgi?seq=$seq&seqDesc=".uri_escape($seqDesc)
                                          . "&locus2=$locusTag"),
                   -style => "text-decoration: none;"},
                  int(0.5 + 100.0*($qEnd-$qBeg+1)/$qLen)."%"))
            );
    $iRow++;
  }
  print qq{</TABLE>}, "\n";
}

my $fnaFile = "../ind/refseq_${gid}.fna";
if (!-e $fnaFile) {
  print p("Not searching for hits againts the genome: could not find the nucleotide sequence");
} else {
  system("cp", $fnaFile, "$tmpPre.fna") == 0 || die "cp failed";
  system("$formatdb -p F -i $tmpPre.fna") == 0 || die "formatdb on $tmpPre.fna failed: $!";
  if ($format eq "tblastn") {
    print p("Searching the genome with tblastn. Or try",
            a({ -href => addOrderToURL("blastGenome.cgi?gid=$gid") },
              "another search")."."), "\n";
  } else {
    print p("Searching the genome for unannotated homologs using tblastn"), "\n";
  }

  my $blastCmd = qq{$blastall -p tblastn -F "m S" -i $tmpPre.faa -d $tmpPre.fna -e 1e-3 -o $tmpPre.hits};
  $blastCmd .= " -m 8" unless $format eq "tblastn";
  system($blastCmd) == 0
  || die "blastall failed: $!\n$blastCmd\n";

  if ($format eq "tblastn") { # raw tblastn output
    print "<pre>\n";
    system("cat", "$tmpPre.hits");
    print "</pre>\n";
  } else { # regular format
    # filter new hits from the tblastn output and show a table
    my @ntHits = ();
    open(my $fhHits, "<", "$tmpPre.hits") || die "Cannot read $tmpPre.hits\n";
    foreach my $line (<$fhHits>) {
      chomp $line;
      my @F = split /\t/, $line;
      push @ntHits, \@F;
    }
    close($fhHits) || die "Error reading $tmpPre.hits";

    # Convert the protein hits into approximate positions, so it is easy to check if the blastx hits are redundant
    my @ntOnlyHits = (); # not redundant with protein-conding hits
    foreach my $ntHit (@ntHits) {
      my (undef, $scaffoldId, $identity, $alen, $mm, $gap, $qBeg, $qEnd, $sBeg, $sEnd, $evalue, $bits) = @$ntHit;
      my $strand = $sBeg < $sEnd ? "+" : "-";
      my $overlap = 0;
      foreach my $geneHit (@geneHits) {
        if ($geneHit->{scaffoldId} eq $scaffoldId
            && $geneHit->{strand} eq $strand) {
          # Compute the approximate coordinates of the gene hit
          my ($geneBeg, $geneEnd); # always beg < end
          if ($strand eq "+") {
            $geneBeg = $geneHit->{begin} + 3 * ($geneHit->{sBegin}-1);
            $geneEnd = $geneHit->{begin} + 3 * ($geneHit->{sEnd}-1);
          } else {
            $geneBeg = $geneHit->{end} - 3 * ($geneHit->{sEnd}-1);
            $geneEnd = $geneHit->{end} - 3 * ($geneHit->{sBegin}-1);
          }
          # xBeg < xEnd
          my $xBeg = min($sBeg,$sEnd);
          my $xEnd = max($sBeg,$sEnd);
          my $overlapBeg = max($xBeg, $geneBeg);
          my $overlapEnd = min($xEnd, $geneEnd);
          # Ignore hits that are almost entirely covered by the protein hit (allow 10 a.a. of wobble)
          if ($overlapBeg < $overlapEnd && $overlapEnd-$overlapBeg+1 >= $xEnd-$xBeg+1 - 10*3) {
            $overlap = 1;
            last;
          }
        }
      }
      push @ntOnlyHits, $ntHit unless $overlap;
    }

    my $tblastnLink = addOrderToURL("blastGenome.cgi?format=tblastn&gid=$gid&query="
                                    . uri_escape($query));
    print p("Found", scalar(@ntHits),  " hits to the genome, of which",
            scalar(@ntHits) - scalar(@ntOnlyHits), "are redundant with hits to annotated proteins.",
            @ntHits > 0 ?
            small("See", a({-href => $tblastnLink}, "raw tblastn output") . ".")
            : "");

    if (@ntOnlyHits > 0) {
      print "<TABLE cellpadding=2 cellspacing=2>",
        Tr(th(["scaffold", "strand", "region",
               a({-title => '%identity'}, "id."),
               a({-title => '%coverage of query'}, "cov."),
               "overlapping genes"])), "\n";
      my $iRow = 0;
      my %hitLocusTag = map { $_->{locusTag} => $_ } @geneHits;
      foreach my $ntHit (@ntOnlyHits) {
        my (undef, $scaffoldId, $identity, $alen, $mm, $gap, $qBeg, $qEnd, $sBeg, $sEnd, $evalue, $bits) = @$ntHit;
        my $strand = $sBeg < $sEnd ? "+" : "-";
        my $overlapping = getDbHandle()->selectall_arrayref(
          "SELECT * from Gene WHERE gid = ? AND scaffoldId = ? AND begin <= ? AND end >= ?"
          . " ORDER BY begin" . ($strand eq "+" ? "" : " DESC"),
          { Slice => {} },
          $gid, $scaffoldId, max($sBeg, $sEnd), min($sBeg, $sEnd));
        my @showOverlap = ();
        foreach my $o (@$overlapping) {
          my $locus = $o->{locusTag};
          my $showOverlap = a({-href => addOrderToURL("gene.cgi?locus=$locus"),
                               -title => encode_entities($o->{desc})
                               . ($o->{proteinId} eq "" ? " (not protein-coding)" : ""),
                               -style => "text-decoration: none;"},
                              $locus);
          if (exists $hitLocusTag{$locus}) {
            my $aaHit = $hitLocusTag{$locus};
            $showOverlap .= " " . small(a({-title => "tblastn aligned to $qBeg:$qEnd of query"
                                           . " instead of $aaHit->{qBegin}:$aaHit->{qEnd}"},
                                          "above"));
          }
          push @showOverlap, $showOverlap;
        }
        my $w = max(500, abs($sBeg-$sEnd)/2);
        my $sBeg2 = max(1, min($sBeg,$sEnd) - $w);
        my $sEnd2 = max($sBeg,$sEnd) + $w;
        print Tr({-bgcolor => $iRow % 2 == 1 ? "white" : "lightgrey"},
                 td({-valign => "top", -align => "left"}, $scaffoldId),
                 td({-valign => "top", -align => "center"}, $strand),
                 td({-valign => "top", -align => "left"},
                    a({-href => "https://www.ncbi.nlm.nih.gov/nuccore/$scaffoldId?report=graph"
                       . "&from=$sBeg2&to=$sEnd2"
                       . " &mk=$sBeg:$sEnd|hit_region|00008f",
                       -title => "NCBI's viewer",
                       -style => "text-decoration: none;"},
                      commify(min($sBeg,$sEnd)) . " : " . commify(max($sBeg,$sEnd)))),
                 td({-valign => "top", -align => "right"},
                    a({-title => "$bits bits, E = $evalue"}, int($identity)."%")),
                 td({-valign => "top", -align => "right"},
                    a({-title => "covers $qBeg:$qEnd/$qLen of query"},
                      int(0.5 + 100.0*($qEnd-$qBeg+1)/$qLen)."%")),
                 td({-valign => "top", -align => "left"},
                    @showOverlap > 0 ? join(", ", @showOverlap) : "&nbsp;"));
        print "\n";
        $iRow++;
      }
      print "</TABLE>\n";
    }
  } # end if regular format
} # end if have nt sequence

unlink("$tmpPre.faa");
unlink("$tmpPre.db");
unlink("$tmpPre.hits");
unlink("$tmpPre.fna");
foreach my $suffix (qw{nhr nin nsq}) {
  unlink("$tmpPre.fna.$suffix");
}
foreach my $suffix (qw{phr pin psq}) {
  unlink("$tmpPre.db.$suffix");
}
finish_page();
