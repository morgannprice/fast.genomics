#!/usr/bin/perl -w
use strict;

use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use IO::Handle; # for autoflush
use HTML::Entities;
use URI::Escape;
use List::Util qw{min max};
use lib "../lib";
use neighbor;
use genesSvg;
use binom qw{hyperLogProbTail};
# neighborWeb.pm relies on various PaperBLAST libraries
use lib "../../PaperBLAST/lib";
use neighborWeb;
use pbweb qw{commify};

# CGI arguments:
# locus (a locus tag in the database) or seqDesc and seq
# either query2, locus2, or seq2 and seqDesc2
# optional: format=tsv

my $cgi = CGI->new;
my $tsv = $cgi->param('format') eq 'tsv';
my ($gene, $seqDesc, $seq) = getGeneSeqDesc($cgi);
my $options = geneSeqDescSeqOptions($gene,$seqDesc,$seq); # for 1st gene
my $hidden = geneSeqDescSeqHidden($gene,$seqDesc,$seq); # for 1st gene

unless ($tsv) {
  start_page(title => 'Compare gene presence/absence');
  autoflush STDOUT 1; # show preliminary results
}

# Show first gene and check that it is protein coding and has homologs
if ($tsv) {
  die "first gene is not protein-coding\n"
    unless defined $seq;
 die "No homologs for first gene yet\n" unless hasMMSeqsHits($seq);
}
else {
  print h3("First gene");
  if ($gene) {
    my $genome = gidToGenome($gene->{gid});
    print p(a({-href => "gene.cgi?$options"}, $gene->{locusTag}) . ":",
            encode_entities($gene->{desc}),
            br(),
            "from", i($genome->{gtdbSpecies}), encode_entities($genome->{strain}));
  } else {
    print p(a({-href => "seq.cgi?$options"}, encode_entities($seqDesc)));
  }
  if (!defined $seq) {
    print p("Sorry, first gene is not protein-coding");
    finish_page();
  }
  unless (hasMMSeqsHits($seq)) {
    print p("Sorry, first gene does not have homologs yet");
    print a({-href => "findHomologs.cgi?$options"}, "Find homologs for the first gene");
    finish_page();
  }
}
my $hits1 = getMMSeqsHits($seq);
if (scalar(@$hits1) == 0) {
  exit(0) if $tsv;
  print p("Sorry, no homologs were found for the first gene");
  finish_page;
}

# Handle the second gene options
my $locus2 = $cgi->param('locus2');
my $seq2 = $cgi->param('seq2');
my $seqDesc2 = $cgi->param('desc2');
my $query2 = $cgi->param('query2') || "";
my $gene2;
print h3("Second gene") unless $tsv;
if (defined $locus2 && $locus2 ne "") {
  $gene2 = locusTagToGene($locus2) || die "Unknown locus tag in locus2 parameter";
} elsif ($seq2) {
  $seq2 =~ m/^[A-Z]+$/ || die "Invalid sequence in seq2";
} else {
  die "query2 option is not supported with tsv\n" if $tsv;
  # parse query, or show the form
  $query2 =~ s/^\s+//;
  $query2 =~ s/\s+$//;
  my %query2 = parseGeneQuery($query2);
  if (!defined $query2{genes} && !defined $query2{seq}) {
    print p(b($query2{error})) if $query2{error};
    print
      start_form(-name => 'input', -method => 'GET', -action => 'compare.cgi'),
      $hidden,
      p("Enter an identifier from UniProt, RefSeq, PDB, or MicrobesOnline,",
        br(),
        "or a protein sequence in FASTA or Uniprot format,",
        br(),
        "or a genus name and a protein description",
        br(),
        textarea( -name  => 'query2', -value => '', -cols  => 70, -rows  => 10 ),
        br(),
        br(),
        submit('Search'), reset()),
      end_form;
    finish_page();
  }
  #else successfuly parsed query
  if (defined $query2{genes}) {
    my $genes = $query2{genes};
    if (scalar(@$genes) > 1) {
      my $URL = "genes.cgi?" . join("&", map "g=" . $_->{locusTag}, @$genes);
      print
        p("Sorry, more than one gene matched", encode_entities($query2)),
        p("See", a{-href => $URL}, "table");
      finish_page();
    } else {
      $gene2 = $genes->[0];
    }
  } elsif (defined $query2{seq}) {
    $seq2 = $query2{seq};
    $seqDesc2 = $query2{seqDesc};
  }
}

if (defined $gene2) {
  $seqDesc2 = $gene2->{locusTag} . " " . $gene2->{desc};
  ($seq2) = getDbHandle()->selectrow_array("SELECT sequence FROM Protein WHERE proteinId = ?",
                                           { Slice => {} }, $gene2->{proteinId})
    if $gene2->{proteinId} ne "";
} elsif (!defined $seqDesc2) {
  $seqDesc2 = length($seq2) . " a.a. beginning with " . substr($seq2, 0, 10);
}
my $options2 = geneSeqDescSeqOptions($gene2,$seqDesc2,$seq2);
if ($tsv) {
  die "Second gene is not protein-coding\n" unless $seq2;
} else {
  if ($gene2) {
    my $genome = gidToGenome($gene2->{gid});
    print p(a({-href => "gene.cgi?$options2"}, $gene2->{locusTag}) . ":",
            encode_entities($gene2->{desc}),
            br(),
            "from", i($genome->{gtdbSpecies}), encode_entities($genome->{strain}));
  } else {
    print p(a({-href => "seq.cgi?$options2"}, encode_entities($seqDesc2)));
  }
  if (!defined $seq2) {
    print p("Sorry, second gene is not protein-coding");
    finish_page();
  }
}
unless (hasMMSeqsHits($seq2)) {
  die "No homologs for second gene yet\n" if $tsv;
  print p("Sorry, second gene does not have homologs yet");
  print a({-href => "findHomologs.cgi?$options2"}, "Find homologs for the second gene");
  finish_page();
}
my $hits2 = getMMSeqsHits($seq2);
if (scalar(@$hits2) == 0) {
  exit(0) if $tsv;
  print p("Sorry, no homologs were found for the second gene");
  finish_page;
}

if ($tsv) {
  print "Content-Type:text/tab-separated-values\n";
  print "Content-Disposition: attachment; filename=compareHomologs.tsv\n\n";
}

# Compare hits1 and hits2
my $geneHits1 = hitsToGenes($hits1);
my $max1 = estimateTopScore($geneHits1->[0], $seq);
my $geneHits2 = hitsToGenes($hits2);
my $max2 = estimateTopScore($geneHits2->[0], $seq2);
print p("Loaded", commify(scalar(@$geneHits1)), "homologs for the first gene",
        "and", commify(scalar(@$geneHits2)), "homologs for the second gene"), "\n"
  unless $tsv;

# Top hits by genome
my %byg1 = ();
foreach my $hit (@$geneHits1) {
  my $gid = $hit->{gid};
  $byg1{$gid} = $hit unless exists $byg1{$gid};
}
my %byg2 = ();
foreach my $hit (@$geneHits2) {
  my $gid = $hit->{gid};
  $byg2{$gid} = $hit unless exists $byg2{$gid};
}

# Compute ranks
foreach my $i (0..(scalar(@$geneHits1)-1)) {
  $geneHits1->[$i]{rank} = $i+1;
}
foreach my $i (0..(scalar(@$geneHits2)-1)) {
  $geneHits2->[$i]{rank} = $i+1;
}

# For each genome that contains a hit to either,
# the two scores, whether they are nearby, whether they are the same gene
my %genomeList = ();
foreach my $gid (keys %byg1) {
  $genomeList{$gid} = [];
}
foreach my $gid (keys %byg2) {
  $genomeList{$gid} = [];
}

my $genomes = getDbHandle()->selectall_hashref("SELECT * FROM Genome", "gid", { Slice => {} });
my $nGenomes = scalar(keys %$genomes);

my $closeKb = 5;
my $closeNt = $closeKb * 1000;
# gid with a hit for either geonme to hash of
# gid, inBoth, same (gene), close, ratio1, ratio2 (0 if no hit), hit1, hit2
my %gidScores = ();
foreach my $gid (keys %genomeList) {
  my $hit1 =  exists $byg1{$gid} ? $byg1{$gid} : undef;
  my $hit2 =  exists $byg2{$gid} ? $byg2{$gid} : undef;
  my $inBoth = defined $hit1 && defined $hit2;
  my $sameGene = $inBoth && $hit1->{locusTag} eq $hit2->{locusTag};
  my $close = $inBoth && ! $sameGene
    && $hit1->{scaffoldId} eq $hit2->{scaffoldId}
    && $hit1->{strand} eq $hit2->{strand}
    && (abs($hit1->{begin} - $hit2->{end}) <= $closeNt
        || abs($hit1->{end} - $hit2->{begin}) <= $closeNt);
  $gidScores{$gid} = { 'gid' => $gid,
                       'inBoth' => $inBoth, 'same' => $sameGene, 'close' => $close,
                       'ratio1' => defined $hit1 ? $hit1->{bits} / $max1 : 0,
                       'ratio2' => defined $hit2 ? $hit2->{bits} / $max2 : 0,
                       'rank1' => defined $hit1 ? $hit1->{rank} : "",
                       'rank2' => defined $hit2 ? $hit2->{rank} : "",
                       'order' => min(defined $hit1 ? $hit1->{rank} : 2*$nGenomes,
                                      defined $hit2 ? $hit2->{rank} : 2*$nGenomes) };
}

# Sort by the better rank
my @gidValues = sort { $a->{order} <=> $b->{order} } values %gidScores;
if ($tsv) {
  my @genomeFields = qw{gtdbDomain gtdbPhylum gtdbClass gtdbOrder gtdbFamily
                        gtdbGenus gtdbSpecies strain};
  print join("\t", "assemblyId",
             @genomeFields,
             "locusTag1", "locusTag2",
             "rank1", "rank2",
             "bits1", "bits2",
             "ratio1", "ratio2",
             "identity1", "identity2",
             "alnLength1", "alnLength2",
             "close")."\n";
  foreach my $row (@gidValues) {
    my $gid = $row->{gid};
    my $tag1 = $byg1{$gid}{locusTag} || "";
    my $tag2 = $byg2{$gid}{locusTag} || "";
    my @genomeOut = map $genomes->{$gid}{$_}, @genomeFields;
    print join("\t", $gid, @genomeOut,
               $tag1, $tag2,
               $row->{rank1}, $row->{rank2},
               $byg1{$gid}{bits} || "", $byg2{$gid}{bits} || "",
               $row->{ratio1}, $row->{ratio2},
               $byg1{$gid}{identity} || "", $byg2{$gid}{identity} || "",
               $byg1{$gid}{alnLength} || "", $byg2{$gid}{alnLength} || "",
               $row->{close} ? 1 : 0)."\n";
  }
  exit(0);
}

#else HTML report

my $nInBoth = scalar(grep $_->{inBoth}, @gidValues);
my $n1 = scalar(grep $_->{ratio1} > 0, @gidValues);
my $n2 = scalar(grep $_->{ratio2} > 0, @gidValues);
my $nClose = scalar(grep $_->{close}, @gidValues);
my $nSame = scalar(grep $_->{same}, @gidValues);

my $n1good = scalar(grep $_->{ratio1} >= 0.3, @gidValues);
my $n2good = scalar(grep $_->{ratio2} >= 0.3, @gidValues);

my @good = grep $_->{ratio1} >= 0.3 && $_->{ratio2} >= 0.3, @gidValues;
my $nGoodBoth = scalar(@good);
my $nGoodClose = scalar(grep $_->{close}, @good);
my $nGoodSame = scalar(grep $_->{same}, @good);

my $nInBothExpected = $nSame == $nGenomes ? $nGenomes
  : $nSame + (($n1-$nSame) * ($n2-$nSame))/($nGenomes-$nSame);
my $nGoodBothE = $nGoodSame == $nGenomes ? $nGenomes
  : $nGoodSame + (($n1good - $nGoodSame) * ($n2good - $nGoodSame))/($nGenomes-$nGoodSame);

print
  h3("Co-occurence of Homologs"),
  p("Found homologs in", commify($n1), "and", commify($n2), "genomes, respectively.",
    "$nInBoth genomes contain homologs of both genes (versus",
    commify(int($nInBothExpected+0.5)), "expected)."),
  p("Considering only the best hit in each genome,",
        "$nClose hits are nearby (within $closeKb kb) and on the same strand, and $nSame are to the same gene."),
  h3("Co-occurence of Good Homologs"),
  p("Found good homologs (above 30% of maximum bit score) in",
    commify($n1good), "and", commify($n2good), "genomes, respectively.",
    "$nGoodBoth genomes contain good homologs of both genes (versus",
    commify(int($nGoodBothE+0.5)), "expected)."),
  p("Among the best hits in those $nGoodBoth genomes,",
    "$nGoodClose are nearby (within $closeKb kb) and on the same strand, and $nGoodSame are to the same gene.");

# Find the most significant threshold (if any)
# for the enrichment of the two being in the same genome
if ($nInBoth > $nSame + 1) {
  my @gidsLeft = map $_->{gid}, grep ! $_->{same}, @gidValues;
  my $nGenomesLeft = $nGenomes - $nSame;

  # Sort the remaining hits by score, and, as we go, keep track of how many
  # are in the same genome
  my @gids1 = grep exists $byg1{$_}, @gidsLeft;
  # Break ties on bits by evalue (mmseqs rounds the bit scores)
  @gids1 = sort { $byg1{$b}{bits} <=> $byg1{$a}{bits}
                    || $byg1{$a}{eValue} <=> $byg1{$b}{eValue} } @gids1;
  my @gids2 = grep exists $byg2{$_}, @gidsLeft;
  @gids2 = sort { $byg2{$b}{bits} <=> $byg2{$a}{bits}
                    || $byg2{$a}{eValue} <=> $byg2{$b}{eValue} } @gids2;


  my $n = min(scalar(@gids1), scalar(@gids2));
  my %in1 = ();
  my %in2 = ();
  my %inBoth = ();
  my @logP = ();
  for (my $i = 0; $i < $n; $i++) {
    my $gid1 = $gids1[$i];
    my $gid2 = $gids2[$i];
    $in1{$gid1} = 1;
    $in2{$gid2} = 1;
    foreach my $gid ($gid1, $gid2) {
      $inBoth{$gid} = 1 if $in1{$gid} && $in2{$gid};
    }
    # log(P) according to a 1-sided fisher exact test on more-extreme contingency tables
    # is the same as the hypergeometric distribution where we fix top hits for genome 1 as being 
    # present in nTop genomes and absent from nGenomesLeft-nTop genomes;
    # then choose nTop genomes (representing those with genome2), and
    # how often do nBoth or more of those have genome1
    my $nTop = $i+1;
    $logP[$i] = hyperLogProbTail(scalar(keys %inBoth), $nTop, $nGenomesLeft-$nTop, $nTop);
  }
  my $bestLogP = min(@logP);
  my ($bestI) = grep $logP[$_] == $bestLogP, (0..($n-1));
  my $correctedLog10P = ($bestLogP + log($n)) / log(10);
  print h3("Strongest Threshold for Co-occurence");
  if ($correctedLog10P < -2) {
    my %gidsThresh1 = map { $_ => 1 } @gids1[0..$bestI];
    my %gidsThresh2 = map { $_ => 1 } @gids2[0..$bestI];
    my @gidsThresh = grep exists $gidsThresh2{$_}, keys %gidsThresh1;
    my $nCloseThresh = scalar(grep $gidScores{$_}{close}, @gidsThresh);
    my $thresh1 = $byg1{ $gids1[$bestI] }{bits};
    my $f1 = int(0.5 + 100 * $thresh1/$max1);
    my $thresh2 = $byg2{ $gids2[$bestI] }{bits};
    my $f2 = int(0.5 + 100 * $thresh2/$max2);
    my $fClose = int(0.5 + 100 * $nCloseThresh/scalar(@gidsThresh));
    print p("The strongest statistical signal of co-occurence is for the top",
            commify($bestI+1), "homologs of each gene, which co-occur in",
            commify(scalar(@gidsThresh)), "genomes,",
            "P = 10<sup>" . sprintf("%.1f", $correctedLog10P) . "</sup>",
            "(Fisher's exact test, 1-sided, with Bonferonni correction).",
            "The corresponding bit score thresholds are $thresh1 ($f1% of max)",
            " and $thresh2 ($f2% of max), respectively.",
           "$nCloseThresh of these co-occuring homologs ($fClose%) are nearby (within $closeKb kb and on the same strand).");
  } else {
    print p("No significant co-occurence");
  }
}

my $seqDesc2e = encode_entities($seqDesc2);
my $downloadURL = join("&",
                       "compare.cgi?${options}",
                       (defined $gene2 ? "locus2=$gene2->{locusTag}"
                        : "seqDesc2=$seqDesc2e&seq2=$seq2"),
                       "format=tsv");
print
  h3("Download"),
  p("Download a",
        a({ -href => $downloadURL }, "table"),
        "(tab-delimited)",
       "of the best homolog(s) in each genome");
finish_page();
