#!/usr/bin/perl -w
# Compare presence/absence pattern of 2 proteins, and
# how often the best hits are near each other and on the same strand
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
use svgPlot;
use binom qw{hyperLogProbTail};
# neighborWeb.pm relies on various PaperBLAST libraries
use lib "../../PaperBLAST/lib";
use neighborWeb;
use pbweb qw{commify};

# required CGI arguments:
# locus (a locus tag in the database) or seqDesc and seq
# optional CGI arguments:
# order -- which subdb to use
# query2, locus2, or seq2 and seqDesc2 (to specify the 2nd locus)
# format -- use tsv to download a table
# For taxon mode:
# taxLevel -- which level to tabulate at
# taxMode -- either "close" or "both"
# Good -- 1 if considering good hits to both loci only
# all -- set to "all" to show results for all taxa, sorted by taxonomy
#   defaults to "freq" -- show only taxa with hits, and sort by frequency of the gene being present

my $closeKb = 5;
my $closeNt = $closeKb * 1000;

my $cgi = CGI->new;
setOrder(param('order'));
my $tsv = ($cgi->param('format') || "") eq  'tsv';
neighborWeb::setQuietMode() if $tsv;
my ($gene, $seqDesc, $seq) = getGeneSeqDesc($cgi);
my $options = geneSeqDescSeqOptions($gene,$seqDesc,$seq); # for 1st gene
my $hidden = geneSeqDescSeqHidden($gene,$seqDesc,$seq); # for 1st gene

my $taxLevel = $cgi->param('taxLevel');
my $taxMode = $cgi->param('taxMode') || "both";

unless ($tsv) {
  my $title = "Compare gene presence/absence";
  my $order = getOrder();
  $title .= " within $order" if $order ne "";
  $title = ($taxMode eq "close" ? "Taxa with both genes nearby" : "Taxa with both genes")
    if $taxLevel;
  start_page(title => $title);
  autoflush STDOUT 1; # show preliminary results
}

# Show first gene and check that it is protein coding and has homologs
if ($tsv) {
  die "first gene is not protein-coding\n"
    unless defined $seq;
 die "No homologs for first gene yet\n" unless hasHits($seq);
}
else {
  print h3("First protein");
  if ($gene) {
    my $genome = gidToGenome($gene->{gid});
    print p(a({-href => "gene.cgi?$options"}, $gene->{locusTag}) . ":",
            encode_entities($gene->{desc}),
            br(),
            "from",
            a({-href => addOrderToURL("genome.cgi?gid=$genome->{gid}"),
               -style => "text-decoration:none;"},
              i($genome->{gtdbSpecies}), encode_entities($genome->{strain})));
  } else {
    print p(a({-href => "seq.cgi?$options"}, encode_entities($seqDesc)));
  }
  if (!defined $seq) {
    print p("Sorry, first gene is not protein-coding");
    finish_page();
  }
  unless (hasHits($seq)) {
    print p("Sorry, homologs for the first protein have not been computed yet");
    print a({-href => "findHomologs.cgi?$options"}, "Compute homologs for the first protein");
    finish_page();
  }
}
my $hits1 = getHits($seq);
if (scalar(@$hits1) == 0) {
  exit(0) if $tsv;
  print p("Sorry, no homologs were found for the first protein");
  finish_page;
}

# Handle the second gene options
my $locus2 = $cgi->param('locus2');
my $seq2 = $cgi->param('seq2');
my $seqDesc2 = $cgi->param('seqDesc2');
my $query2 = $cgi->param('query2') || "";
my $gene2;
print h3("Second protein") unless $tsv;
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
      my $URL = addOrderToURL("genes.cgi?" . join("&", map "g=" . $_->{locusTag}, @$genes));
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
            "from",
            a({-href => addOrderToURL("genome.cgi?gid=".$genome->{gid}),
               -style => "text-decoration: none;"},
              i($genome->{gtdbSpecies}), encode_entities($genome->{strain})));
  } else {
    print p(a({-href => "seq.cgi?$options2"}, encode_entities($seqDesc2)));
  }
  if (!defined $seq2) {
    print p("Sorry, second gene is not protein-coding");
    finish_page();
  }
}

my $baseURL = "compare.cgi?${options}";
if ($gene2) {
  $baseURL .= "&locus2=$gene2->{locusTag}";
} else {
  my $seqDesc2E = encode_entities($seqDesc2);
  $baseURL .= "&seqDesc2=$seqDesc2E&seq2=$seq2";
}

unless (hasHits($seq2)) {
  die "No homologs for second protein yet\n" if $tsv;
  print p("Sorry, homologs for the second protein have not been computed yet");
  print a({-href => "findHomologs.cgi?${options2}&compare1=" . uri_escape($options) }, "Compute homologs for the second protein");
  finish_page();
}
my $hits2 = getHits($seq2);
if (scalar(@$hits2) == 0) {
  exit(0) if $tsv;
  print p("Sorry, no homologs were found for the second protein");
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
print p("Loaded", commify(scalar(@$geneHits1)), "homologs for the first protein",
        "and", commify(scalar(@$geneHits2)), "homologs for the second protein."), "\n"
  unless $tsv;

if ($taxLevel) { # taxon distribution mode
  die "Cannot set tsv with taxLevel" if $tsv;
  my @levelsAllowed = taxLevels();
  shift @levelsAllowed; # don't allow the top level (domain or order)
  if (getOrder() eq "") {
    # no genus or species
    pop @levelsAllowed;
    pop @levelsAllowed;
  }
  my %levelsAllowed = map { $_ => 1 } @levelsAllowed;
  die "Invalid taxLevel $taxLevel\n" unless exists $levelsAllowed{$taxLevel};

  # Which taxon levels to include in the analysis
  my @levels = taxLevels(); # all the ones in the db
  # Remove levels after $taxLevel
  my ($iLevelAt) = grep $levels[$_] eq $taxLevel, (0..(scalar(@levels)-1));
  splice @levels, $iLevelAt + 1;
  # Remove order if subdb (the order is always the same)
  shift @levels if getOrder() ne "" && $levels[0] eq "order";
  my @levelsShow = map capitalize($_), @levels;

  die "Invalid taxMode" unless $taxMode eq "both" || $taxMode eq "close";
  my $goodOnly = $cgi->param('Good') || 0;
  my $all = $cgi->param('all') || "freq";
  die "Invalid all" unless $all eq "all" || $all eq "freq";

  # All hits by genome
  my %byg1 = ();
  foreach my $hit (@$geneHits1) {
    push @{ $byg1{ $hit->{gid} } }, $hit if !$goodOnly || $hit->{bits} >= $max1 * 0.3;
  }
  my %byg2 = ();
  foreach my $hit (@$geneHits2) {
    push @{ $byg2{ $hit->{gid} } }, $hit if !$goodOnly || $hit->{bits} >= $max2 * 0.3;
  }

  # genome to list of relevant hits
  my %gid1 = ();
  my %gid2 = ();
  if ($taxMode eq "both") {
    foreach my $gid (keys %byg1) {
      if (exists $byg2{$gid}) {
        push @{ $gid1{$gid} }, @{ $byg1{$gid} };
        push @{ $gid2{$gid} }, @{ $byg2{$gid} };
      }
    }
  } else { # keep genome only if there's a close pair
    foreach my $gid (keys %byg1) {
      next unless exists $byg2{$gid};
      # Keep the two seen lists separate to ensure that byg2 is non-empty whenever
      # byg1 is added to. This was causing a bug before I added the check that
      # hit1 and hit2 have different locusTag; with that check, I'm not sure
      # if the two seen lists could actually overlap, but maybe it is possible.
      my %seen1 = (); # locusTag => 1 if already listed in byg1
      my %seen2 = (); # locusTag => 1 if already listed in byg2
      my $list1 = $byg1{$gid};
      my $list2 = $byg2{$gid};
      my $keep = 0;
      foreach my $hit1 (@$list1) {
        foreach my $hit2 (@$list2) {
          if ($hit1->{scaffoldId} eq $hit2->{scaffoldId}
              && $hit1->{locusTag} ne $hit2->{locusTag}
              && $hit1->{strand} eq $hit2->{strand}
              && (abs($hit1->{begin} - $hit2->{end}) <= $closeNt
                  || abs($hit1->{end} - $hit2->{begin}) <= $closeNt)) {
            push @{ $gid1{$gid} }, $hit1 unless exists $seen1{ $hit1->{locusTag} };
            $seen1{ $hit1->{locusTag} } = 1;
            push @{ $gid2{$gid} }, $hit2 unless exists $seen2{ $hit2->{locusTag} };
            $seen2{ $hit2->{locusTag} } = 1;
            $keep = 1;
            last;
          }
        }
        last if $keep;
      }
    }
  }

  my $hitsString = "homologs";
  $hitsString = a({-title => "At least 30% of maximum possible bit score" }, "good homologs")
    if $goodOnly;
  my $whatString = "for both genes";
  $whatString = "nearby" if $taxMode eq "close";
  my %taxonPlural = qw{phylum phyla class classes order orders family families genus genera species species};
  print p("Showing which $taxonPlural{$taxLevel} have", $hitsString, $whatString . ".",
          "Genomes with",
          ($goodOnly ? "good homologs for" : ""),
          ($taxMode eq "both" ? "both genes" : "the genes within $closeKb kb and on the same strand").":",
          scalar(keys %gid1)."."), "\n";
  my $hidden2;
  if ($gene2){
    $hidden2 = qq[<INPUT type="hidden" name="locus2" value="$gene2->{locusTag}">];
  } else {
    my $seqDesc2E = encode_entities($seqDesc2);
    $hidden2 = qq[<INPUT type="hidden" name="seqDesc2" value="$seqDesc2E">]
      . qq[<INPUT type="hidden" name="seq2" value="$seq2">];
  }
  my %modeLabels = ('close' => 'nearby' , 'both' => 'both present');
  print
    start_form( -name => 'input', -method => 'GET', -action => 'compare.cgi'),
    $hidden,
    $hidden2,
    p("Genes are:",
      popup_menu(-name => 'taxMode', -values => [qw(close both)],
                 -default => $taxMode, -labels => \%modeLabels),
      "&nbsp;",
      "Level:",
      popup_menu(-name => 'taxLevel', -values => \@levelsAllowed,
                 -default => $taxLevel),
      "&nbsp;",
      a({-title => "a good homolog has at least 30% of the maximum possible bit score"},
        checkbox(-name => 'Good', -checked => $goodOnly),
        "homologs only?"),
      "&nbsp;",
      "Show:",
      popup_menu(-name => 'all', -values => [ "all", "freq" ],
                 -labels => { "all" => "all taxa", "freq" => "taxa with homologs" },
                 -default => $all),
      "&nbsp;",
      submit('Change')),
  end_form,
  "\n";

  if (scalar(keys %gid1) > 0) {
    my $genomes = getDbHandle()->selectall_hashref("SELECT * from Genome", "gid");
    # taxString is each taxon in this list of levels, joined by ";;;"
    my %taxStringN = ();
    my %taxStringGid = (); # taxString to list of gid
    foreach my $genome (values %$genomes) {
      $genome->{taxString} = join(";;;", map $genome->{"gtdb".$_}, @levelsShow);
      $taxStringN{ $genome->{taxString} }++;
      push @{ $taxStringGid{ $genome->{taxString} } }, $genome->{gid};
    }
    my %taxHits = (); # taxString to list of relevant hits
    while (my ($gid, $hits) = each %gid1) {
      my $genome = $genomes->{ $gid } || die $gid;
      push @{ $taxHits{$genome->{taxString}} }, @$hits;
      push @{ $taxHits{$genome->{taxString}} }, @{ $gid2{$gid} };
    }
    # Each row includes taxString, taxLevels, nHitGenomes, nGenomes
    my @rows = ();
    if ($all eq "freq") {
      while (my ($taxString, $tHits) = each %taxHits) {
        my $row = { 'taxString' => $taxString };
        $row->{nGenomes} = $taxStringN{$taxString} || die $taxString;
        my %gidThis = (); # genome id to 1 if has a hit
        foreach my $hit (@{ $taxHits{$taxString} }) {
          $gidThis{ $hit->{gid} } = 1;
        }
        $row->{nHitGenomes} = scalar(keys %gidThis);
        push @rows, $row;
      }
      @rows = sort { $b->{nHitGenomes} <=> $a->{nHitGenomes}
                       || $a->{taxString} cmp $b->{taxString} } @rows;
    } else {
      # all taxa mode
      my $taxa = getTaxa();
      my @levels = map lc, @levelsShow;
      foreach my $tax (values %{ $taxa->{$taxLevel} }) {
        my $levels = taxToParts($tax, $taxa);
        my $taxString = join(";;;", map $levels->{$_}, @levels);
        my $nHits = 0;
        my $nHitGenomes = 0;
        if (exists $taxHits{$taxString}) {
          my %gids = map { $_->{gid} => 1 } @{ $taxHits{$taxString} };
          $nHitGenomes = scalar(keys %gids);
        }
        die $taxString unless defined $taxStringN{$taxString};
        push @rows, { 'taxString' => $taxString,
                      'nGenomes' => $taxStringN{$taxString},
                      'nHitGenomes' => $nHitGenomes };
      }
      @rows = sort { $a->{taxString} cmp $b->{taxString} } @rows;
    }

    my @header = @levelsShow;
    $header[0] = "&nbsp;" if getOrder() eq "";
    push @header, "#Genomes";
    $hitsString = $goodOnly ? "good homologs" : "homologs";
    if ($taxMode eq "close") {
      push @header, a({-title => "#Genomes with $hitsString nearby"}, "#Nearby");
    } else {
      push @header, a({-title => "#Genomes with $hitsString for both genes"}, "#Both");
    }
    my $homologString = $goodOnly ? "good homologs" : "homologs";
    push @header,
      map(a({-title => "Bit score ratio for the best homolog of protein $_, among "
             . ($taxMode eq "close" ? "nearby pairs of $homologString"
                : "genomes with $homologString for both") },
         "Max ratio$_"), 1..2);
    print qq[<TABLE cellpadding=1 cellspacing=1>], "\n";
    print Tr(map th($_), @header), "\n";
    my $iRow = 0;
    foreach my $row (@rows) {
      my @lineage = split /;;;/, $row->{taxString};
      my @out = @lineage;
      for (my $i = 0; $i < scalar(@levels); $i++) {
        if ($levels[$i] eq "domain") {
          $out[$i] = domainHtml($out[$i]);
        } else {
          $out[$i] = a({-href => addOrderToURL("taxon.cgi?level=".$levels[$i]
                                               ."&taxon=".uri_escape($out[$i])),
                        -style => "text-decoration:none;"},
                       encode_entities($out[$i]));
        }
      }
      my $showNHit;
      my $showMax1 = "&nbsp;";
      my $showMax2 = "&nbsp;";
      if (exists $taxHits{ $row->{taxString} }) {
        my @tHits = @{ $taxHits{ $row->{taxString} } };
        my @maxBits; # 1 => max bits1, 2 => max bits2
        foreach my $i (1,2) {
          my $gidI = $i == 1 ? \%gid1 : \%gid2;
          $maxBits[$i] = 0;
          foreach my $gid (@{ $taxStringGid{$row->{taxString}} }) {
            foreach my $hit (@{ $gidI->{$gid} }) {
              $maxBits[$i] = $hit->{bits} if $hit->{bits} > $maxBits[$i];
            }
          }
        }
        $showMax1 = sprintf("%.2f", $maxBits[1] / $max1);
        $showMax2 = sprintf("%.2f", $maxBits[2] / $max2);
        my $truncate = 0;
        my $maxShow = 250;
        if (@tHits > $maxShow) {
          $truncate = 1;
          splice @tHits, $maxShow;
        }
        my $URL = addOrderToURL("genes.cgi?" . join("&", map "g=" . $_->{locusTag}, @tHits));
        my $truncateString = $truncate ? " top $maxShow" : "";
        my $modeString = $taxMode eq "close" ? "close-by" : "co-occurring";
        $showNHit = a({ -href => $URL,
                        -style => "text-decoration: none;",
                        -title => "see$truncateString $modeString $hitsString in $row->{nHitGenomes} "
                        . encode_entities($lineage[-1]) },
                      commify($row->{nHitGenomes}));
      } else {
        $showNHit = "&nbsp;";
      }
      
      my $bgColor = $iRow++ % 2 == 0 ? "lightgrey" : "white";
      print Tr({-style => "background-color: $bgColor;"},
               td(\@out),
               td({-style => "text-align: right;"},
                  [ commify($row->{nGenomes}), $showNHit, $showMax1, $showMax2 ])) . "\n";
    }
    print "</TABLE>\n";
  } # end has genomes to show
  print p("Or see", a({-href => $baseURL}, "presence/absence plots"));
  finish_page();
} # end taxon distribution mode

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

# Put all the statistics on the left and the svg (if there is one) on the right, using a table
# with 1 row and 2 columns
print
  qq{<TABLE cellpadding=2 cellspacing=2><TR><TD valign="top">},
  h4("Co-occurence: all homologs"),
  p("Found homologs in", commify($n1), "and", commify($n2), "genomes, respectively.",
    "$nInBoth genomes contain homologs of both genes (versus",
    commify(int($nInBothExpected+0.5)), "expected)."),
  p("Considering only the best hit in each genome,",
        "$nClose hits are nearby (within $closeKb kb) and on the same strand, and $nSame are to the same gene."),
  h4("Co-occurence: good homologs"),
  p("Found good homologs (above 30% of maximum bit score) in",
    commify($n1good), "and", commify($n2good), "genomes, respectively.",
    "$nGoodBoth genomes contain good homologs of both genes (versus",
    commify(int($nGoodBothE+0.5)), "expected).");
print
  p("Among the best hits in those $nGoodBoth genomes,",
    "$nGoodClose are nearby (within $closeKb kb) and on the same strand, and $nGoodSame are to the same gene.")
  if $nGoodBoth > 0;

# Find the most significant threshold (if any)
# for the enrichment of the two being in the same genome,
# and plot the similarity in gene presence
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
  print h4("Co-occurrence: optimal threshold");
  my ($thresh1,$thresh2);
  if ($correctedLog10P < -3) {
    my %gidsThresh1 = map { $_ => 1 } @gids1[0..$bestI];
    my %gidsThresh2 = map { $_ => 1 } @gids2[0..$bestI];
    my @gidsThresh = grep exists $gidsThresh2{$_}, keys %gidsThresh1;
    my $nCloseThresh = scalar(grep $gidScores{$_}{close}, @gidsThresh);
    $thresh1 = $byg1{ $gids1[$bestI] }{bits};
    my $f1 = int(0.5 + 100 * $thresh1/$max1);
    $thresh2 = $byg2{ $gids2[$bestI] }{bits};
    my $f2 = int(0.5 + 100 * $thresh2/$max2);
    my $fClose = int(0.5 + 100 * $nCloseThresh/scalar(@gidsThresh));
    print p("The strongest statistical signal of co-occurence is for the top",
            commify($bestI+1), "homologs of each gene, which co-occur in",
            commify(scalar(@gidsThresh)), "genomes,",
            "P = 10<sup>" . sprintf("%.1f", $correctedLog10P) . "</sup>",
            "(Fisher's exact test, 1-sided, with Bonferonni correction).",
            "The corresponding bit score thresholds are $thresh1 ($f1% of max)",
            " and $thresh2 ($f2% of max), respectively.",
           "$nCloseThresh of these co-occurring homologs ($fClose%) are nearby (within $closeKb kb and on the same strand).");
    print p({-style => "font-size: 90%;"}, "Warning: the P value is based on the assumption that the two genes appear in genomes independently. If both genes are conserved within a group of related genera, then they will have significant co-occurrence, even if there is no functional relationship between them. Nevertheless, the P value is useful for selecting a threshold.");
  } else {
    print p("No significant co-occurence");
  }

  # Show a scatterplot of ratio1 vs. ratio2 (or below 0 if no hit)
  my $below = 0.1; # space at the bottom of non-hits, so you can see them better
  my $plot = svgPlot->new(xRange => [-$below,1], yRange => [-$below,1],
                          width => 525, height => 500,
                          margin => [3,3,3,4]);
  my @ticks = (0,0.2,0.4,0.6,0.8,1.0);
  print qq{</TD><TD valign="top">};
  print $plot->start();
  my $x0 = $plot->convertX(0);
  my $y0 = $plot->convertY(0);
  my $drawWidth = $plot->{drawRight} - $plot->{drawLeft};
  my $drawHeight = $plot->{drawBottom} - $plot->{drawTop};
  my $below0Width = $x0 - $plot->{drawLeft};
  my $below0Height = $plot->{drawBottom} - $y0;
  my $greyTitle = "genomes with homologs for only one of the genes are shown in the grey region";
  print join("\n",
             qq{<rect x="$plot->{drawLeft}" y="$y0" width="$drawWidth" height="$below0Height" fill="lightgrey" stroke="none">},
             qq{<title>$greyTitle</title></rect>},
             qq{<rect x="$plot->{drawLeft}" y="$plot->{drawTop}" width="$below0Width" height="$drawHeight" fill="lightgrey" stroke="none">},
             qq{<title>$greyTitle</title></rect>});
  {
    my $x = $plot->convertX(0.3);
    my $y = $plot->convertY(0.3);
    my $params = qq{ stroke="orange" stroke-dasharray="4" };
    my $title = qq{The 'good' threshold (both score ratios >= 0.3)};
    print qq{<line x1="$x" y1="$y" x2="$x" y2="$plot->{drawTop}" $params ><title>$title</title></line>};
    print qq{<line x1="$x" y1="$y" x2="$plot->{drawRight}" y2="$y" $params ><title>$title</title></line>};
    print qq{<text x="$plot->{drawRight}" y="$y" dominant-baseline="bottom" fill="orange"><title>$title</title>good</text>};
  }
  if (defined $thresh1) {
    my $x = $plot->convertX($thresh1/$max1);
    my $y = $plot->convertY($thresh2/$max2);
    my $params = qq{ stroke="black" stroke-dasharray="4" };
    my $title = qq{The 'optimal' threshold, from maximizing -log(P)};
    print qq{<line x1="$x" y1="$y" x2="$x" y2="$plot->{drawTop}" $params ><title>$title</title></line>};
    print qq{<line x1="$x" y1="$y" x2="$plot->{drawRight}" y2="$y" $params ><title>$title</title></line>};
    print qq{<text x="$plot->{drawRight}" y="$y" dominant-baseline="bottom" fill="black"><title>$title</title>optimal</text>};
  }
  my $max1Show = int(0.5 + $max1);
  my $max2Show = int(0.5 + $max2);
  print
    $plot->axes(),
    $plot->axisTicks("x", \@ticks),
    $plot->axisTicks("y", \@ticks),
    $plot->marginText("Score ratios in each genome", "top",
                      title => "For each genome with a homolog of either gene, the score ratios (bits/max) for the best hits",
                      style => "font-size: larger; font-weight: bold;"),
    $plot->marginText("Score ratio for protein 1", "bottom",
                      title => "Score (bits) of best homolog of protein 1 / max ($max1Show bits)"),
    $plot->marginText("Score ratio for protein 2", "left",
                      title => "Score (bits) of best homolog of protein 2 / max ($max2Show bits)");

  # legend for color-coding, at top
  my $top = $plot->{drawTop} - $plot->{lineSize} / 2;
  my $labelLeft = $plot->convertX(0);
  my $labelRight = $plot->convertX(0.8);
  my $labelMid = ($labelLeft+$labelRight)/2;
  my @labelX = ($labelLeft, $labelMid, $labelRight);
  my @col = qw{green blue black};
  my @fill = ("green", "blue", "#EEEEEE");
  my @symbol = ("o", "+", "o");
  my @labelTitles = ("Within $closeKb kb and on the same strand",
                     "Both proteins have the same best hit",
                     "The best hits are not near each other");
  my @labels = ("close by", "same", "other");
  foreach my $class (0..2) {
    print $plot->pointAbsolute($labelX[$class] - 6, $top - $plot->{lineSize}/3,
                               'size' => 2.5,
                               'color' => $col[$class], 'fill' => $fill[$class], 'symbol' => $symbol[$class]);
    print qq{<text x="$labelX[$class]" y="$top" text-anchor="left" alignment-baseline="center" fill="$col[$class]">},
      qq{<title>$labelTitles[$class]</title>$labels[$class]</text>}, "\n";
  }

  # points, with random moderately-negative value instead of 0
  srand(01312023);
  foreach my $v (reverse @gidValues) {
    my $gid = $v->{gid};
    my $genome = $genomes->{$gid};
    my $x = $v->{ratio1};
    $x = -$below/10 - rand(0.8*$below) if $x == 0;
    my $y = $v->{ratio2};
    $y = -$below/10 - rand(0.8*$below) if $y == 0;
    my @locusTags = ();
    push @locusTags, $byg1{$gid}{locusTag} if exists $byg1{$gid};
    push @locusTags, $byg2{$gid}{locusTag} if exists $byg2{$gid};
    my $lineage = join(" ",
                       $genome->{gtdbPhylum}, $genome->{gtdbClass}, $genome->{gtdbOrder},
                       $genome->{gtdbFamily});
    my $title = join(" and ", @locusTags);
    $title .= " (close by)" if $v->{close};
    $title = $locusTags[0] . " (similar to both)" if $v->{same};
    my $URL = addOrderToURL("genes.cgi?" .  join("&", map "g=$_", @locusTags));
    $URL = addOrderToURL("gene.cgi?locus=$locusTags[0]") if $v->{same};
    my $class = $v->{close} ? 0 : ($v->{same} ? 1 : 2);
    print $plot->point($x, $y,
                       size => 2.5,
                       color => $col[$class],
                       fill => $fill[$class],
                       symbol => $symbol[$class],
                       URL => $URL,
                       title => "$title in " . $genome->{gtdbSpecies}
                       . " ($lineage)"), "\n";
  }
  print $plot->end();
}
print "</TD></TR></TABLE>";

my $downloadURL = "$baseURL&format=tsv";
my $defaultTaxLevel = getOrder() eq "" ? "phylum" : "family";
print
  p("Or see which taxa have",
    a({-href => "$baseURL&taxLevel=${defaultTaxLevel}&taxMode=both"},
      "both genes"),
    "or have",
    a({-href => "$baseURL&taxLevel=${defaultTaxLevel}&taxMode=close"},
      "the genes nearby,"),
    "or",
    a({ -href => $downloadURL }, "download a table"),
    "(tab-delimited)",
    "of the best homolog(s) in each genome");
finish_page();
