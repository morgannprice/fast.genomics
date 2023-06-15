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

# required CGI arguments:
#   g -- list of locus tags
# optional CGI arguments:
#   order (which subdb to use)
#   format=faa or blastp

my $cgi = CGI->new;
setOrder(param('order'));
my $format = param('format') || "";
my @g = $cgi->param('g');
splice @g, 250 if @g > 250;
my @genes = map { locusTagToGene($_) || die "Unknown locusTag $_" } @g;
my @proteinGenes = grep $_->{proteinId} ne "", @genes;

if ($format eq "faa") {
  print "Content-Type:text\n";
  print "Content-Disposition: attachment; filename=fast_genomics.faa\n\n";
  foreach my $g (@proteinGenes) {
    my ($seq) = getDbHandle()->selectrow_array("SELECT sequence FROM Protein WHERE proteinId = ?",
                                               {}, $g->{proteinId});
    die "No sequence for " . $g->{locusTag} unless $seq;
    my $genome = gidToGenome($g->{gid});
    die "Unknown gid $g->{gid}" unless $genome;
    print ">" . $g->{locusTag} . " " . $g->{proteinId} . " " . $g->{desc}
      . " [" . $genome->{gtdbSpecies} . " " . $genome->{strain} . "]"
      . "\n" . $seq . "\n";
  }
  exit(0);
}

# If they all in the same genome, show the genome info at the top,
# and give more details

my %gid = map { $_->{gid} => 1 } @genes;
my $constGenome = gidToGenome( (keys %gid)[0] )
  if scalar(keys %gid) == 1;

my $title = scalar(@g) . " genes";
$title .= " in $constGenome->{gtdbSpecies} $constGenome->{strain}"
  if $constGenome;

start_page('title' => $title);
autoflush STDOUT 1; # show preliminary results

if ($format eq "blastp" && @proteinGenes > 1) {
  if (@proteinGenes > 100) {
    print p("Only considering the first 100 proteins"), "\n";
    splice @proteinGenes, 100;
  }
  # Fetch all the protein sequences and run BLASTp on them
  my $tmpFaa = "/tmp/genes_blastp_$$.faa";
  open (my $fhFaa, ">", $tmpFaa) || die "Cannot write to $tmpFaa\n";
  foreach my $g (@proteinGenes) {
    $g->{seq} = getDbHandle()->selectrow_array("SELECT sequence FROM Protein WHERE proteinId = ?",
                                               {}, $g->{proteinId});
    die "No sequence for " . $g->{locusTag} unless $g->{seq};
    $g->{genome} = gidToGenome($g->{gid}) || die "Unknown gid $g->{gid}";
    print $fhFaa ">" . $g->{locusTag} . "\n" . $g->{seq} . "\n";
  }
  close($fhFaa) || die "Error writing to $tmpFaa";
  system("../bin/blast/formatdb -p T -i $tmpFaa") == 0 || die "formatdb failed on $tmpFaa: $!";
  print p("Comparing", scalar(@proteinGenes), " proteins with BLAST, or see",
          a({-href => addOrderToURL("genes.cgi?" . join("&", map "g=$_", @g))}, "list")),
    "\n";
  my $cmd = qq{../bin/blast/blastall -p blastp -i $tmpFaa -d $tmpFaa -m 8 -o $tmpFaa.hits -e 1e-3 -F "m S" -a 10};
  system($cmd) == 0 || die "blastall failed\n$cmd\n$!";
  open(my $fhHits, "<", "$tmpFaa.hits") || die "Cannot read $tmpFaa.hits";
  my @hits;
  while (my $line = <$fhHits>) {
    chomp $line;
    my @F = split /\t/, $line;
    push @hits, \@F;
  }
  close($fhHits) || die "Error reading $tmpFaa.hits";
  unlink("$tmpFaa");
  foreach my $suffix (qw{hits phr pin psq}) {
    unlink("$tmpFaa.$suffix");
  }
  # Remove self htis
  @hits  = grep $_->[0] ne $_->[1], @hits;
  if (@hits == 0) {
    print p("No significant similarities (E < 0.001) were found");
  } else {
    my %proteinGenes = map { $_->{locusTag} => $_ } @proteinGenes;
    print
      "<TABLE cellpadding=2 cellspacing=2>",
      Tr(th(["query", "query genome", "subject", "subject genome",
             a({-title => '%identity'}, '%id.'),
             a({-title => '%coverage of query'}, '%cov.')])),
      "\n";
    my $iRow = 0;
    foreach my $hit (@hits) {
      my ($query, $subject, $identity, $alen, $mm, $gap, $qbeg, $qend, $sbeg, $send, $evalue, $bits) = @$hit;
      die "Error in BLASTp output" unless defined $bits && $bits =~ m/ *\d+$/;
      my $qGene = $proteinGenes{$query} || die "Unknown query $query";
      my $sGene = $proteinGenes{$subject} || die "Unknown subject $subject";
      my $qOrg = $qGene->{genome}{gtdbSpecies} . " " . $qGene->{genome}{strain};
      my $sOrg = $sGene->{genome}{gtdbSpecies} . " " . $sGene->{genome}{strain};
      my $qLen = length($qGene->{seq});
      my $sLen = length($sGene->{seq});
      print Tr({-bgcolor => $iRow++ % 2 == 1 ? "white" : "lightgrey"},
               td({-valign => "top", -align => "left"},
                  a({ -href => addOrderToURL("gene.cgi?locus=$query"),
                      -title => encode_entities($qGene->{desc}),
                      -style => "text-decoration: none;" },
                    $query)),
               td({-valign => "top", -align => "left"},
                  a({ -href => addOrderToURL("genome.cgi?gid=".$qGene->{gid}),
                      -style => "text-decoration: none;" },
                    $qOrg)),
               td({-valign => "top", -align => "left"},
                  a({ -href => addOrderToURL("gene.cgi?locus=$subject"),
                      -title => encode_entities($sGene->{desc}),
                      -style => "text-decoration: none;" },
                    $subject)),
               td({-valign => "top", -align => "left"},
                  a({ -href => addOrderToURL("genome.cgi?gid=".$sGene->{gid}),
                      -style => "text-decoration: none;" },
                    $sOrg)),
               td({-valign => "top", -align => "right"},
                  a({ -title => "$bits bits, E = $evalue",
                      -style => "text-decoration: none;" },
                    $identity."%")),
               td({-valign => "top", -align => "right"},
                  a({ -title => "$qbeg:$qend/$qLen of $query aligns to $sbeg:$send/$sLen of $subject",
                      -href => addOrderToURL("alignPair.cgi?locus=$query&locus2=$subject"),
                      -style => "text-decoration: none;" },
                    int(0.5 + 100 * ($qend-$qbeg+1)/$qLen)."%"))
              ), "\n";
    }
    print "</TABLE>\n";
  }
  finish_page();
}

# else
if (@proteinGenes > 1) {
  my $gSpec = join("&", map "g=$_", @g);
  print
    p("Tools:",
      a({ -href => addOrderToURL("genes.cgi?format=blastp&$gSpec") },
        "compare protein sequences"),
      "or",
      a({ -href => addOrderToURL("genes.cgi?format=faa&$gSpec"),
          -title => "fasta format" },
        "download protein sequences"));
}


if ($constGenome) {
  my @lineage = ();
  foreach my $level (qw{phylum class order family}) {
    my $tax = $constGenome->{"gtdb" . capitalize($level)};
    my $URL = "taxon.cgi?level=$level&taxon=$tax";
    if (getOrder() ne "" && ($level eq "phylum" || $level eq "class")) {
      push @lineage, a({-href => $URL,
                        -title => "See $level $tax in the main database",
                        -style => "text-decoration: none; color: black;" },
                       $tax);
    } else {
      push @lineage, a({-href => addOrderToURL($URL),
                        -title => $level,
                        -style => "text-decoration: none;"},
                       $tax);
    }
  }
  print p("Lineage:",
          domainHtml($constGenome),
          @lineage,
          a({-href => addOrderToURL("genome.cgi?gid=$constGenome->{gid}"),
             -style => "text-decoration:none;"},
            i($constGenome->{gtdbSpecies})),
          $constGenome->{strain});
  print "<TABLE cellpadding=1 cellspacing=1>\n";
  my @header = qw{Locus Protein Description Begin End Strand};
  print Tr(th([ @header ])), "\n";
  my $iRow = 0;
  foreach my $gene (@genes) {
    my $bgColor = $iRow % 2 == 0 ? "lightgrey" : "white";
    print Tr({-style => "background-color: $bgColor;"},
             td([ a({-href => addOrderToURL("gene.cgi?locus=$gene->{locusTag}") }, $gene->{locusTag}),
                  small($gene->{proteinId}),
                  encode_entities($gene->{desc}),
                  commify($gene->{begin}),
                  commify($gene->{end}) ])
             . td({ -align => "center"}, $gene->{strand})), "\n";
    $iRow++;
  }
  print "</TABLE>\n";
} else {
  foreach my $gene (@genes) {
    my $genome = gidToGenome($gene->{gid}) || die $gene->{gid};
    my @lineage = ();
    foreach my $level (qw{phylum class order family}) {
      my $field = "gtdb" . capitalize($level);
      my $tax = $genome->{$field};
      my $URL = "taxon.cgi?level=$level&taxon=".uri_escape($tax);
      if (getOrder() ne "" && ($level eq "phylum" || $level eq "class")) {
        # Cannot link to this taxon in the subdb.
        push @lineage, a( {-href => $URL,
                           -title => "See $level $tax in the main database",
                           -style => "text-decoration: none; color: black;" },
                          $tax);
      } else {
        push @lineage, a( {-href => addOrderToURL($URL),
                           -title => $level,
                           -style => "text-decoration: none;" },
                          $tax);
      }
    }
    print p(domainHtml($genome),
            ":",
            join(" : ", @lineage),
            ":",
            a({-href => addOrderToURL("genome.cgi?gid=$genome->{gid}"),
               -style => "text-decoration:none;"},
              i($genome->{gtdbSpecies})),
            $genome->{strain},
            br(),
            a({-href => addOrderToURL("gene.cgi?locus=$gene->{locusTag}")}, $gene->{locusTag}),
            small($gene->{proteinId}),
            $gene->{desc});
  }
}

finish_page;
