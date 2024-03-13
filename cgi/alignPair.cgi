#!/usr/bin/perl -w
use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use FindBin qw{$RealBin};
use HTML::Entities;
use URI::Escape;
use lib "../lib";
use neighbor;
use svgPlot;
# neighborWeb.pm relies on various PaperBLAST libraries
use lib "../../PaperBLAST/lib";
use neighborWeb;

# CGI arguments:
# locus (a locus tag in the database), or, seqDesc and seq
# locus2, or, seq2 and seqDesc2
#   (and optionally the order)

my $cgi = CGI->new;
setOrder(param('order'));
my ($gene, $seqDesc, $seq) = getGeneSeqDesc($cgi);
my $options = geneSeqDescSeqOptions($gene,$seqDesc,$seq); # for 1st gene

my $locus2 = $cgi->param('locus2');
my $seq2 = $cgi->param('seq2');
my $seqDesc2 = $cgi->param('seqDesc2');
my $query2 = $cgi->param('query2') || "";
my $gene2;

if (defined $locus2 && $locus2 ne "") {
  $gene2 = locusTagToGene($locus2) || die "Unknown locus tag in locus2 parameter";
  $seqDesc2 = $gene2->{locusTag} . " " . $gene2->{desc};
  ($seq2) = getDbHandle()->selectrow_array("SELECT sequence FROM Protein WHERE proteinId = ?",
                                           { Slice => {} }, $gene2->{proteinId})
    if $gene2->{proteinId} ne "";
} elsif ($seq2) {
  $seq2 =~ m/^[A-Z]+$/ || die "Invalid sequence in seq2";
}
die "No second protein specified" unless defined $seq2;
my $options2 = geneSeqDescSeqOptions($gene2,$seqDesc2,$seq2);

start_page(title => "Compare two protein sequences");

unless (defined $seq && $seq ne "") {
  print p("Sorry, the first gene is not protein-coding");
  finish_page();
}
unless (defined $seq2 && $seq2 ne "") {
  print p("Sorry, the second gene is not protein-coding");
  finish_page();
}

print h3("First protein");
if ($gene) {
  my $genome = gidToGenome($gene->{gid});
  print p(a({-href => "gene.cgi?$options"}, $gene->{locusTag}) . ":",
          encode_entities($gene->{desc}),
          br(),
          "from", i($genome->{gtdbSpecies}), encode_entities($genome->{strain}));
} else {
  print p(a({-href => "seq.cgi?$options"}, encode_entities($seqDesc)));
}

print h3("Second protein");
if ($gene2) {
  my $genome = gidToGenome($gene2->{gid});
  print p(a({-href => "gene.cgi?$options2"}, $gene2->{locusTag}) . ":",
          encode_entities($gene2->{desc}),
          br(),
          "from", i($genome->{gtdbSpecies}), encode_entities($genome->{strain}));
} else {
  print p(a({-href => "seq.cgi?$options2"}, encode_entities($seqDesc2)));
}

print h3("Blast comparison");
# Compare the two proteins by running bl2seq
my $bl2seq = "$RealBin/../bin/blast/bl2seq";
die "No such executable: $bl2seq" unless -x $bl2seq;
# Sanitize identifiers for bl2seq
$seqDesc =~ s/[^a-zA-Z0-9._'" -]/./g;
$seqDesc2 =~ s/[^a-zA-Z0-9._'" -]/./g;
my $tmpPre = "/tmp/alignPair.$$";
my $tmp1 = "$tmpPre.1";
my $tmp2 = "$tmpPre.2";
open(my $fh1, ">", $tmp1) || die "Cannot write to $tmp1";
print $fh1 ">", $seqDesc, "\n", $seq, "\n";
close ($fh1) || die "Error writing to $tmp1";
open(my $fh2, ">", $tmp2) || die "Cannot write to $tmp2";
print $fh2 ">", $seqDesc2, "\n", $seq2, "\n";
close ($fh2) || die "Error writing to $tmp2";

print "\n<pre>\n";
my $cmd = "$bl2seq -p blastp -i $tmp1 -j $tmp2 -e 1e-5 '-F m S'";
system($cmd) == 0
  || die "bl2seq failed -- $!\nCommand; $cmd\n";
unlink($tmp1);
unlink($tmp2);
print "</pre>\n";
print p({-style => "font-size: 90%;"},
        getOrder() eq "" ? "BLAST bit scores are not exactly the same as the bit scores from mmseqs2, which are shown on other pages of this site." : "",
        "E-values on this page will be much lower than on other pages because it is comparing to just one protein instead of a large database.");
finish_page();
