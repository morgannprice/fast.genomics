#!/usr/bin/perl -w
# CGI to find similarities for a 16S sequence
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

# Optional CGI arguments:
# query
#
# If no sequence is specified, shows an input form. Unlike other pages, locus tags are not supported as inputs.

my $cgi = CGI->new;
my $minAln = 100; # minimum alignment length
my $minIdentity = 0.85; # minimum %identity
my $maxHits = 200;

my $query = param('query') || "";
my ($seqDesc, $seq);

my $inputError;
if ($query ne "" && $query !~ m/^\s+$/) {
   $query =~ s/^\s+//;
   $query =~ s/\s+$//;
   my @lines = split /\n/, $query;
   if ($query =~ m/^>/ && @lines > 1) {
     $seqDesc = shift @lines;
     $seqDesc =~ s/^>//;
     $seq = join("", @lines);
   } else {
     $seqDesc = "";
     $seq = $query;
   }
   $seq =~ s/\s//g;
   $seq = uc($seq);
   $seq =~ s/U/T/g;
   if ($seq eq "") {
     $inputError = "No sequence";
   } elsif ($seq !~ m/^[A-Z]+$/) {
     $inputError = "Input sequence must contain characters only";
   } elsif (length($seq) < $minAln) {
     $inputError = "Input sequence must be at least $minAln characters";
   } else {
     my @nt = $seq =~ m/[ACGTN]/g;
     my $fACGTN = scalar(@nt) / length($seq);
     if ($fACGTN < 0.95) {
       $inputError = "Input sequence must be a nucleotide sequence";
     } else {
       my @n = $seq =~ m/N/g;
       my $fN = scalar(@n) / length($seq);
       if ($fN > 0.5) {
         $inputError = "Input sequence is majority Ns!";
       }
     }
   }
 }

$seq = "" if $inputError;

start_page('title' => "Search for similar 16S sequences");
autoflush STDOUT 1; # show preliminary results

print p("Invalid input:", $inputError) if $inputError;

if (! $seq) {
  print
    p(start_form(-name => 'inputForm', -method => 'GET', -action => '16Ssim.cgi'),
      "Enter a 16S sequence, by itself or in fasta format",
      br(),
      textarea( -name  => 'query', -value => '', -cols  => 90, -rows  => 8, -override => 1 ),
      br(),
      submit('16S search'),
      end_form);
    finish_page();
} # end no query

if ($seqDesc eq "") {
  $seqDesc = substr($seq, 0, 10) . "... (" . length($seq) . " characters)";
}

my $usearch = "../bin/usearch";
die "No such executable: $usearch\n" unless -x $usearch;
my $udb = "../data/filtered16S.udb";
die "No such file: $udb\n" unless -e $udb;

print p("Searching for similar 16S sequences for", encode_entities($seqDesc));
print showSequence($seqDesc, $seq);
print p(sprintf("Minimum identity: %d%%, minimum alignment: %d nucleotides, query length: %d",
                int($minIdentity*100+0.5), $minAln, length($seq)));
print "\n";

my $tmpDir = $ENV{TMPDIR} || "/tmp";
my $tmpPre = "$tmpDir/16Ssim.$$";
my $tmpFna = "$tmpPre.fna";
my $tmpHits = "$tmpPre.hits";
open(my $fh, ">", $tmpFna) || die "Cannot write to $tmpFna";
print $fh ">query\n$seq\n";
close($fh) || die "Error writing to $tmpFna";

my $cmd = "$usearch -usearch_local $tmpFna -db $udb -id $minIdentity -mincols $minAln -strand both"
  . " -maxaccepts $maxHits -maxrejects $maxHits -blast6out $tmpHits -quiet --threads 1";
system($cmd) == 0 || die "usearch failed\n" . join(" ", $cmd) . "\n -- $!";

my @hits = ();
open($fh, "<", $tmpHits) || die "Cannot read $tmpHits";
while(my $line = <$fh>) {
  chomp $line;
  my @F = split /\t/, $line;
  push @hits, \@F;
}
close($fh) || die "Error reading $tmpHits";

unlink($tmpFna);
unlink($tmpHits);
splice @hits, $maxHits if @hits > $maxHits;
my $found = "Found " . scalar(@hits) . " hits";
$found .=  " (the maximum)" if @hits == $maxHits;
$found .= ", or try " . a({-href => "16Ssim.cgi"}, "another query");
print $found, "\n";

if (@hits > 0) {
  print
    qq{<TABLE cellpadding=2 cellspacing=2>},
    Tr(th({-style => "text-align: left"}, ['%id.', 'len.', 'Taxonomy', 'Strain', 'Locus'])),
    "\n";
  my $iRow = 0;
  foreach my $row (@hits) {
    my (undef, $subject, $identity, $alen, $mm, $gap, $qbeg, $qend, $sbeg, $send, $bits) = @$row;
    my ($gid, $locusTag) = split /:/, $subject;
    die unless defined $locusTag;
    my $genome = getTopDbHandle()->selectrow_hashref("SELECT * from AllGenome WHERE gid = ?", {}, $gid);
    die "Unknown gid $gid" unless defined $genome;
    my $orderParam = "";
    $orderParam = $genome->{prefix} unless $genome->{inTop};
    my @out = ();
    push @out, td({-style => 'text-align: right;'}, $identity);
    push @out, td({-style => 'text-align: right;'},
                  a({-title => "$qbeg:$qend of query aligns to $sbeg:$send of $locusTag"}, $alen));
    my @taxParts;
    foreach my $level (qw{Phylum Class Order Family}) {
      my $taxon = $genome->{"gtdb".$level};
      my $URL = "taxon.cgi?taxon=${taxon}&level=" . lc($level);
      push @taxParts, a({-href => $URL, -style => 'text-decoration: none;'}, $taxon);
    }
    push @taxParts, a({-href => "taxon.cgi?level=genus&order=$orderParam&taxon=$genome->{gtdbGenus}",
                       -style => 'text-decoration: none;'},
                      i($genome->{gtdbSpecies}));
    push @out, td(join(", ", @taxParts));
    push @out, td(a({-style => 'text-decoration: none;', -href => "genome.cgi?gid=$gid&order=$orderParam"},
                  encode_entities($genome->{strain}) ));
    push @out, td(a({-style => 'text-decoration: none;', -href => "gene.cgi?locus=$locusTag&order=$orderParam"},
                    $locusTag));
    my $bgColor = $iRow++ % 2 == 0 ? "white" : "lightgrey";
    print Tr({-style => "background-color: $bgColor"}, @out), "\n";
  }
  print "</TABLE>";
}
finish_page();
