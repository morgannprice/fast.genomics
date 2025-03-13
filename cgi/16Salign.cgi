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

# Required CGI arguments:
# query -- fasta format header and sequence
# subject -- a locusTag that is in the All16S table

my $cgi = CGI->new;
my $query = param('query') || die "Must specify query";
my $subject = param('subject') || die "Must specify subject";
my $sObj = getTopDbHandle()->selectrow_hashref("SELECT * from All16S WHERE locusTag = ?",
                                               {}, $subject);
die "Unknown subject" unless $sObj;
my @lines = split /\n/, $query;
my $header = shift @lines;
die "No header" unless $header =~ m/^>/;
$header =~ s/^>//;
my $seq = join("", @lines);
$seq =~ s/\s//g;
die "Invalid sequence" unless $seq =~ m/^[A-Z]+$/;

start_page('title' => "16S Alignment to $subject");

my $sURL = "gene.cgi?locus=$subject";
my $sg = getTopDbHandle()->selectrow_hashref("SELECT * from AllGenome WHERE gid = ?",
                                             {}, $sObj->{gid});
$sURL .= "&order=" . $sg->{prefix} unless $sg->{inTop};
print h3("Aligning", encode_entities($header), "to",
         a({-href => $sURL}, $subject));
print p("Query length:", length($seq));
print showSequence($header, $seq);

my $tmpDir = $ENV{TMPDIR} || "/tmp";
my $tmpPre = "$tmpDir/16Salign.$$";
my $bl2seq = "../bin/blast/bl2seq";
die "No such executable: $bl2seq" unless -x $bl2seq;

my $tmpQ = "$tmpPre.query";
open(my $fhQ, ">", $tmpQ) || die "Cannot write to $tmpQ";
print $fhQ ">query\n${seq}\n";
close($fhQ) || die "Error writing to $tmpQ";

my $tmpS = "$tmpPre.subject";
open(my $fhS, ">", $tmpS) || die "Cannot write to $tmpS";
print $fhS ">$subject\n" . $sObj->{sequence} . "\n";
close($fhS) || die "Error writing to $tmpS";

my $tmpOut = "$tmpPre.bl2seq";
my $cmd = "bl2seq -p blastn -e 1e-5 -i $tmpQ -j $tmpS -o $tmpOut";
system($cmd) == 0 || die "bl2seq failed\n$cmd\n$!";

my @aln;
open(my $fh, "<", $tmpOut) || die "Cannot read $tmpOut";
@aln = <$fh>;
close($fh) || die "Error reading $tmpOut";

unlink($tmpQ);
unlink($tmpS);
unlink($tmpOut);

print "<PRE>", "\n";
my $writing = 0;
foreach my $line (@aln) {
  $writing = 1 if $line =~ m/^>/;
  $writing = 0 if $line =~ m/^Lambda/;
  print $line if $writing;
}
print "</PRE>\n";
finish_page();
