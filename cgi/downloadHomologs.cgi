#!/usr/bin/perl -w
use strict;
#!/usr/bin/perl -w
use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use lib "../lib";
# neighborWeb.pm relies on various PaperBLAST libraries
use lib "../../PaperBLAST/lib";
use neighborWeb;

my $cgi = CGI->new;
my ($gene, $seqDesc, $seq) = getGeneSeqDesc($cgi);
die "Not protein sequence\n" unless defined $seq;

unless (hasMMSeqsHits($seq)) {
  print redirect(-url => "findHomologs.cgi?" . geneSeqDescSeqOptions($gene,$seqDesc,$seq));
  exit(0);
}

my $hits = getMMSeqsHits($seq);
if (scalar(@$hits) == 0) {
  start_page('title' => 'Error');
  print p("Sorry, no homologs were found for this sequence");
  finish_page;
}
my $geneHits = hitsToGenes($hits);
my $genomes = getDbHandle()->selectall_hashref("SELECT * from Genome", "gid");

my $filename = defined $gene ? "homologs_$gene->{locusTag}.tsv"
  : "homologs_" . length($seq) . ".tsv";
print "Content-Type:text/tab-separated-values\n";
print "Content-Disposition: attachment; filename=homologs.tsv\n\n";

my @fields = qw{locusTag proteinId assemblyId
                scaffoldId geneBegin geneEnd strand
                gtdbDomain gtdbPhylum gtdbClass gtdbOrder gtdbFamily gtdbGenus gtdbSpecies strain
                identity alnLength nGapOpens
                qBegin qEnd sBegin sEnd eValue bits};
print join("\t", @fields)."\n";
foreach my $hit (@$geneHits) {
  $hit->{assemblyId} = $hit->{gid};
  $hit->{geneBegin} = $hit->{begin};
  $hit->{geneEnd} = $hit->{end};
  my $genome = $genomes->{$hit->{gid}} || die $hit->{gid};
  foreach my $field (qw{gtdbDomain gtdbPhylum gtdbClass gtdbOrder gtdbFamily gtdbGenus gtdbSpecies strain}) {
    $hit->{$field} = $genome->{$field};
  }
  print join("\t", map $hit->{$_}, @fields) . "\n";
}
