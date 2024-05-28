#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use DBI;
use FindBin qw{$RealBin};
use lib "$RealBin/../../PaperBLAST/lib";
use pbutils qw{ReadFastaEntry};

my @outFields = qw{query subject identity alnLength mismatch nGaps qBegin qEnd sBegin sEnd
                   sSequence};

my $usage =<<END
searchKmerSeekDb.pl -query in.fasta -db db.fasta -kmerdb kmer.db -out out.tsv

Finds similar sequences via matching first or last kmers and reports the top one.
The input file should have a single sequence. The output is tab-delimited with fields
  @outFields
where "q" and "s" are short for query and subject, and identity ranges from 0 to 1.

Optional arguments:
-debug
END
;

my ($queryFile, $dbFasta, $dbFile, $outFile, $debug);
die $usage
  unless GetOptions('query=s' => \$queryFile,
                    'db=s' => \$dbFasta,
                    'kmerdb=s' => \$dbFile,
                    'out=s' => \$outFile,
                    'debug' => \$debug)
  && @ARGV ==0
  && defined $queryFile
  && defined $dbFasta
  && defined $dbFile
  && defined $outFile;

open(my $fhQuery, "<", $queryFile) || die "Cannot read $queryFile\n";
open(my $fhFaa, "<", $dbFasta) || die "Cannot read $dbFasta\n";
die "No such file: $dbFile\n" unless -e $dbFile;
my $dbh = DBI->connect("dbi:SQLite:dbname=$dbFile","","",{ RaiseError => 1 }) || die $DBI::errstr;

my $usearch = "$RealBin/../bin/usearch";
die "No such executable: $usearch\n" unless -x $usearch;

my $state = {};
my ($queryHeader, $querySeq) = ReadFastaEntry($fhQuery, $state);
die "No input sequence\n" unless $querySeq;

my ($kmerLen) = $dbh->selectrow_array("SELECT LENGTH(kmer) FROM KmerLast limit 1;");
die "No kmer length\n" unless $kmerLen > 0;

# Strip leading methionines, since that is how the index works. Also strip any *
my $querySeqM = $querySeq;
$querySeqM =~ s/^[Mm]+//;
$querySeqM =~ s/[*]//g;
my $kmer1 = substr($querySeqM, 0, $kmerLen);
my $hits1 = $dbh->selectcol_arrayref("SELECT seek FROM KmerFirst WHERE kmer = ? ORDER BY seek", {}, $kmer1);

my $hits2 = [];
if (length($querySeqM) > $kmerLen) {
  my $kmer2 = substr($querySeqM, length($querySeqM) - $kmerLen);
  $hits2 = $dbh->selectcol_arrayref("SELECT seek FROM KmerLast WHERE kmer = ? ORDER BY seek", {}, $kmer2);
}

my %hits = ();
foreach my $seek (@$hits1) { $hits{$seek}++; }
foreach my $seek (@$hits2) { $hits{$seek}++; }

my @fetch = keys %hits;
if (@fetch >= 100) {
  # Select only those that match both kmers
  my %hits2 = map { $_ => 1 } @$hits2;
  @fetch = grep exists $hits2{$_}, @$hits1;
}
# Fetch in order to reduce seek times
@fetch = sort { $a <=> $b } @fetch;
print STDERR sprintf("Hits1 %d Hits2 %d Fetching %d\n",
                     scalar(@$hits1), scalar(@$hits2), scalar(@fetch))
  if defined $debug;

my @headers = ();
my @hitSeqs = ();
my %headerToSeq = ();
foreach my $i (0..(scalar(@fetch)-1)) {
  my $seek = $fetch[$i];
  seek($fhFaa, $seek, 0) || die "Cannot seek to $seek\n";
  my $header = <$fhFaa>;
  die "No header at $seek\n" unless defined $header && $header =~ m/^>/;
  chomp $header;
  $header =~ s/^>//;
  my $hitSeq = "";
  while(my $line = <$fhFaa>) {
    chomp $line;
    if ($line =~ m/^>/) {
      last;
    } else {
      $hitSeq .= $line;
    }
  }
  die "Empty seq at $seek with $header\n" unless $hitSeq;
  $headers[$i] = $header;
  $hitSeqs[$i] = $hitSeq;
  $headerToSeq{$header} = $hitSeq;
}

open(my $fhOut, ">", $outFile) || die "Cannot write to $outFile";
print $fhOut join("\t", @outFields)."\n";

if (@fetch == 0) {
  close($fhOut) || die "Error writing to $outFile\n";
  exit(0);
}

my $tmpDir = $ENV{TMPDIR} || "/tmp";
my $tmpPre = "$tmpDir/searchKmerSeekDb.$$";
open(my $fhHitFaa, ">", "$tmpPre.faa") || die "Cannot write to $tmpPre.faa\n";
foreach my $i (0..(scalar(@fetch)-1)) {
  print $fhHitFaa ">".$headers[$i]."\n";
  print $fhHitFaa $hitSeqs[$i] . "\n";
}
close($fhHitFaa) || die "Error writing to $tmpPre.hits\n";
my $cmd = "$usearch -search_local $queryFile -db $tmpPre.faa -id 0.8 -blast6out $tmpPre.hits -quiet";
print STDERR "Wrote tmpPre.faa\nRunning $cmd\n" if defined $debug;
system($cmd) == 0 || die "usearch failed:\n$cmd\n$!";



open (my $fhHits, "<", "$tmpPre.hits") || die "Cannot read $tmpPre.hits";
while (my $line = <$fhHits>) {
  chomp $line;
  my (undef, $header, $identity, $aLen, $mm, $nGaps, $qBeg, $qEnd, $sBeg, $sEnd) = split /\t/, $line;
  my $seq = $headerToSeq{$header};
  die "Unknown header from usearch in $tmpPre.hits"
    unless defined $seq;
  print $fhOut join("\t", $queryHeader, $header, $identity, $aLen, $mm, $nGaps, 
                    $qBeg, $qEnd, $sBeg, $sEnd, $seq) . "\n";
}
close($fhHits) || die "Error reading $tmpPre.hits\n";
close($fhOut) || die "Error writing to $outFile\n";
unlink("$tmpPre.faa") unless defined $debug;
unlink("$tmpPre.hits") unless defined $debug;
