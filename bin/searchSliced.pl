#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};

my $nThreads = 1; # 1 thread per slice
# These defaults are as in mmseqs itself
my $eValue = 0.001;
my $maxSeqs = 300;
my $dbLoadMode = 1;
my $sens = 5.7;

my $usage = <<END
searchSliced.pl -in fasta -sliced prefix -out hits
Runs mmseqs easy-search on each slice (prefix0, prefix1, ...)
and combines the results.

Optional arguments:
-e $eValue
-s $sens -- controls the sensitivity of mmseqs2
-max-seqs $maxSeqs -- controls the pre-filter of mmseqs2,
  and also limits the maximum #sequences returned
-limit -- specify the number of hits to return
  (defaults to max-seqs)
-threads 1 -- #threads for each slice
  (each slice is run in parallel)
-db-load-mode $dbLoadMode -- use 2 if you ran touchSliced.pl
END
;

my ($inFaa, $pre, $outFile, $limit, $debug);
die $usage
  unless GetOptions('in=s' => \$inFaa, 'sliced=s' => \$pre, 'out=s' => \$outFile,
                    'debug' => \$debug,
                    's=f' => \$sens,
                    'e=f' => \$eValue,
                    'max-seqs=i' => \$maxSeqs,
                    'limit=i' => \$limit,
                    'threads=i' => \$nThreads,
                    'db-load-mode=i' => \$dbLoadMode)
  && defined $inFaa && defined $pre;
die "#threads per slice must be in 1 to 100\n" unless $nThreads >= 1 && $nThreads <= 100;
die "e value must be positive\n" unless $eValue > 0;
die "max-seqs must be positive\n" unless $maxSeqs > 0;
$limit = $maxSeqs if !defined $limit;
die "limit must be positive\n" unless $limit > 0;

die "No such file: $inFaa\n" unless -e $inFaa;

open (my $fh, "<", "$pre.nSlices") || die "Cannot read $pre.nSlices\n";
my $nSlices = <$fh>;
close($fh) || die "Error reading $pre.nSlices\n";
chomp $nSlices;
die "Invalid #slices in $pre.nSlices\n"
  unless $nSlices >= 1 && $nSlices <= 100;

my $mmseqs = "$RealBin/mmseqs";
system("$mmseqs -h > /dev/null") == 0
  || die "Cannot run $mmseqs\n";

my $tmp = $ENV{TMPDIR} || "/tmp";
die "Not a directory: $tmp\n" unless -d $tmp;
my $tmpDir = "$tmp/buildSliced.$$";
mkdir($tmpDir) || die "Cannot mkdir $tmpDir";

my $eValueSlice = $eValue / $nSlices;
my $maxSeqsSlice = int(0.5 + 1.1 * $maxSeqs / $nSlices);
$maxSeqsSlice = 1 if $maxSeqsSlice < 1;
for (my $i = 0; $i < $nSlices; $i++) {
  # Each job needs its own temporary directory
  my $tmpSub = "$tmpDir/$i";
  mkdir($tmpSub) || die "Cannot mkdir $tmpSub\n";
  my $cmd = "$mmseqs easy-search $inFaa $pre$i $tmpDir/hits$i $tmpSub --db-load-mode $dbLoadMode -e $eValueSlice -s $sens --max-seqs $maxSeqs --threads $nThreads > $tmpDir/log$i";
  print STDERR "Running: $cmd\n" if defined $debug;
  if (my $pid = fork()) {
    ; # do nothing in the parent process
  } else {
    if (system($cmd) != 0) {
      unlink("$tmpDir/hits$i");
      die "$cmd failed: $!";
    }
    exit(0);
  }
}

my $success = 1;
for(;;) {
  my $wait = wait;
  if ($wait < 0) {
    last;
  } else {
    my $err = $? >> 8;
    $success = 0 if $err;
  }
}

if (! $success) {
  system("rm -Rf $tmpDir");
  die "mmseqs failed\n";
}

open (my $fhOut, ">", $outFile) || die "Cannot write to $outFile\n";
print STDERR "All $nSlices mmseqs jobs succeeded\n" if defined $debug;
my @hits = ();
for (my $i = 0; $i < $nSlices; $i++) {
  open (my $fhHits, "<", "$tmpDir/hits$i") || die "Cannot read $tmpDir/hits$i";
  while (my $line = <$fhHits>) {
    chomp $line;
    my @F = split /\t/, $line;
    # adjust evalues
    $F[10] = sprintf("%.4g", $nSlices * $F[10]);
    push @hits, \@F;
  }
  close($fhHits) || die "Error reading $tmpDir/hits$i\n";
}
print STDERR "Read all hits\n" if defined $debug;

# sort by query (index 0) and then by bits (index 11) descending
@hits = sort { $a->[0] cmp $b->[0] || $b->[11] <=> $a->[11] } @hits;
# Reduce length to maxSeqs
my %nShown = (); # query to #hits shown so far
foreach my $row (@hits) {
  my $query = $row->[0];
  if (!exists $nShown{$query} || $nShown{$query} < $limit) {
    $nShown{$query}++;
    print $fhOut join("\t",  @$row) . "\n";
  }
}
close($fhOut) || die "Error writing to $outFile\n";
system("rm -Rf $tmpDir");
print STDERR "Wrote $outFile\n" if defined $debug;
