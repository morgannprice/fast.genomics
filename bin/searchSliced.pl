#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};

my $nThreads = 1; # 1 thread per slice
# These defaults are as in mmseqs itself
my $eValue = 0.001;
my $limit = 300;
my $nAlignFactor = 10;
my $dbLoadMode = 1;
my $sens = 5.7;

my $usage = <<END
searchSliced.pl -in fasta -sliced prefix -out hits
Runs mmseqs easy-search on each slice (prefix0, prefix1, ...)
and combines the results.

Optional arguments:
-s $sens -- controls the sensitivity of mmseqs2
-e $eValue -- minimum evalue for results to report.
   (Each slice is run at E <= eValue / nSlices, and then those e values
    are then multiplied by the number of slices, and the hit is
    discarded if the corrected e value is above the threshold)
-limit $limit -- maximum number of hits to return
-max-align -- the maximum #alignments to consider per query
  defaults to limit * $nAlignFactor
  (Each slice is run with a limit of max-align / nSlices)
-threads $nThreads -- #threads for each slice
  (the slices are run in parallel)
-db-load-mode $dbLoadMode -- use 2 if you ran touchSliced.pl
END
;

my ($inFaa, $pre, $outFile, $maxAlign, $debug);
die $usage
  unless GetOptions('in=s' => \$inFaa, 'sliced=s' => \$pre, 'out=s' => \$outFile,
                    'debug' => \$debug,
                    's=f' => \$sens,
                    'e=f' => \$eValue,
                    'limit=i' => \$limit,
                    'max-align=i' => \$maxAlign,
                    'threads=i' => \$nThreads,
                    'db-load-mode=i' => \$dbLoadMode)
  && defined $inFaa && defined $pre;
die "#threads per slice must be in 1 to 100\n" unless $nThreads >= 1 && $nThreads <= 100;
die "e value must be positive\n" unless $eValue > 0;
die "limit must be positive\n" unless $limit > 0;
$maxAlign = $nAlignFactor * $limit if !defined $maxAlign;
die "max-align must be >= limit\n" unless $maxAlign >= $limit;

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

my $maxAlignSlice = int(0.5 + $maxAlign  / $nSlices);
$maxAlignSlice = 1 if $maxAlignSlice < 1;
my $eValueSlice = $eValue / $nSlices;
for (my $i = 0; $i < $nSlices; $i++) {
  # Each job needs its own temporary directory
  my $tmpSub = "$tmpDir/$i";
  mkdir($tmpSub) || die "Cannot mkdir $tmpSub\n";
  my $cmd = "$mmseqs easy-search $inFaa $pre$i $tmpDir/hits$i.tmp $tmpSub -e $eValueSlice --db-load-mode $dbLoadMode -s $sens --max-seqs $maxAlignSlice --threads $nThreads > $tmpDir/log$i";
  print STDERR "Running: $cmd\n" if defined $debug;
  if (my $pid = fork()) {
    ; # do nothing in the parent process
  } else {
    system($cmd) == 0 || die "$cmd failed: $!";
    rename("$tmpDir/hits$i.tmp", "$tmpDir/hits$i");
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

if (defined $debug) {
  for(my $i = 0; $i < $nSlices; $i++) {
    my $nLines = `wc -l < $tmpDir/hits$i`;
    chomp $nLines;
    print STDERR "slice $i hits $nLines\n";
  }
}

open (my $fhOut, ">", $outFile) || die "Cannot write to $outFile\n";
print STDERR "All $nSlices mmseqs jobs succeeded\n" if defined $debug;
my @hits = ();
for (my $i = 0; $i < $nSlices; $i++) {
  open (my $fhHits, "<", "$tmpDir/hits$i") || die "Cannot read $tmpDir/hits$i";
  while (my $line = <$fhHits>) {
    chomp $line;
    my @F = split /\t/, $line;
    push @hits, \@F;
  }
  close($fhHits) || die "Error reading $tmpDir/hits$i\n";
}
print STDERR "Read all hits\n" if defined $debug;

# sort by query (index 0) and then by bits (index 11) descending
@hits = sort { $a->[0] cmp $b->[0] || $b->[11] <=> $a->[11] } @hits;
my %nShown = (); # query to #hits shown so far
foreach my $row (@hits) {
  my $query = $row->[0];
  # Reduce #hits to limit when printing
  if (!exists $nShown{$query} || $nShown{$query} < $limit) {
    $nShown{$query}++;
    # adjust eValue and check it
    $row->[10] = sprintf("%.4g", $nSlices * $row->[10]);
    print $fhOut join("\t",  @$row) . "\n" if $row->[10] <= $eValue;
  }
}
close($fhOut) || die "Error writing to $outFile\n";
system("rm -Rf $tmpDir");
print STDERR "Wrote $outFile\n" if defined $debug;
