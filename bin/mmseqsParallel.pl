#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw{$RealBin};
use File::Copy;

my $mmseqs = "$RealBin/mmseqs";
my $nAlignFactor = 8;
my $limit = 300;
my $nCPUs = 10;
my $dbLoadMode = 2;
my $minPerJob = 20;
# These defaults are as in mmseqs itself
my $eValue = 0.001;
my $sens = 5.7;
my ($logFile, $debug);

sub run;

my $usage = <<END
Usage: mmseqsParallel.pl -in faa -db mmseqsdb -out hits

Runs mmseqs to find homologs for the input protein (fasta format) in
an mmseqs database, with the "align" step done in parallel. Intended
for single protein queries. The results should be identical to those
of mmseqs easy-search (except for the sorting of ties in both bits and
e-value) if you set -limit and -max-align to the same value as mmseq's
--max-seqs.

The mmseqs database can be created from a fasta file using mmseqs
createdb. For good performance, the database should be indexed (mmseqs
createindex) and loaded into memory (mmseqs touchdb).

Optional arguments:
-nCPUs $nCPUs -- how many CPUs to use
-sens $sens -- controls the sensitivity of mmseqs2
-eValue $eValue -- minimum evalue for results to report.
-limit $limit -- maximum number of hits to return
-max-align -- the maximum #alignments to consider per query
  defaults to limit * $nAlignFactor
-db-load-mode $dbLoadMode -- change to 0 if you did not run mmseqs touchdb
-mmseqs mmseqs_executable
  default: $mmseqs
-log logFile -- write logs from mmseqs to this file
-debug -- write debugging information and do not
  remove temporary files
END
;
my ($inFile, $dbFile, $outFile, $nMaxAlign);
die $usage
  unless GetOptions('in=s' => \$inFile,
                    'db=s' => \$dbFile,
                    'outFile=s' => \$outFile,
                    'mmseqs=s' => \$mmseqs,
                    'nCPUs=i' => \$nCPUs,
                    'eValue=f' => \$eValue,
                    'limit=i' => \$limit,
                    'max-align=i' => \$nMaxAlign,
                    'sens=f' => \$sens,
                    'db-load-mode=i' => \$dbLoadMode,
                    'log=s' => \$logFile,
                    'debug' => \$debug)
  && defined $inFile && defined $dbFile && defined $outFile;
die "No such executable: $mmseqs\n" unless -x $mmseqs;
require Time::HiRes if defined $debug;
foreach my $file ($inFile, $dbFile) {
  die "No such file: $file\n" unless -e $file;
}
die "No .idx file for $dbFile -- is this really an mmseqs database?\n" unless -e "$dbFile.idx";
die "Invalid eValue\n" unless $eValue > 0;
die "Invalid sens\n" unless $sens >= 1 && $sens <= 7.5;
die "Invalid nCPUs\n" unless $nCPUs >= 1;
$nCPUs = 50 if $nCPUs > 50;
die "Invalid limit\n" unless $limit > 0;
$nMaxAlign = $limit * $nAlignFactor if !defined $nMaxAlign;
die "Invalid max-align\n" unless $nMaxAlign > 0;

my $tmpDir = ($ENV{TMPDIR} || "$RealBin/../tmp") . "/mmseqsParallel.$$";
mkdir($tmpDir) || die "Cannot mkdir $tmpDir\n";

$logFile = "$tmpDir/log" if !defined $logFile;
open(STDOUT, ">", "$logFile") || die "Cannot write to $logFile\n";
print "Temporary directory: $tmpDir\n";
print STDERR "Logging to $logFile\n" if $debug;

# These steps mirror what mmseqs easy-search does, except that the
# results of the prefilter are split into nCPUs pieces

run($mmseqs, "createdb", $inFile, "$tmpDir/query",
       "--dbtype", "1", "--shuffle", "0", "--write-lookup", "0");
run($mmseqs, "prefilter", "$tmpDir/query", "$dbFile.idx", "$tmpDir/pref",
    "-s", $sens, "--max-seqs", $nMaxAlign, "--db-load-mode", $dbLoadMode, "--threads", "1");

# Read the prefilter candidates
my @prefilter = (); # a list of lines, except the final null string
open(my $fhPre, "<", "$tmpDir/pref") || die "Cannot read $tmpDir/pref";
while(my $line = <$fhPre>) {
  if ($line eq "\0") {
    last;
  } else {
    push @prefilter, $line;
  }
}
close($fhPre) || die "Error reading $tmpDir/pref";
my $nPrefilter = scalar(@prefilter);
print STDERR "Read " . $nPrefilter . " entries from prefilter\n" if defined $debug;
if ($nCPUs > $nPrefilter) {
  $nCPUs = $nPrefilter;
  $nCPUs = 1 if $nCPUs < 1;
}

for (my $iCPU = 0; $iCPU < $nCPUs; $iCPU++) {
  if (my $pid = fork()) {
    ; # do nothing in the parent process
  } else {
    # Write pref$i and .dbtype and .index files
    my $fh;
    # Write entries iCPU, iCPU + nCPUs, iCPU + 2*nCPUs... from pref to pref{$iCPU}
    open($fh, ">", "$tmpDir/pref${iCPU}") || die "Cannot write to $tmpDir/pref${iCPU}";
    for (my $i = $iCPU; $i < $nPrefilter; $i += $nCPUs) {
      print $fh $prefilter[$i];
    }
    print $fh "\0";
    close($fh) || die "Error running prefilter";
    # Write pref${iCPU}.index as tab-delimited 0, 0, filesize
    my $size = -s "$tmpDir/pref${iCPU}";
    open ($fh, ">", "$tmpDir/pref${iCPU}.index") || die "Cannot write to $tmpDir/pref${iCPU}.index";
    print $fh join("\t", 0, 0, $size)."\n";
    close($fh) || die "Error writing to $tmpDir/pref${iCPU}.index";
    # Copy pref.dbtype (which is small) to pref${iCPU}.dbtype
    copy("$tmpDir/pref.dbtype", "$tmpDir/pref${iCPU}.dbtype");
    run($mmseqs, "align", "$tmpDir/query", "$dbFile.idx", "$tmpDir/pref${iCPU}", "$tmpDir/result${iCPU}",
        "--alignment-mode", "3", "-e", $eValue, "--db-load-mode", $dbLoadMode, "--threads", "1");
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
die "mmseqs failed, see $logFile\n" unless $success;

# Combine the results into one result file, and limit to #limit entries
my @results = ();
for (my $iCPU = 0; $iCPU < $nCPUs; $iCPU++) {
  open(my $fhR, "<", "$tmpDir/result${iCPU}") || die "Cannot read $tmpDir/result${iCPU}";
  while(my $line = <$fhR>) {
    if ($line eq "\0") {
      last;
    } else {
      chomp $line;
      my @F = split /\t/, $line;
      die "Invalid line\n$line\n in $tmpDir/result${iCPU}" unless @F == 10;
      push @results, \@F;
    }
  }
  close($fhR) || die "Error reading $tmpDir/result{$iCPU}";
}
print STDERR "Read " . scalar(@results) . " results from align\n" if defined $debug;

# Sort and write out the alignment results database
# [3] is the e value, [1] is bit score,  and [0] is the targetId (integer)
#   (use evalue first because it has more precision)
@results = sort { $a->[3] <=> $b->[3]  || $b->[1] <=> $a->[1] || $a->[0] <=> $b->[0] } @results;
splice @results, $limit; # limit the result size
open(my $fhR, ">", "$tmpDir/result") || die "Cannot write to $tmpDir/result";
foreach my $result (@results) {
  print $fhR join("\t", @$result) . "\n";
}
print $fhR "\0";
close($fhR) || die "Error writing $tmpDir/result";
copy("$tmpDir/result0.dbtype", "$tmpDir/result.dbtype");
my $resultSize = -s "$tmpDir/result";
open(my $fhIndex, ">", "$tmpDir/result.index") || die "Cannot write $tmpDir/result.index";
print $fhIndex join("\t", 0, 0, $resultSize)."\n";
close($fhIndex) || die "Error writing to $tmpDir/result.index";

run($mmseqs, "convertalis", "$tmpDir/query", "$dbFile.idx", "$tmpDir/result", $outFile,
    "--db-load-mode", "2", "--search-type", "1", "--threads", "1");

if (defined $debug) {
 print STDERR "Done. Did not remove $tmpDir\n";
} else {
  # Ensure that tmpDir can be deleted
  open(STDOUT, ">", "/dev/null");
  system("rm -Rf $tmpDir");
}

sub run {
  my (@cmd) = @_;
  print join(" ", "Running: ", @cmd) . "\n";
  print STDERR join(" ", "Running: ", @cmd) . "\n" if defined $debug;
  my $t0 = [Time::HiRes::gettimeofday()] if defined $debug;
  system(@cmd) == 0
    || die "Command failed: $!\n" . join(" ", @cmd) . "\nSee $logFile\n";
  print STDERR "Took " . Time::HiRes::tv_interval($t0) . " s\n" if defined $debug;
}
