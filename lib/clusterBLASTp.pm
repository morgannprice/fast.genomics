# Cluster a set of protein sequences and then
# search against them.
package clusterBLASTp;
require Exporter;
use strict;
use DBI;
our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw{cluster parseClustering expandByClusters clusteringIntoDb
             proteinsToFaa formatBLASTp runBLASTp proteinDbSize
             parseBLASTpHits saveBLASTpHits
             clusteredBLASTp
             removeDuplicates
          };

# Writes the reduced fasta to out and the clustering to out.clstr
# Saves the cd-hit log in out.log
# The input is a hash, which must include faa, cdhit, and out
sub cluster {
  my (%param) = @_;
  my $cdhit = $param{cdhit} || die "Must specify cdhit";
  die "Not an executable: $cdhit\n" unless -x $cdhit;
  my $faa = $param{faa};
  die "Must specify input faa" unless defined $param{faa};
  my $out = $param{out};
  die "Must specify out" unless defined $out && $out ne "";
  my $nCPUs = $param{CPUs} || $ENV{MC_CORES} || (`egrep -c '^processor' /proc/cpuinfo`);
  my $minIdentity = $param{minIdentity} || 0.8;
  my $minCoverage = $param{minCoverage} || 0.8;

  my $clusterCmd = "$cdhit -i $faa -o $out -c $minIdentity -M 0 -T $nCPUs -n 5 -d 0"
    . " -aS $minCoverage -aL $minCoverage >& $out.log";
  system($clusterCmd) == 0 || die "cd-hit failed:\n$clusterCmd\n$!\n";
  return parseClustering("$out.clstr");
}

# Returns a reference to a hash of clusterId => list of members,
# where the clusterId is the id of the representative sequence
sub parseClustering($) {
  my ($clusterFile) = @_;
  my @clusters = ();
  my @clusterIds = ();
  my %protSeen = ();
  my $hasClusterId = 0;
  open(my $fh, "<", $clusterFile) || die "Cannot read $clusterFile\n";
  while (my $line = <$fh>) {
    if ($line =~ m/^>/) {
      die "No representative for cluster\n"
        if ! $hasClusterId && @clusters > 0;
      push @clusters, []; # start an emtpy cluster
      $hasClusterId = 0;
    } else {
      $line =~ m/>(\S+)[.][.][.] (.)/ || die "Cannot parse cluster line\n$line\nin $clusterFile\n";
      my $proteinId = $1;
      my $rep = $2 eq "*";
      die "Protein $proteinId occurs more than once in\n$clusterFile\n"
        if exists $protSeen{$proteinId};
      $protSeen{$proteinId} = 1;
      die "Cluster member before cluster number\n" if @clusters == 0;
      push @{ $clusters[-1] }, $proteinId;
      if ($rep) {
        die "Cluster already has a representative sequence\n$line\nin $clusterFile\n"
          if $hasClusterId;
        $hasClusterId = 1;
        push @clusterIds, $proteinId;
      }
    }
  }
  close($fh) || die "Error reading $clusterFile\n";
  die "Error parsing $clusterFile" if scalar(@clusters) != scalar(@clusterIds);
  my %clusters = ();
  foreach my $i (0..(scalar(@clusters)-1)) {
    $clusters{ $clusterIds[$i] } = $clusters[$i];
  }
  return \%clusters;
}

# Import the clusering into the ClusterProtein table
# (The table must exist already, and it adds to the table.)
sub clusteringIntoDb($$) {
  my ($clusters, $dbFile) = @_;
  die "No such file: $dbFile\n" unless -e $dbFile;
  my $tmpFile = ($ENV{TMPDIR} || "/tmp") . "/clusterBLASTp.$$";
  open (my $fh, ">", $tmpFile) || die "Cannot write to $tmpFile\n";
  foreach my $clusterId (sort keys %$clusters) {
    foreach my $proteinId (@{ $clusters->{$clusterId} }) {
      print $fh join("\t", $clusterId, $proteinId)."\n";
    }
  }
  close($fh) || die "Error writing to $tmpFile\n";
  open(my $fhSql, "|-", "sqlite3", "$dbFile") || die "Cannot run sqlite3 on $dbFile";
  print $fhSql ".mode tabs\n",
    ".import $tmpFile ClusterProtein\n";
  close($fhSql) || die "Error loading $tmpFile into $dbFile\n$!";
  unlink($tmpFile);
}

# The directory should contain makeblastdbd, if $blastPlus, or formatdb otherwise.
# Uf $blastPlus, puts the database in $db.plusdb
sub formatBLASTp($$$) {
  my ($exeDir, $faa, $blastPlus) = @_;
  die "Not a directory: $exeDir\n" unless -d $exeDir;
  die "No such file: $faa\n" unless -e $faa;
  my $formatCmd;
  if ($blastPlus) {
    die "No such executable: $exeDir/makeblastdb\n" unless -x "$exeDir/makeblastdb";
    $formatCmd = "$exeDir/makeblastdb -in $faa -dbtype prot -out $faa.plusdb -logfile /dev/null";
  } else {
    die "No such executable: $exeDir/formatdb\n" unless -x "$exeDir/formatdb";
    $formatCmd = "$exeDir/formatdb -i $faa -p T";
  }
  system($formatCmd) == 0 || die "makeblastdb failed:\n$formatCmd\n$!";
}

# Required members of the input hash:
# exeDir (containing either original BLAST blastall or BLAST+'s blastp),
# query (a sequence),
# db (the fasta database).
# Optional arguments:
#  blastPlus -- if not specifiedif db.plusdb.pin exists, then
#    it uses blast+, otherwise it uses blastall.
#  eValue -- default 1e-3
#  nCPUs -- default 8
#  dbSize -- default 0 (the actual size)
#  maxHits -- defaults to 10,000
#
# Returns a reference to list of hashes, each with subject, identity (as a
# fraction, not a percentage), alnLength, nMismatch, nGapOpens,
# qBegin, qEnd, sBegin, sEnd, eValue, bits
sub runBLASTp {
  my (%param) = @_;
  my $exeDir = $param{exeDir} || die "Must specify exeDir";
  die "Not a directory: $exeDir" unless -d $exeDir;
  my $query = $param{query} || die "Must specify query";
  die "Invalid query sequence" unless $query =~ m/^[A-Z]+$/;
  my $db = $param{db};
  die "Must specify db" if !defined $db || $db eq "";
  my $maxHits = $param{maxHits};
  $maxHits = 10000 if !defined $maxHits;
  die "Invalid maxHits" unless $maxHits =~ m/^\d+$/;

  my $blastPlus = $param{blastPlus};
  $blastPlus = -e "$db.plusdb.pin" if !defined $blastPlus;

  if ($blastPlus) {
    die "No such executable: $exeDir/blastp" unless -x "$exeDir/blastp";
    die "No such file: $db.plusdb.pin" unless -e "$db.plusdb.pin";
  } else {
    die "No such executable: $exeDir/blastall" unless -x "$exeDir/blastall";
    die "No such file: $db.pin" unless -e "$db.pin";
  }

  my $minEValue = $param{eValue} || 1e-3;
  my $nCPUs = $param{nCPUs} || 8;
  my $dbSize = $param{dbSize} || 0;

  my $tmpPre = ($ENV{TMPDIR} || "/tmp") . "/clusterBLASTp.$$";
  my $tmpFaa = "$tmpPre.faa";
  open (my $fhFaa, ">", $tmpFaa) || die "Cannot write to $tmpFaa";
  print $fhFaa ">query\n$query\n";
  close($fhFaa) || die "Error writing to $tmpFaa";

  my $tmpHits = "$tmpPre.hits";
  my $blastCmd;
  if ($blastPlus) {
    $blastCmd = "$exeDir/blastp -evalue $minEValue -num_threads $nCPUs -query $tmpFaa -db $db.plusdb"
      . " -max_target_seqs $maxHits -outfmt 6 -out $tmpHits";
    $blastCmd .= " -dbsize $dbSize" if $dbSize;
    $blastCmd .= " >& /dev/null";
  } else {
    $blastCmd = qq{$exeDir/blastall -p blastp -e $minEValue -z $dbSize -a $nCPUs -i $tmpFaa -d $db -F "m S" -m 8}
      . " -b $maxHits -v $maxHits -o $tmpHits >& /dev/null";
  }
  die "Run: $blastCmd\n" if ! $blastPlus;
  system($blastCmd) == 0 || die "runBLASTp() failed:\n$blastCmd\n$!";
  unlink($tmpFaa);

  my $hits = parseBLASTpHits($tmpHits);
  unlink($tmpHits);
  return($hits);
}

# Returns a reference to a list of rows, each as a hash
# This also works on last BLAST-like format, by skipping lines that start with #
sub parseBLASTpHits($) {
  my ($file) = @_;
  my @out = ();
  open(my $fhHits, "<", $file) || die "Cannot read $file";
  while (my $line = <$fhHits>) {
    chomp $line;
    next if $line eq "#" || $line =~ m/^# /;
    my ($query, $subject, $identity, $alen, $mm, $gap, $qBegin, $qEnd, $sBegin, $sEnd, $eValue, $bits)
      = split /\t/, $line;
    die "Cannot parse\n$line\nfrom blastp" unless defined $bits && $bits =~ m/\d/;
    $bits =~ s/^ *//;
    push @out, { 'query' => $query, 'subject' => $subject,
                 'identity' => $identity/100,
                 'alnLength' => $alen, 'nMismatch' => $mm, 'nGapOpens' => $gap,
                 'qBegin' => $qBegin, 'qEnd' => $qEnd,
                 'sBegin' => $sBegin, 'sEnd' => $sEnd,
                 'eValue' => $eValue, 'bits' => $bits };
  }
  close($fhHits) || die "Error reading $file";
  return \@out;
}

sub saveBLASTpHits($$) {
  my ($hits, $file) = @_;
  open(my $fh, ">", $file) || die "Cannot write to $file";
  foreach my $row (@$hits) {
    print $fh join("\t", "query", $row->{subject},
                   $row->{identity} * 100,
                   map $row->{$_}, qw{alnLength nMismatch nGapOpens
                                      qBegin qEnd
                                      sBegin sEnd
                                      eValue bits})."\n";
  }
  close($fh) || die "Error writing to $file";
}

# Given a reference to a list of clusterIds, and a handle to a database with the ClusterProtein table,
# expand to a list of all members. Returns a reference to a list.
sub expandByClusters($$) {
  my ($clusterIds, $dbh) = @_;
  my @out = ();
  foreach my $clusterId (@$clusterIds) {
    my $proteinIds = $dbh->selectcol_arrayref("SELECT proteinId FROM ClusterProtein WHERE clusterId = ?",
                                              {}, $clusterId);
    die "Cannot query ClusterProtein" unless defined $proteinIds;
    die "No members for cluster $clusterId"
      unless @$proteinIds > 0;
    push @out, @$proteinIds;
  }
  return \@out;
}

# Return the total number of amino acids in the fasta file
sub proteinDbSize($) {
  my ($faa) = @_;
  open(my $fh, "<", $faa) || die "Cannot read $faa\n";
  my $n = 0;
  while (my $line = <$fh>) {
    next if $line =~ m/^>/;
    chomp $line;
    $n += length($line);
  }
  close($fh) || die "Error reading $faa\n";
  return $n;
}

# Given a list of proteins and a database handle with the Protein table,
# make a fasta file of their sequences
# Returns the total #a.a.
sub proteinsToFaa($$$) {
  my ($proteinIds, $faaOut, $dbh) = @_;
  my $nAA = 0;
  open(my $fh, ">", $faaOut) || die "Cannot write to $faaOut";
  foreach my $proteinId (@$proteinIds) {
    my ($seq) = $dbh->selectrow_array("SELECT sequence FROM Protein WHERE proteinId = ?",
                                      {}, $proteinId);
    die "No sequence for $proteinId" unless $seq;
    $nAA += length($seq);
    print $fh ">" . $proteinId . "\n" . $seq . "\n";
  }
  close($fh) || die "Error writing to $faaOut";
  return $nAA;
}

# The input hash must include the arguments
# query (a sequence),
# clusterDb (the BLAST+ database file name, without the .plusdb suffix),
# dbh (the database handle, for expanding clusters and fetching protein sequences),
# nCPUs,
# bin (the directory with lastdb, lastall, and subdirectories blast+ or blast),
# maxHits (a list with the maximum #hits for the first round and for the second round)
#
# Optional argument:
# dbSize -- amino acids in the full (unclustered) database, or 0 to not correct evalues
# eValue -- default 1e-3
# quiet -- if not set, outputs HTML comments about what it is doing.
#
# Returns a list of hits, as from parseBLASTpHits
sub clusteredBLASTp {
  my (%in) = @_;
  my $query = $in{query} || die "Must specify query";
  die "Invalid query" unless $query =~ m/^[A-Z]+$/;
  my $clusterDb = $in{clusterDb} || die "Must specify clusterDb";
  my $dbh = $in{dbh} || die "Must specify dbh";
  my $nCPUs = $in{nCPUs} || die "Must specify nCPUs > 0";
  my $bin = $in{bin} || die "Must specify bin";
  my $maxHits = $in{maxHits} || die "Must specify maxHits";
  my ($nMaxHits1, $nMaxHits2) = @$maxHits;
  die "Invalid maxHits" unless $nMaxHits1 > 0 && $nMaxHits2 > 0;
  my $dbSize = $in{dbSize} || 0;
  my $eValue = $in{eValue} || 1e-3;
  my $quiet = $in{quiet} || 0;

  print "<P><small>Running BLASTp.</small>\n" unless $quiet;
  my $hits1 = runBLASTp('exeDir' => "$bin/blast+",
                        'blastPlus' => 1,
                        'query' => $query, 'db' => $clusterDb,
                        'eValue' => $eValue, 'maxHits' => $nMaxHits1, 'nCPUs' => $nCPUs);
  return [] if @$hits1 == 0;
  my @clusterIds = removeDuplicates(map $_->{subject}, @$hits1);
  splice @clusterIds, $nMaxHits1;
  print "<small>Processing " . scalar(@clusterIds) . " clusters and running lastal.</small>\n"
      unless $quiet;
  my $proteinIds = expandByClusters(\@clusterIds, $dbh);
  splice @$proteinIds, $nMaxHits2;
  print "<!-- expanded to " . scalar(@$proteinIds) . " -->\n" unless $quiet;
  my $tmpPre = ($ENV{TMPDIR} || "/tmp") . "/clusterBLASTp.$$";
  my $nAA = proteinsToFaa($proteinIds, "$tmpPre.faa", $dbh);
  print "<!-- fetched protein sequences -->\n" unless $quiet;

  my $inFile = "$tmpPre.in";
  open(my $fh, ">", $inFile) || die "Cannot write to $inFile";
  print $fh ">query\n$query\n";
  close($fh) || die "Error writing $inFile";

  my $lastal = "$bin/lastal";
  my $lastdb = "$bin/lastdb";
  foreach my $x ($lastal, $lastdb) {
    die "No such executable: $x\n" unless -x $x;
  }
  my $maxCand = $nMaxHits2;
  $maxCand = scalar(@$proteinIds) if scalar(@$proteinIds) < $nMaxHits2;
  my @lastCmds = ("$lastdb -P $nCPUs -p $tmpPre.lastdb $tmpPre.faa",
                  "$lastal -m $maxCand -f BlastTab -P $nCPUs $tmpPre.lastdb $inFile > $tmpPre.hits");
  foreach my $cmd (@lastCmds) {
    system($cmd) == 0
      || die "last failed:\n$cmd\n$!";
  }
  my $hits = parseBLASTpHits("$tmpPre.hits");
  unlink("$tmpPre.faa");
  unlink("$inFile");
  unlink("$tmpPre.hits");
  foreach my $suffix (qw{bck des prj sds ssp suf tis}) {
    unlink("$tmpPre.lastdb.$suffix");
  }
  # Limit the number of hits
  splice @$hits, $nMaxHits2;
  if ($dbSize > 0) {
    # Scale the evalues
    my $scale = $dbSize / $nAA;
    foreach my $hit (@$hits) {
      $hit->{eValue} = sprintf("%.0g", $scale * $hit->{eValue});
    }
  }
  my @hits2 = grep $_->{eValue} <= $eValue, @$hits;
  return \@hits2;
}

# Remove duplicates from a list, and return the new list
sub removeDuplicates {
  my (@in) = @_;
  my %seen = ();
  my @out = ();
  foreach my $value (@in) {
    next if exists $seen{$value};
    $seen{$value} = 1;
    push @out, $value;
  }
  return @out;
}

1;
