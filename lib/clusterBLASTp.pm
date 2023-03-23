# Cluster a set of protein sequences and then
# search against them.
package clusterBLASTp;
require Exporter;
use strict;
our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw{cluster parseClustering expandByClusters clusteringIntoDb
             proteinsToFaa formatBLASTp runBLASTp proteinDbSize};

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

sub formatBLASTp($$) {
  my ($formatdb, $faa) = @_;
  die "No such executable: $formatdb\n" unless -x $formatdb;
  die "No such file: $faa\n" unless -e $faa;
  my $formatCmd = "$formatdb -p T -i $faa";
  system($formatCmd) == 0 || die "formatdb failed:\n$formatCmd\n$!";
}

# The input parameters must include blastall (the path tot he
# executable), query (a sequence), and db (the fasta database) Returns
# a reference to list of hashes, each with subject, identity (as a
# fraction, not a percentage), alnLength, nMismatch, nGapOpens,
# qBegin, qEnd, sBegin, sEnd, eValue, bits
sub runBLASTp {
  my (%param) = @_;
  my $blastall = $param{blastall} || die "Must specify blastall";
  die "No such executable: $blastall\n" unless -x $blastall;
  my $query = $param{query} || die "Must specify query";
  die "Invalid query sequence" unless $query =~ m/^[A-Z]+$/;
  my $db = $param{db};
  die "Must specify db" if !defined $db || $db eq "";
  die "No such file: $db.pin" unless -e "$db.pin";
  my $minEValue = $param{eValue} || 1e-3;
  my $nCPUs = $param{nCPUs} || 8;
  my $dbSize = $param{dbSize} || 0;

  my $tmpPre = ($ENV{TMPDIR} || "/tmp") . "/clusterBLASTp.$$";
  my $tmpFaa = "$tmpPre.faa";
  open (my $fhFaa, ">", $tmpFaa) || die "Cannot write to $tmpFaa";
  print $fhFaa ">query\n$query\n";
  close($fhFaa) || die "Error writing to $tmpFaa";

  my $tmpHits = "$tmpPre.hits";
  my $blastCmd = qq{$blastall -p blastp -e $minEValue -z $dbSize -a $nCPUs -i $tmpFaa -d $db -F "m S" -m 8}
    . " -b 10000 -v 10000 -o $tmpHits >& /dev/null";
  system($blastCmd) == 0 || die "blastall failed:\n$blastCmd\n$!";
  unlink($tmpFaa);

  my @out = ();
  open(my $fhHits, "<", $tmpHits) || die "Cannot read $tmpHits";
  while (my $line = <$fhHits>) {
    chomp $line;
    my (undef, $subject, $identity, $alen, $mm, $gap, $qBegin, $qEnd, $sBegin, $sEnd, $eValue, $bits)
      = split /\t/, $line;
    die "Cannot parse\n$line\nfrom blastp" unless defined $bits && $bits =~ m/\d/;
    $bits =~ s/^ *//;
    push @out, { 'subject' => $subject, 'identity' => $identity/100,
                 'alnLength' => $alen, 'nMismatch' => $mm, 'nGapOpens' => $gap,
                 'qBegin' => $qBegin, 'qEnd' => $qEnd,
                 'sBegin' => $sBegin, 'sEnd' => $sEnd,
                 'eValue' => $eValue, 'bits' => $bits };
  }
  close($fhHits) || die "Error reading $tmpHits";
  unlink($tmpHits);
  return(\@out);
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
sub proteinsToFaa($$$) {
  my ($proteinIds, $faaOut, $dbh) = @_;
  open(my $fh, ">", $faaOut) || die "Cannot write to $faaOut";
  foreach my $proteinId (@$proteinIds) {
    my ($seq) = $dbh->selectrow_array("SELECT sequence FROM Protein WHERE proteinId = ?",
                                      {}, $proteinId);
    die "No sequence for $proteinId" unless $seq;
    print $fh ">" . $proteinId . "\n" . $seq . "\n";
  }
  close($fh) || die "Error writing to $faaOut";
}

1;
