#!/usr/bin/perl -w
use strict;
use DBI;
use Getopt::Long;
use FindBin qw{$RealBin};

my $dataDir = "$RealBin/../data";
my $fraction = 0.0002;

my $usage = <<END
touchSubDbs.pl -subdb subdb
touchSubDbs.pl -all

For each sub database, touchSubDbs.pl puts the ClusterProtein tables
into memory, and loads a small fraction of the proteins to keep the
index of the Protein table in memory

Optional arguments:
-data data_directory
  default: $dataDir
-fraction $fraction -- fraction of proteins to load
END
;

my ($subDbSpec, $all, $debug);
die $usage
  unless GetOptions('all' => \$all, 'subdb=s' => \$subDbSpec,
                    'fraction=f' => \$fraction,
                    'data=s' => \$dataDir,
                    'debug' => \$debug)
    && @ARGV == 0;
die "No such directory: $dataDir\n" unless -d $dataDir;
die "Invalid fraction\n" unless $fraction >= 0 && $fraction <= 1;
die "Must specify -subdb or -all\n"
  unless defined $all || defined $subDbSpec;

my @subDbs = ();
if (defined $all) {
  my $topDb = "$dataDir/neighbor.db";
  my $topDbh = DBI->connect("dbi:SQLite:dbname=$topDb", "", "", { RaiseError => 1 })
    || die $DBI::errstr;
  push @subDbs, @{  $topDbh->selectcol_arrayref("SELECT prefix FROM SubDb;") };
  $topDbh->disconnect()
} else {
  @subDbs = ( $subDbSpec );
}
print STDERR "Touching " . scalar(@subDbs) . " subdbs with protein fraction = $fraction\n";

foreach my $subDb (@subDbs) {
  die "No such directory: $dataDir/$subDb\n" unless -d "$dataDir/$subDb";
  my $dbFile = "$dataDir/$subDb/sub.db";
  my $dbh = DBI->connect("dbi:SQLite:dbname=$dbFile", "", "", { RaiseError => 1 })
    || die $DBI::errstr;
  my $clusters = $dbh->selectall_arrayref("SELECT proteinId,clusterId from ClusterProtein;");
  print "Loaded clusters from $dbFile\n" if defined $debug;
  foreach my $row (@$clusters) {
    if (rand() < $fraction) {
      my $rows = $dbh->selectall_arrayref("SELECT sequence from Protein WHERE proteinId = ?",
                                          {}, $row->[0]);
    }
  }
  print "Loaded proteins from $dbFile\n" if defined $debug;
  $dbh->disconnect();
}
print STDERR "Finished touching\n";
