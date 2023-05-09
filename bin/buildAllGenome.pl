#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use DBI;

my $usage = <<END
buildAllGenome.pl -dir data [ -test ]

Empties the AllGenome table from data/neighbor.db and repopulates it
based on the Genome table and all the sub-databases. Writes temporary
files to data/. In -test mode, does not alter the database.
END
;

my $dir;
my $test;

die $usage
  unless GetOptions('dir=s' => \$dir, 'test' => \$test)
  && @ARGV == 0
  && defined $dir;
die "Not a directory: $dir\n" unless -d $dir;
my $dbFile = "$dir/neighbor.db";
die "No such file: $dbFile\n" unless -e $dbFile;

my $dbhTop = DBI->connect("dbi:SQLite:dbname=$dbFile","","",{ RaiseError => 1 }) || die $DBI::errstr;
my $topGenomes = $dbhTop->selectall_arrayref("SELECT * from Genome", { Slice => {} });
my %gidToGenome = map { $_->{gid} => $_ } @$topGenomes;
print STDERR "Found " . scalar(keys %gidToGenome) . " genomes in the top level database\n";

foreach my $genome (values %gidToGenome) {
  $genome->{inTop} = 1;
  $genome->{prefix} = "";
}

my $subDbs = $dbhTop->selectcol_arrayref("SELECT prefix FROM SubDb");
$dbhTop->disconnect();

foreach my $subDb (@$subDbs) {
  my $subDbFile = "$dir/$subDb/sub.db";
  die "No such file: $subDbFile\n" unless -e $subDbFile;
  my $dbhSub = DBI->connect("dbi:SQLite:dbname=$subDbFile","","",{ RaiseError => 1 }) || die $DBI::errstr;
  my $subdbGenomes = $dbhSub->selectall_arrayref("SELECT * FROM Genome", { Slice => {} });
  foreach my $genome (@$subdbGenomes) {
    my $gid = $genome->{gid};
    if (exists $gidToGenome{$gid}) {
      print STDERR "Warning: $gid is in more than one SubDb\n"
        if $gidToGenome{$gid}{prefix} ne "";
      $gidToGenome{$gid}{prefix} = $subDb;
    } else {
      $genome->{inTop} = 0;
      $genome->{prefix} = $subDb;
      $gidToGenome{$gid} = $genome;
    }
  }
  $dbhSub->disconnect();
}
print STDERR "Found " . scalar(keys %gidToGenome) . " total genomes\n";

my $allFile = "$dir/AllGenome.tsv";
open(my $fhAll, ">", $allFile) || die "Cannot write to $allFile\n";
foreach my $gid (sort keys %gidToGenome) {
  my $genome = $gidToGenome{$gid};
  my @out = map $genome->{$_}, qw{gid
                                gtdbDomain gtdbPhylum gtdbClass gtdbOrder
                                gtdbFamily gtdbGenus gtdbSpecies
                                strain gtdbAccession assemblyName
                                ncbiTaxonomy nGenes nProteins inTop prefix};
  foreach my $value (@out) {
    die "Unset field in $gid" unless defined $value;
  };
  print $fhAll join("\t", @out)."\n";
}
close($fhAll) || die "Error writing to $allFile";

if (defined $test) {
  print STDERR "Wrote $allFile\n";
  print STDERR "Done (test mode)\n";
  exit(0);
}
# else
print STDERR "Clearing and loading AllGenome\n";
open(my $fhSql, "|-", "sqlite3", $dbFile) || die "Cannot run sqlite3 on $dbFile";
print $fhSql <<END
.mode tabs
DELETE FROM AllGenome;
.import $allFile AllGenome
SELECT 'AllGenome rows:', COUNT(*) FROM AllGenome;
END
;
close($fhSql) || die "sqlite3 failed: $!";
unlink($allFile);
print STDERR "Done\n";

