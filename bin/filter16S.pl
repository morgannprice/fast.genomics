#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use DBI;
use FindBin qw{$RealBin};
use lib "$RealBin/../../PaperBLAST/lib";
use pbutils qw{ReadFastaEntry};

sub filter($$$);
my $minLength = 250;

my $usage = <<END
filter16S.pl -fna 16S.fna -dir data > 16Sfiltered.fna

Given a collection of 16S sequences, such as from get16S.pl, with
sequence names of the form genomeId:locusTag, filter it to keep only
those that are in the databases in data/neighbor.db and its
sub-databases.

Optional arguments:
-minLength $minLength -- the minimum length of a 16S gene to include
END
  ;

my ($dir, $fna);
die $usage
  unless GetOptions('fna=s' => \$fna, 'dir=s' => \$dir,
                    'minLength=i' => \$minLength)
  && @ARGV == 0
  && defined $dir && defined $fna;
die "No such file: $fna\n" unless -e $fna;
die "No such directory: $dir\n" unless -d $dir;

my $dbhTop = DBI->connect("dbi:SQLite:dbname=$dir/neighbor.db","","",{ RaiseError => 1 }) || die $DBI::errstr;
my $allgenomes = $dbhTop->selectall_hashref("SELECT * from AllGenome", "gid");

my %bygid = (); # gid => locusTag => sequence
open(my $fh, "<", $fna) || die "Cannot read $fna\n";
my $state = {};
my $nNoTag = 0;
my $nShort = 0;
my $nKeep = 0;
while (my ($id, $seq) = ReadFastaEntry($fh, $state)) {
  $id =~ s/ .*//;
  my ($gid, $locusTag) = split /:/, $id;
  die "Invalid identifier $id\n" unless defined $locusTag;
  if ($locusTag eq "") {
    $nNoTag++;
  } elsif (!exists $allgenomes->{$gid}) {
    print STDERR "Skipping $id -- unknown gid\n";
  } elsif (exists $bygid{$gid}{$locusTag}) {
    print STDERR "Skipping duplicate $id\n";
  } elsif (length($seq) < $minLength) {
    $nShort++;
  } else {
    $bygid{$gid}{$locusTag} = $seq;
    $nKeep++;
  }
}
close($fh) || die "Error reading $fna\n";
print STDERR "Skipped $nShort short sequences, $nNoTag with no locus tags, and kept $nKeep\n";


my @gids = grep exists $bygid{$_}, sort keys %$allgenomes;
my @gidsTop = grep $allgenomes->{$_}{inTop}, @gids;
my @gidsNotop = grep ! $allgenomes->{$_}{inTop}, @gids;

foreach my $gid (@gidsTop) {
  filter($gid, $dbhTop, $bygid{$gid});
}

my %gidsByPrefix = ();
foreach my $gid (@gidsNotop) {
  my $g = $allgenomes->{$gid} || die $gid;
  push @{ $gidsByPrefix{$g->{prefix}} }, $gid;
}

foreach my $prefix (sort keys %gidsByPrefix) {
  my $subdbh = DBI->connect("dbi:SQLite:dbname=$dir/$prefix/sub.db","","",{ RaiseError => 1 }) || die $DBI::errstr;
  foreach my $gid (@{ $gidsByPrefix{$prefix} }) {
    filter($gid, $subdbh, $bygid{$gid});
  }
  $subdbh->disconnect();
}

# Output 16S if the locus tag is known, and warn otherwise
sub filter($$$) {
  my ($gid, $dbh, $hash) = @_;
  foreach my $locusTag (sort keys %$hash) {
    my $gene = $dbh->selectrow_hashref("SELECT * from Gene WHERE gid = ? AND locusTag = ?",
                                       {}, $gid, $locusTag);
    if (defined $gene) {
      print ">$gid:$locusTag\n" . $hash->{$locusTag} . "\n";
    } else {
      print STDERR "Skipping unknown locusTag $gid:$locusTag\n";
    }
  }
}
