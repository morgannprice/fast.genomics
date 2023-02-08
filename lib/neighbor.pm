# Utilities for the gene neighbor viewer
package neighbor;
require Exporter;
use strict;

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(parseTaxString readMMSeqsHits estimateTopScore cumsum);

# Returns a hash of d/p/c/o/f/g/s to value
# (any of which could be missing),
# or undef if it cannot parse it
sub parseTaxString($) {
  my ($tax) = @_;
  my @pieces = split /;/, $tax;
  my %taxa = ();
  foreach my $piece (@pieces) {
    $piece =~ m/^([a-z])__(.*)$/ || return undef;
    $taxa{$1} = $2;
  }
  return \%taxa;
}

# Given a file with mmseqs hits, returns a reference to a list of
# hashes, each with query, subject, identity (as a fraction, not a
# percentage), alnLength, nMismatch, nGapOpens, qBegin, qEnd, sBegin,
# sEnd, eValue, bits
sub readMMSeqsHits($) {
  my ($file) = @_;
  open (my $fh, "<", $file) || die "Cannot read $file";
  my @colNames = qw{query subject identity alnLength nMismatch nGapOpens qBegin qEnd sBegin sEnd eValue bits};
  my @out;
  while (my $line = <$fh>) {
    chomp $line;
    my @F = split /\t/, $line;
    die "Wrong number of columns in mmseqs hits file $file"
      unless scalar(@F) == scalar(@colNames);
    my %row = ();
    foreach my $i (0..(scalar(@F)-1)) {
      die "Empty column in mmseqs hits file $file" unless $F[$i] =~ m/\S/;
      $row{ $colNames[$i] } = $F[$i];
    }
    foreach my $field (qw{alnLength nMismatch nGapOpens qBegin qEnd sBegin sEnd}) {
      die "Invalid $field in mmseqs hits file $file"
        unless $row{$field} =~ m/^\d+$/;
    }
    die "Invalid bits in mmseqs file $file"
      unless $row{bits} =~ m/^[0-9.]+$/;
    die "Invalid E-value in mmseqs file $file"
      unless $row{eValue} =~ m/^[0-9.eE+-]+$/;
    push @out, \%row;
  }
  close($fh) || die "Error reading $file";
  return \@out;
}

# Given the top hit and the sequence, estimate what the
# maximum possible bit score would be.
# (If the top hit is exact, this will just return the top bit score.)
sub estimateTopScore($$) {
  my ($top, $seq) = @_;
  my $maxBits = $top->{bits};
  my $maxCov = ($top->{qEnd} - $top->{qBegin} + 1) / length($seq);
  my $maxIdentity = $top->{identity};
  return $maxBits / ($maxCov * $maxIdentity);
}
1;

# cumulative sum
sub cumsum {
  my (@in) = @_;
  my @out = ();
  my $tot = 0;
  for (my $i = 0; $i < scalar(@in); $i++) {
    $tot += $in[$i];
    $out[$i] = $tot;
  }
  return @out;
}
