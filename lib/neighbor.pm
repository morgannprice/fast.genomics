# Utilities for the gene neighbor viewer
package neighbor;
require Exporter;
use strict;

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw!parseTaxString!;

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

1;

