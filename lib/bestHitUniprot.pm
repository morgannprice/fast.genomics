#!/usr/bin/perl -w
package bestHitUniprot;
# Find the best hit in uniprot, if there is a close hit, using
# the SANSParallel web service
# (http://ekhidna2.biocenter.helsinki.fi/cgi-bin/sans/sans.cgi)

require Exporter;
use strict;
use LWP::Simple qw{get};
our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw{bestHitUniprot};

my $parseError =  {"error" => "cannot parse the response from SANSparallel" };

# Given a protein sequence, returns a hash with
# uniprotId, geneName, species, desc,
# prefix (tr or sp)
# identity (0-1), eValue, bits, qBegin, qEnd, sBegin, sEnd
# or a hash with error
# or an empty hash if there is no hit
sub bestHitUniprot($) {
  my ($seq) = @_;
  return { 'error' => 'no input' } unless defined $seq && $seq ne "";
  $seq =~ s/[*]//g; # occasionally occurs from translating stop codons
  $seq = uc($seq);
  die "Invalid input to bestHitUniprot" unless $seq =~ m/^[A-Z]+$/;
  # H  =1 means only the best hit
  my $base = "http://ekhidna2.biocenter.helsinki.fi/cgi-bin/sans/sans.cgi";
  my $URL = "${base}?db=uniprot&mode=table&H=1&seq=$seq";
  my $result = get($URL);
  return $parseError unless $result;
  my @lines = split /\n/, $result;
  my $iHeader;
  for (my $i = 0; $i < scalar(@lines); $i++) {
    if ($lines[$i] =~ m!^<TR><TH>Rank.*</TR>$!i) {
      $iHeader = $i;
      last;
    }
  }
  return { 'error' => 'cannot parse response from SANSparallel' } unless
    defined $iHeader;
  my $row = $lines[$iHeader+1];
  return { 'service' => 'SANSparallel', 'serviceURL' => $base}
    unless $row =~ m!<TR><TD>1</TD><TD>(.*)</TD></TR>!; # no hits
  my (undef, $identity, $ranges, $alnLength, $bits, $eValue, $idSpec, $desc, $species, $geneName) = split "</TD><TD>", $1;
  $species =~ s/ *$//;
  return $parseError unless defined $geneName;
  $idSpec = $1 if $idSpec =~ m!^<A HREF=.*>(.*)</A>$!i;
  return $parseError unless $idSpec =~ m/([a-z]+)[|]([^|]+)/;
  my ($prefix, $uniprotId) = ($1,$2);
  # parse ranges
  return $parseError
    unless $ranges =~ m/^(\d+)-(\d+):(\d+)-(\d+)$/;
  my ($qBegin, $qEnd, $sBegin, $sEnd) = ($1,$2,$3,$4);
  return { 'uniprotId' => $uniprotId, 'geneName' => $geneName, 'species' => $species,
           'prefix' => $prefix, desc => $desc,
           'identity' => $identity, 'eValue' => $eValue, 'bits' => $bits,
           'qBegin' => $qBegin, 'qEnd' => $qEnd, 'sBegin' => $sBegin, 'sEnd' => $sEnd,
           'service' => 'SANSparallel', 'serviceURL' => $base };
}
