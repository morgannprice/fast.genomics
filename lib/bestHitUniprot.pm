#!/usr/bin/perl -w
package bestHitUniprot;
# Find the best hit in uniprot, if there is a close hit, using
# the SANSParallel web service
# (http://ekhidna2.biocenter.helsinki.fi/cgi-bin/sans/sans.cgi)

require Exporter;
use strict;
use LWP::Simple qw{get};
use HTML::Entities qw{encode_entities};
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
  my @hits = ();

  # The parsing code below also supports the AlphaFold database (consider using that in case their
  # copy of uniprot is broken)

  foreach my $db ("uniprot") {
    my $URL = "${base}?db=${db}&mode=table&H=1&seq=$seq";
    print "<!-- SANSparallel URL $URL -->\n";
    my $result = get($URL);
    return $parseError unless $result;
    my @lines = split /\n/, $result;
    my $iHeader;
    for (my $i = 0; $i < scalar(@lines); $i++) {
      if ($lines[$i] =~ m!^<TR><TH>Rank.*</TR>!i) {
        $iHeader = $i;
        last;
      }
    }
    print "<P>No header\n" unless defined $iHeader;
    return $parseError unless defined $iHeader;
    my $row = $lines[$iHeader+1];
    next unless $row =~ m!<TR><TD>1</TD><TD>(.*)</TD>!; # no hits

    # was getting geneName and species as well from the split
    my (undef, undef, $identity, $ranges, $alnLength, $bits, $eValue, $idSpec, $desc,
        $species, $geneName) = split "</TD>[ \t]*<TD[^>]*>", $row;
    return $parseError
      unless $ranges =~ m/^(\d+)-(\d+):(\d+)-(\d+)$/;
    my ($qBegin, $qEnd, $sBegin, $sEnd) = ($1,$2,$3,$4);
    my ($prefix, $uniprotId);
    if ($db eq "swiss" || $db eq "uniprot") {
      $idSpec = $1 if $idSpec =~ m!^<A HREF=.*>(.*)</A>$!i;
      return $parseError unless $idSpec =~ m/^([a-z]+)[|]([^|]+)/;
      ($prefix, $uniprotId) = ($1,$2);
      $species =~ s/ *$//;
    } else { # AlphaFold
      $desc =~ m/AFDB:AF-([^-]+)-[A-Z][0-9]+ (.*)/ || return $parseError;
      $uniprotId = $1;
      $desc = $2;
      $prefix = "tr"; # a hack
      $geneName = $idSpec;
      $species = "";
    }
    my $hit = { 'uniprotId' => $uniprotId, 'geneName' => $geneName, 'species' => $species,
                'prefix' => $prefix, desc => $desc,
                'identity' => $identity, 'eValue' => $eValue, 'bits' => $bits,
                'qBegin' => $qBegin, 'qEnd' => $qEnd, 'sBegin' => $sBegin, 'sEnd' => $sEnd,
                'service' => 'SANSparallel', 'serviceURL' => $base };
    return $hit if $identity >= 0.99;
    push @hits, $hit;
  } # end loop over db
  return { 'service' => 'SANSparallel', 'serviceURL' => $base}
    if @hits == 0;
  @hits = sort { $b->{identity} <=> $a->{identity} } @hits;
  return $hits[0];
}
