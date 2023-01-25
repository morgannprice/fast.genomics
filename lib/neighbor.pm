# Utilities for the gene neighbor viewer
package neighbor;
require Exporter;
use strict;

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw!parseTaxString parseLast readMMSeqsHits!;

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

# Given a file handle, returns a list of alignments, each containing
# query, subject, qBegin, qEnd, qLength, sBegin, sEnd, sLength, bits, score, eValue, match, (fraction) identity
# (1-based coordinates)
# Only protein alignments are supported.
sub parseLast($) {
  my ($fh) = @_;
  my $scoreLine = undef;
  my @alnLines = ();
  my @out = ();
  my ($lambda, $K);
  while (my $line = <$fh>) {
    chomp $line;
    if ($line =~ m/^# lambda=([0-9.+-]+) K=([0-9.+-]+)$/) {
      $lambda = $1;
      $K = $2;
    } elsif ($line =~ m/^#/) {
      next;
    } elsif ($line =~ m/^a score=/) {
      die "Unexpected score line $line in MAF file" if defined $scoreLine;
      $scoreLine = $line;
    } elsif ($line =~ m/^s /) {
      die "Unexpected alignment line $line in MAF file" unless defined $scoreLine;
      push @alnLines, $line;
    } elsif ($line eq "") {
      if (defined $scoreLine) {
        die "No alignment lines for $scoreLine in MAF file" unless @alnLines > 0;
        die "Wrong number of alignment lines for $scoreLine in MAF file" unless @alnLines == 2;
        $scoreLine =~ m/^a score=([0-9.]+) EG2=([0-9e.+-]+) E=([0-9.e+-]+)$/
          || die "Cannot parse score line $scoreLine in MAF file";
        my ($score, $EG2, $evalue) = ($1, $2, $3);
        my $bits = ($lambda * $score - log($K))/log(2);
        # Another way to compute bits, except that high-scoring alignments have $EG2 = 0
        # $bits = log(1e18/$EG2)/log(2) if $EG2 > 0;
        my @alnParsed = (); # id, length, begin, end, alignment string
        foreach my $alnLine (@alnLines) {
          my @pieces = split / +/, $alnLine;
          my ($sString, $id, $begin, $alnLength, $strand, $seqLength, $alnString) = @pieces;
          die "Cannot parse aligment line $alnLine in MAF file"
            unless @pieces == 7 && $sString eq "s"
              && $begin =~ m/^\d+$/ && $alnLength =~ m/^\d+$/ && $seqLength =~ m/^\d+$/
              && $strand eq "+"
              && $alnString =~ m/^[A-Za-z-]+$/;
          $begin++; # 1-based
          my $end = $begin + $alnLength - 1;
          die if $end > $seqLength;
          push @alnParsed, [$id, $seqLength, $begin, $end, $alnString];
        }
        my ($query, $qLength, $qBegin, $qEnd, $qAln) = @{ $alnParsed[0] };
        my ($subject, $sLength, $sBegin, $sEnd, $sAln) = @{ $alnParsed[1] };
        $qAln = uc($qAln);
        $sAln = uc($sAln);
        die "Alignment lengths do not match in MAF file for $scoreLine"
          unless length($qAln) == length($sAln);
        my $nMatch = 0;
        my $nGapOpen = 0;
        my $nMismatch = 0;
        my @qPos = split //, $qAln;
        my @sPos = split //, $sAln;
        for (my $i = 0; $i < length($qAln); $i++) {
          if ($qPos[$i] eq $sPos[$i]) {
            $nMatch++;
          } elsif ($qPos[$i] ne "-" && $sPos[$i] ne "-") {
            $nMismatch++;
          } elsif ($qPos[$i] eq "-" || $sPos[$i] eq "-") {
            $nGapOpen++ if $i == 0 || ($qPos[$i-1] ne "-" && $sPos[$i-1] ne "-");
          }
        }
        my $identity = 100 * $nMatch/length($qAln);
        my $out = { 'query' => $query, 'qLength' => $qLength, 'qBegin' => $qBegin, 'qEnd' => $qEnd,
                    'subject' => $subject, 'sLength' => $sLength, 'sBegin' => $sBegin, 'sEnd' => $sEnd,
                    'score' => $score, 'eValue' => $evalue, 'bits' => $bits, 'EG2' => $EG2,
                    'alnLength' => length($qAln),
                    'nMatch' => $nMatch, 'nMismatch' => $nMismatch, 'nGapOpen' => $nGapOpen,
                    'identity' => $identity };
        push @out, $out;
        $scoreLine = undef;
        @alnLines = ();
      }
    } else {
      die "Unexpected line $line in MAF file";
    }
  }
  die "Ended MAF file without finishing alignment" if defined $scoreLine;
  return \@out;
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
      unless $row{eValue} =~ m/^[0-9.E+-]+$/;
    push @out, \%row;
  }
  close($fh) || die "Error reading $file";
  return \@out;
}

1;

