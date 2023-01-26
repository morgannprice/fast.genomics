#!/usr/bin/perl -w
# Using LAST, find all similarities between a est of protein sequences
# and cluster them
# Assumes that the LAST executables lastal and lastdb are in $RealBin/../bin/
# uses $ENV{TMPDIR} or /tmp as the temporary directory for the last database and hits
package clusterProteins;
require Exporter;
use strict;
use Digest::MD5 qw{md5_hex};
use FindBin qw{$RealBin};

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(parseLast clusterProteins proteinsToSimilarity);

# Given a file handle, returns a list of alignments, each containing
# query, subject, qBegin, qEnd, qLength, sBegin, sEnd, sLength, bits, score, eValue,
# match, alnLen, and (fraction) identity.
# All coordinates are 1-based.
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

# Given a reference to a hash of proteinId to sequence,
# runs lastal, and returns a reference to
# a list of hits from neighbor::parseLast()
sub proteinsToSimilarity {
  my ($proteinSeq) = @_;
  return [] if scalar(keys %$proteinSeq) == 0;
  my $lastal = "$RealBin/../bin/lastal";
  my $lastdb = "$RealBin/../bin/lastdb";
  foreach my $x ($lastal, $lastdb) {
    die "No such executable: $x\n" unless -x $x;
  }

  my $md5 = md5_hex(join(",", sort keys %$proteinSeq));
  my $tmpDir = $ENV{TMPDIR} || "/tmp";
  my $tmpPre = "$tmpDir/$md5.$$";
  my $tmpFaa = "$tmpPre.faa";
  open(my $fh, ">", $tmpFaa) || die "Cannot write to $tmpFaa";
  foreach my $proteinId (sort keys %$proteinSeq) {
    print $fh ">$proteinId\n$proteinSeq->{$proteinId}\n";
  }
  close($fh) || die "Error writing $tmpFaa";
  my $tmpLast = "$tmpPre.last";
  my $cmd;
  $cmd = "$lastdb -p $tmpLast $tmpFaa";
  system($cmd) == 0 || die "$cmd -- failed: $!";
  $cmd = "$lastal -P 2 $tmpLast $tmpFaa > $tmpLast.out";
  system($cmd) == 0 || die "$cmd -- failed: $!";
  open(my $fhHits, "<", "$tmpLast.out") || die "Cannot read $tmpLast.out";
  my $hits = parseLast($fhHits);
  close($fhHits) || die "Error reading $tmpLast.out";
  unlink($tmpFaa);
  foreach my $suffix (qw{bck des prj sds ssp suf tis out}) {
    unlink("$tmpLast.$suffix");
  }
  return $hits;
}

# Given a reference to a hash of proteinId => protein sequence,
# and minimum coverage fraction such as 0.5,
# cluster proteins by similarity,
# and returns a list of clusters, each being a list of at least two protein ids.
# Uses proteinsToSimilarity() to find similar proteins.
sub clusterProteins {
  my ($proteinSeq, $minCoverage) = @_;
  my $hits = proteinsToSimilarity($proteinSeq);
  my @hits = grep $_->{query} ne $_->{subject}, @$hits;
  @hits = sort { $b->{score} <=> $a->{score} } @hits;
  my @clusters = (); # index to list of proteinIds
  my %cluster = (); # proteinId to index
  foreach my $hit (@hits) {
    my $query = $hit->{query};
    my $subject = $hit->{subject};
    next unless $hit->{qEnd} - $hit->{qBegin} + 1 >= $minCoverage * $hit->{qLength}
      || $hit->{sEnd} - $hit->{sBegin} + 1 >= $minCoverage * $hit->{sLength};
    if (exists $cluster{$query} && exists $cluster{$subject}) {
      # nothing ot do if already in the same cluster
      unless ($cluster{$query} eq $cluster{$subject}) {
        # move subjects' cluster into query's cluster
        my $oldI = $cluster{$subject};
        my $newI = $cluster{$query};
        my $oldCluster = $clusters[$oldI];
        push @{ $clusters[$newI] }, @$oldCluster;
        foreach my $proteinId (@$oldCluster) {
          $cluster{$proteinId} = $newI;
        }
        $clusters[$oldI] = [];
      }
    } elsif (exists $cluster{$query}) {
      # put subject in query's cluster
      my $iCluster = $cluster{$query};
      $cluster{$subject} = $iCluster;
      push @{ $clusters[$iCluster] }, $subject;
    } elsif (exists $cluster{$subject}) {
      # put query in subject's cluster
      my $iCluster = $cluster{$subject};
      $cluster{$query} = $iCluster;
      push @{ $clusters[$iCluster] }, $query;
    } else {
      # put both in a new cluster
      my $iCluster = scalar(@clusters);
      $cluster{$query} = $iCluster;
      $cluster{$subject} = $iCluster;
      push @clusters, [ $query, $subject ];
    }
  }
  @clusters = grep scalar(@$_) > 0, @clusters;
  return \@clusters;
}
