#!/usr/bin/perl -w
package binom;
# Compute log probabilities for the binomial and hypergeometric distributions
# (All values are natural logarithms, i.e. ln(p) not log2(p) or log10(p).)
require Exporter;
use strict;

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw{binomLogProb binomLogTail hyperLogProb hyperLogProbTail
             logFactorial logChoose combineLogP};

# cache log factorial values; start with just 0 => 0
my @logFactorial = (0);

sub logFactorial($) {
  my ($n) = @_;
  die $n if $n < 0;
  if (!defined $logFactorial[$n]) {
    for (my $i = length(@logFactorial); $i <= $n; $i++) {
      $logFactorial[$i] = log($i) + $logFactorial[$i-1];
    }
  }
  return $logFactorial[$n];
}

sub logChoose($$) {
  my ($n, $k) = @_;
  return logFactorial($n) - logFactorial($k) - logFactorial($n-$k);
}

# log(probability of getting k sucesses out of n trials with probability prob)
sub binomLogProb($$$) {
  my ($k, $n, $prob) = @_;
  die "k $k n $n" unless $k >= 0 && $k <= $n;
  die "prob $prob" unless $prob > 0 && $prob < 1;

  my $logChoose = logFactorial($n) - logFactorial($k) - logFactorial($n-$k);
  return $logChoose + $k * log($prob) + ($n-$k) * log(1-$prob);
}

sub combineLogP($$) {
  my ($logP, $logQ) = @_;
  # log(p + q) = log(p * (1 + q/p)) = log(p) + log(1 + exp(log(q)-log(p)))
  # Ensure logQ <= logP to prevent overflow (underflow is harmless)
  ($logP, $logQ) = ($logQ, $logP) if $logQ > $logP;
  return $logP + log(1 + exp($logQ - $logP));
}

# Probability of $k or larger
sub binomLogTail($$$) {
  my ($k, $n, $prob) = @_;
  my $logP = binomLogProb($k, $n, $prob);
  for (my $j = $k+1; $j <= $n; $j++) {
    my $logQ = binomLogProb($j, $n, $prob);
    last if $logQ < $logP - 100; # exp(-50) = 2e-22, tail will no longer matter
    $logP = combineLogP($logP, $logQ);
  }
  return $logP;
}

# log(probability of getting x white balls when drawing from an urn
# with m white balls and n black balls and drawing k times
sub hyperLogProb($$$$) {
  my ($x, $m, $n, $k) = @_;
  die unless $k >= 0 && $k <= $m + $n;
  return logChoose($m, $x) + logChoose($n, $k-$x) - logChoose($m+$n, $k);
}

# log(probability of getting >= x white balls when drawing from an urn
# with m white balls and n black balls and drawing k times
sub hyperLogProbTail($$$$) {
  my ($x, $m, $n, $k) = @_;
  die unless $k >= 0 && $k <= $m + $n;
  my $logP = hyperLogProb($x, $m, $n, $k);
  for (my $j = $x+1; $j <= $m && $j <= $k; $j++) {
    my $logQ = hyperLogProb($j, $m, $n, $k);
    last if $logQ < $logP - 100; # exp(-50) = 2e-22, tail will no longer matter
    $logP = combineLogP($logP, $logQ);
  }
  return $logP;
}

