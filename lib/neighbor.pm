# Utilities for the gene neighbor viewer
package neighbor;
require Exporter;
use strict;

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw(parseTaxString readMMSeqsHits estimateTopScore cumsum featuresToGenes orderToSubName
             csvQuote lengthToMMSeqsSens);

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

# Given a reference to a list of features, from FetchAssembly::ParseNCBIFeatureFile,
# return a list of genes.
# In particular, combine gene and CDS entries into protein-coding genes.
# Returns a reference to a list of genes, each a hash with the fields
# scaffoldId, begin, end, strand, locusTag, proteinId, desc;
# and a reference to a list of warnings.
sub featuresToGenes($) {
  my ($features) = @_;

  # The feature file usually has pairs of rows, first a gene entry, then usually a CDS
  # (even if a pseudogene), or rRNA or ncRNA type entry, which will have the description under "name".
  # But sometimes the rows are not in order (i.e., two rows for an antisense RNA between
  # the gene and CDS entries, see BSU_32469 in GCF_000009045.1). Or sometimes
  # there's no CDS entry at all.
  # So first, group the entries by locus_tag.
  my %locusTagFeatures = ();
  foreach my $feature (@$features) {
    push @{ $locusTagFeatures{ $feature->{locus_tag} } }, $feature;
  }

  my @genes;
  my @warnings;

  foreach my $locus_tag (sort keys %locusTagFeatures) {
    if ($locus_tag eq "") {
      push @warnings, "Empty locus_tag";
      next;
    }
    my $list = $locusTagFeatures{$locus_tag};
    my @geneFeatures = grep $_->{"# feature"} eq "gene", @$list;
    my @others = grep $_->{"# feature"} ne "gene", @$list;
    if (@geneFeatures == 0) {
      push @warnings, "No gene entry for locus_tag $locus_tag";
      next;
    }
    push @warnings, "More than one gene entry for locus_tag $locus_tag"
      if @geneFeatures > 1;
    print STDERR "Warning: more than one non-gene entry for locus_tag $locus_tag\n"
      if @others > 1;
    my $feature = $geneFeatures[0];
    my $other = $others[0]; # or undef

    my $locusTag = $feature->{"locus_tag"};
    my $proteinId = "";
    my $desc = "?";
    if (defined $other) {
      $desc = $other->{name};
      $desc = $desc . " (pseudogene)" if $other->{class} eq "without_protein";
    } else {
      $desc = "pseudogene" if $feature->{attributes} =~ m/pseudo/i;
    }
    # The CDS may have class="with_protein" or empty
    if ($feature->{class} eq "protein_coding" && defined $other
        && $other->{"# feature"} eq "CDS"
        && ($other->{class} eq "with_protein" || $other->{class} eq "")) {
      $proteinId = $other->{"product_accession"};
    }
    push @genes, { 'scaffoldId' => $feature->{"genomic_accession"},
                   'start' => $feature->{start},
                   'end' => $feature->{end},
                   'strand' => $feature->{strand},
                   'locusTag' => $locusTag,
                   'proteinId' => $proteinId,
                   'desc' => $desc };
  }
  return (\@genes, \@warnings);
}

# Given the name of an order, return a similar string suitable for use as the directory name
# (As of Marc 2023, none of the orders in GTDB need changing.)
sub orderToSubName($) {
  my ($order) = @_;
  $order =~ s/[^A-Za-z0-9_.-]/../g;
  return($order);
}

# sqlite3 expects CVS format, not exactly tab delimited format
# So, need to replace any " with "" and surround the field with quotes.
sub csvQuote($) {
  my ($in) = @_;
  return $in unless $in =~ m/"/;
  $in =~ s/"/""/g;
  return '"' . $in . '"';
}

# Choose sensitivity for MMseqs2 given the number of amino acids in the query
sub lengthToMMSeqsSens($) {
  my ($len) = @_;
  return 7.5 if $len <= 100;
  return 7.0 if $len <= 150;
  return 6.0 if $len <= 800;
  return 5.7;
}

1;
