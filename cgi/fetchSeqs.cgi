#!/usr/bin/perl -w
# Given a file with locus tags or protein identifiers, fetch their sequences
# Required CGI arguments: none
use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use DBI;
use HTML::Entities;
use lib "../lib";
use lib "../../PaperBLAST/lib";
use neighborWeb;

# Optional CGI arguments:
# order (which subdb to use)
# idfile -- uploaded file with identifiers

my $cgi = CGI->new;
setOrder(param('order'));
my $idFile = $cgi->upload('idfile');

unless ($idFile) {
  my $title = 'Fetch protein sequences';
  if (getOrder() eq "") {
    $title .= " from representative genomes";
  } else {
    $title .= " from " . encode_entities(getOrder());
  }
  start_page(title => $title);
  print p("Please upload a file with one locus tag (like A0F92_RS00005) or one protein identifier (like AAR38856.1) per line.",
          "Only the first field in each line will be considered.",
          "Identifiers that do not match will be ignored.",
          "Protein sequences will be returned in FASTA format.",
          "If you get fewer sequences than expected, please check that you are in the correct database",
          "(i.e., representative genomes versus a specific order such as Bacteroidales).");
  print
    start_form( -name => 'input', -method => 'POST', -action => 'fetchSeqs.cgi' ),
    orderToHidden(),
    p('Upload:', filefield(-name => 'idfile', -size => 50),
      submit('Go'));
  finish_page();
  exit(0);
}

my $fh = $idFile->handle;
my @ids = ();
while (my $line = <$fh>) {
  $line =~ s/\s.*//;
  next unless $line =~ m/\d/; # not an identifier
  push @ids, $line;
}
close($fh) || die "Error reading uploaded file";

print "Content-Type:text\n";
my $fileName = "fetched" . (getOrder() eq "" ? "" : "_" . getOrder()) . ".faa";
print "Content-Disposition: attachment; filename=$fileName\n\n";

my %seen = ();
foreach my $id (@ids) {
  next if exists $seen{$id};
  $seen{$id} = 1;
  my $gene;
  $gene = getDbHandle()->selectrow_hashref(qq{SELECT * from Gene JOIN Protein USING (proteinId) JOIN Genome USING (gid)
                                                WHERE locusTag = ? },
                                           {}, $id)
    || getDbHandle()->selectrow_hashref(qq{SELECT * from Gene JOIN Protein USING (proteinId) JOIN Genome USING (gid)
                                                WHERE proteinId = ? LIMIT 1 },
                                        {}, $id);
  if (defined $gene && exists $gene->{sequence}) {
    print ">" . $id . " " . $gene->{gtdbSpecies} . " " . $gene->{strain} . "\n"
      . $gene->{sequence} . "\n";
  }
}
