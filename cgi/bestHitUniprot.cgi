#!/usr/bin/perl -w
use strict;
use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use IO::Handle; # for autoflush
use HTML::Entities;
use URI::Escape;
use lib "../lib";
use bestHitUniprot;
# neighborWeb.pm relies on various PaperBLAST libraries
use lib "../../PaperBLAST/lib";
use neighborWeb;
use pbweb qw{analysisLinks UniProtToFasta};

# CGI arguments:
# query (optional) -- an identifier or locus tag, or a sequence in fasta format
#	see neighborWeb::parseGeneQuery for details.
my $cgi = CGI->new;
my $query = $cgi->param('query') || "";

start_page('title' => "Find the best match in UniProt",
           'banner' => i("fast.genomics"));
print p({-style => "font-size: smaller;"}, "Due to technical problems with the SANSparallel server, the best match is currently computed using smaller databases (SwissProt and AlphaFoldDB v2) instead of against all of UniProt. Many proteins may lack close hits.")."\n";

autoflush STDOUT 1; # show preliminary results

my %query = parseGeneQuery($query);

if (!defined $query{seq}) {
  print p(b($query{error})) if $query{error};
  print start_form( -name => 'input', -method => 'GET', -action => 'bestHitUniprot.cgi' ),
    p("Enter an identifier from UniProt, RefSeq, PDB, or MicrobesOnline,",
      br(),
      "or a protein sequence in FASTA or Uniprot format",
      br(),
      textarea( -name  => 'query', -value => '', -cols  => 80, -rows  => 4 ),
      br(),
      br(),
      submit('Search')),
    end_form;
  finish_page();
}

print p("Searching for", encode_entities($query{seqDesc}), "(" . length($query{seq}) . " amino acids)")
  if ($query{seqDesc});
print showSequence($query{seqDesc}, $query{seq}), "\n";
my $hit = bestHitUniprot($query{seq});
if (exists $hit->{error}) {
  print p("Error:", encode_entities($hit->{error}));
} elsif (exists $hit->{uniprotId}) {
  my $uniprotId = $hit->{uniprotId};
  my $uniprotURL = "https://www.uniprot.org/uniprot/" . $uniprotId;
  my $interproURL = "https://www.ebi.ac.uk/interpro/protein/UniProt/" . $uniprotId;
  my $identityShow = int($hit->{identity} * 100 + 0.5);
  my $descShow = encode_entities($hit->{desc});
  $descShow .= " (" . i(encode_entities($hit->{species})) . ")" if $hit->{species};
  print
    p("${identityShow}% identical to",
      $hit->{prefix} eq "sp" ? "Swiss-Prot" : "UniProt",
      a({-href => $uniprotURL}, $uniprotId),
      $hit->{geneName},
      "(".a({-href => $interproURL}, "InterPro").")"),
   "\n",
   p("Description:", $descShow),
   "\n";
  # Fetch the target's sequence and build link to alignment
  my $fasta = UniProtToFasta($uniprotId);
  my @lines = split /\n/, $fasta;
  my $header = shift @lines; $header =~ s/^>//;
  my $hitSeq = join("", @lines);
  my $genes = $query{genes};
  my $gene;
  $gene = $genes->[0] if defined $genes && scalar(@$genes) == 1;
  my $alnURL = "alignPair.cgi?" . geneSeqDescSeqOptions($gene, $query{seqDesc}, $query{seq})
    . "&seq2=$hitSeq&seqDesc2=" . uri_escape($header);
  print
    p($hit->{qBegin}.":".$hit->{qEnd}."/".length($query{seq}),
      "of query aligns to",
      $hit->{sBegin}.":".$hit->{sEnd}."/".length($hitSeq),
      "of",  encode_entities($uniprotId).",",
      "see",
      a({-href => $alnURL }, "alignment"));
} else {
  # no hits
  print p("Sorry, no close match to this sequence was found in UniProt");
}

print p(a({-href => "bestHitUniprot.cgi"}, "Try another query"));

print p("Powered by", a({-href => $hit->{serviceURL}}, $hit->{service}))
  if exists $hit->{service};

print
  h3("Other sequence analysis tools"),
  start_ul,
  join("\t", map li($_), analysisLinks('desc' => $query{seqDesc},
                                       'seq' => $query{seq},
                                       'fbLoad' => 1,
                                       'skip' => {"UniProt" => 1 })),
  end_ul;
finish_page();
