# Draw genes on a track, in svg format
use strict;
package genesSvg;
require Exporter;
use List::Util qw{min max};
use HTML::Entities;

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw{genesSvg};

# genesSvg() produces the SVG for a track (not an entire SVG object).
#
# The first argument is a list of genes to show
#	each gene should include begin, end, strand, and optionally
#	color, label, desc, and URL
# The remaining arguments are parameters, which may include
#       begin/end -- tracks will be clipped to these coordinates
#	   defaults to minimum and maximum begin/end of genes
#	yTop -- top of the track (default 0)
#	kbWidth -- width of one kilobase (default 150)
#	xLeft -- left of the track (default 0)
#	invert -- set to 1 to flip the strand
#	arrowWidth -- width in SVG units of the triangular part of each gene
#	geneHeight -- how high each gene should be
#	showLabel -- whether or not to draw the gene labels
#
# Returns a hash which includes
#	svg (string), xMax, yMax
sub genesSvg {
  my $genes = shift;
  die unless defined $genes;
  my %param = @_;
  my $begin = $param{begin};
  $begin = min(map $_->{begin}, @$genes) if !defined $begin;
  my $end = $param{end};
  $end = max(map $_->{end}, @$genes) if !defined $end;
  my $yTop = $param{yTop} || 0;
  my $kbWidth = $param{kbWidth} || 200;
  my $ntWidth = $kbWidth/1000;
  my $xLeft = $param{xLeft} || 0;
  my $xRight = $xLeft + ($end-$begin) * $ntWidth;
  my $invert = $param{invert} || 0;
  my $arrowWidth = $param{arrowWidth} || 15;
  my $geneHeight = $param{geneHeight} || 18;
  my $showLabel = $param{showLabel} || 0;
  # XXX no scale bar yet

  my @svgLines = ();
  my $geneYMid = $yTop + $geneHeight/2;
  push @svgLines,
    qq[<line x1="$xLeft" y1="$geneYMid" x2="$xRight" y2="$geneYMid" style="stroke:black; stroke-width:1;"/>];
  foreach my $gene (@$genes) {
    my ($x1, $x2, $showStrand);
    if ($invert) {
      $x1 = $xLeft + ($end - $gene->{end}) * $ntWidth;
      $x2 = $xLeft + ($end - $gene->{begin}) * $ntWidth;
      $showStrand = $gene->{strand} eq "+" ? "-" : "+";
    } else {
      $x1 = $xLeft + ($gene->{begin} - $begin) * $ntWidth;
      $x2 = $xLeft + ($gene->{end} - $begin) * $ntWidth;
      $showStrand = $gene->{strand};
    }
    my $y1 = $yTop;
    my $y2 = $yTop + $geneHeight; # SVG has +y axis going down
    my $geneYMid = ($y1+$y2)/2;
    my ($xStart, $xStop) = ($x1,$x2);
    ($xStart,$xStop) = ($x2,$x1) if $showStrand eq "-";
    my @points; # list of x,y pairs
    if (abs($xStart-$xStop) < $arrowWidth) {
      @points = ([$xStart,$y2], [$xStart,$y1], [$xStop,$geneYMid]);
    } else {
      my $xMid = $xStart < $xStop ? $xStop - $arrowWidth : $xStop + $arrowWidth;
      @points = ([$xStart,$y2], [$xStart,$y1], [$xMid,$y1], [$xStop,$geneYMid], [$xMid, $y2]);
    }
    my $pointstr = join(" ", map { $_->[0].",".$_->[1] } @points);
    my $color = $gene->{color} || "white";
    my $poly = qq{<polygon points="$pointstr" style="fill:$color; stroke:black; stroke-width:1;" />};
    my $URL = encode_entities( $gene->{URL} || "");
    my $title = encode_entities( join(": ", $gene->{label}, $gene->{desc}) );
    push @svgLines, qq{<a xlink:href="$URL">},
      "<title>$title</title>",
      $poly;
    if ($showLabel && defined $gene->{label} && $gene->{label} ne "") {
      my $tag = $gene->{label};
      $tag =~ s/^.*_/_/; # shorten so that it is likely to fit
      $tag =~ s/[<>'"]//g;
      if (abs($xStop-$xStart) >= length($tag) * 9 && $gene->{begin} >= $begin && $gene->{end} <= $end) {
        my $xMid = ($xStart+$xStop)/2;
        push @svgLines,
          qq[<text x="$xMid" y="$geneYMid" text-anchor="middle" dominant-baseline="middle" font-size="smaller" >],
          qq[$tag</text>];
      }
    }
    push @svgLines, "</a>";
  }
  return (svg => join("\n", @svgLines) . "\n",
          xMax => int(0.5 + $xRight), yMax => int(0.5 + $yTop + $geneHeight));
}
