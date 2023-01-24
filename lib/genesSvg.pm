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
#
# Returns a hash which includes
#	svg (string), xMax, yMax
sub genesSvg {
  my $genes = shift;
  die unless defined $genes;
  my %param = @_;
  my $begin = $param{begin} || min(map $_->{begin}, @$genes);
  my $end = $param{end} || max(map $_->{begin}, @$genes);
  my $yTop = $param{yTop} || 0;
  my $kbWidth = $param{kbWidth} || 200;
  my $ntWidth = $kbWidth/1000;
  my $xLeft = $param{xLeft} || 0;
  my $xRight = $xLeft + ($end-$begin) * $ntWidth;
  my $invert = $param{invert} || 0; # XXX not implemented yet
  my $arrowWidth = $param{arrowWidth} || 15;
  my $geneHeight = $param{geneHeight} || 18;
  # XXX no scale bar yet

  my @svgLines = ();
  my $geneYMid = $yTop + $geneHeight/2;
  push @svgLines,
    qq[<line x1="$xLeft" y1="$geneYMid" x2="$xRight" y2="$geneYMid" style="stroke:black; stroke-width:1;"/>];
  foreach my $gene (@$genes) {
    my $x1 = $xLeft + ($gene->{begin} - $begin) * $ntWidth;
    my $x2 = $xLeft + ($gene->{end} - $begin) * $ntWidth;
    my $y1 = $yTop;
    my $y2 = $yTop + $geneHeight; # SVG has +y axis going down
    my $geneYMid = ($y1+$y2)/2;
    my ($xStart, $xStop) = ($x1,$x2);
    ($xStart,$xStop) = ($x2,$x1) if $gene->{strand} eq "-";
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
      $poly, "</a>";
  }
  return (svg => join("\n", @svgLines) . "\n",
          xMax => int(0.5 + $xRight), yMax => int(0.5 + $yTop + $geneHeight));
}
