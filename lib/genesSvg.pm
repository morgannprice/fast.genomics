# Draw genes on a track, in svg format
use strict;
package genesSvg;
require Exporter;
use List::Util qw{min max};
use HTML::Entities;

our (@ISA,@EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw{genesSvg scaleBarSvg};

my $defaultKbWidth = 150;

# Helper routine for position on chromosome to SVG x coordinate (not clipped)
sub posToX($$$$$$) {
  my ($pos, $xLeft, $begin, $end, $ntWidth, $invert) = @_;
  return $xLeft + ($end-$pos) * $ntWidth if $invert;
  #else
  return $xLeft + ($pos-$begin) * $ntWidth;
}

# genesSvg() produces the SVG for a track (not an entire SVG object).
#
# The first argument is a list of genes to show
#	each gene should include begin, end, strand, and optionally
#	color, label, desc, URL, and bar, which if present should contain
#	beginFraction, endFraction, and optionally color and title
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
#	scaffoldLength -- the length of the scaffold (affects extent of center line)
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
  my $kbWidth = $param{kbWidth} || $defaultKbWidth;
  my $ntWidth = $kbWidth/1000;
  my $xLeft = $param{xLeft} || 0;
  my $xRight = $xLeft + ($end-$begin) * $ntWidth;
  my $invert = $param{invert} || 0;
  my $arrowWidth = $param{arrowWidth} || 15;
  my $geneHeight = $param{geneHeight} || 18;
  my $showLabel = $param{showLabel} || 0;

  # ensure genes are sorted
  my @genes = sort { $a->{begin} <=> $b->{begin} } @$genes;

  my @svgLines = ();
  my $geneYMid = $yTop + $geneHeight/2;
  my ($lineLeft, $lineRight); # if invert, the direction might be swapped
  if (defined $param{scaffoldLength}) {
    my $scaffoldLength = $param{scaffoldLength};
    $lineLeft = max($xLeft, min($xRight, posToX(1, $xLeft, $begin, $end, $ntWidth, $invert)));
    $lineRight = max($xLeft, min($xRight, posToX($scaffoldLength, $xLeft, $begin, $end, $ntWidth, $invert)));
  } else {
    $lineLeft = $xLeft;
    $lineRight = $xRight;
  }
  push @svgLines,
    qq[<line x1="$lineLeft" y1="$geneYMid" x2="$lineRight" y2="$geneYMid" style="stroke:darkgrey; stroke-width:1;"/>];
  my $showPlus = $invert ? "-" : "+";
  foreach my $gene (@genes) {
    my $x1 = posToX($gene->{begin}, $xLeft, $begin, $end, $ntWidth, $invert);
    my $x2 = posToX($gene->{end}, $xLeft, $begin, $end, $ntWidth, $invert);

    my $truncated = 0;
    if ($x1 < $xLeft) {
      $x1 = $xLeft;
      $truncated = 1;
    }
    if ($x2 < $xLeft) {
      $x2 = $xLeft;
      $truncated = 1;
    }
    if ($x1 > $xRight) {
      $x1 = $xRight;
      $truncated = 1;
    }
    if ($x2 > $xRight) {
      $x2 = $xRight;
      $truncated = 1;
    }

    my $y1 = $yTop;
    my $y2 = $yTop + $geneHeight; # SVG has +y axis going down
    my ($xStart, $xStop) = ($x1,$x2);
    ($xStart,$xStop) = ($x2,$x1) if $gene->{strand} eq "-";
    my @points; # list of x,y pairs
    if (abs($xStart-$xStop) < $arrowWidth) {
      @points = ([$xStart,$y2], [$xStart,$y1], [$xStop,$geneYMid]);
    } else {
      my $xMid = $xStart < $xStop ? $xStop - $arrowWidth : $xStop + $arrowWidth;
      @points = ([$xStart,$y2], [$xStart,$y1], [$xMid,$y1],
                 [$xStop,$geneYMid], [$xMid, $y2]);
    }
    my $pointstr = join(" ", map { $_->[0].",".$_->[1] } @points);
    my $color = $gene->{color} || "white";
    my $poly = qq{<polygon points="$pointstr" style="fill:$color; stroke:black; stroke-width:1;" />};
    if ($color =~ m/url[(]/) {
      # hatch shading -- put a background behind it
      $poly = qq{<polygon points="$pointstr" style="fill:grey;" />}
        . $poly;
    }
    my $URL = encode_entities( $gene->{URL} || "");
    my $title = encode_entities( join(": ", $gene->{label}, $gene->{desc}) );
    $title = "(extends beyond this view) $title" if $truncated;
    push @svgLines, qq{<a xlink:href="$URL">},
      "<title>$title</title>",
      $poly;
    if ($showLabel && defined $gene->{label} && $gene->{label} ne "") {
      my $tag = $gene->{label};
      $tag =~ s/^.*_/_/; # shorten so that it is likely to fit
      $tag =~ s/[<>'"]//g;
      if (abs($xStop-$xStart) >= length($tag) * 9) {
        my $xMid = ($xStart+$xStop)/2;
        push @svgLines,
          qq[<text x="$xMid" y="$geneYMid" text-anchor="middle" dominant-baseline="middle" font-size="smaller" >],
          qq[$tag</text>];
      }
    }
    push @svgLines, "</a>";
    if (exists $gene->{bar} && ! $truncated) {
      my $bar = $gene->{bar};
      my ($xBar1, $xBar2);
      $xBar1 = $xStart + $bar->{beginFraction} * ($xStop-$xStart);
      $xBar2 = $xStart + $bar->{endFraction} * ($xStop-$xStart);
      #my $yBar = $yTop + $geneHeight + 4;
      my $yBar = ($y1+$y2)/2;
      my $barColor = $bar->{color} || "darkgrey";
      my $line  = qq[<line x1="$xBar1" x2="$xBar2" y1="$yBar" y2="$yBar" style="stroke: $barColor; stroke-width: 3;">];
      $line .= "<title>" . encode_entities($bar->{title}) . "</title>"
        if defined $bar->{title} && $bar->{title} ne "";
      $line .= "</line>";
      if (exists $bar->{URL} && $bar->{URL} ne "") {
        my $URL = encode_entities($bar->{URL});
        $line = qq[<a xlink:href="$URL">$line</a>];
      }
      push @svgLines, $line;
      # re-draw the lines on top
      push @svgLines, qq{<polygon points="$pointstr" style="fill:none; stroke:black; stroke-width:1;" />};
    }
  }

  # add invisible boxes in between adjacent genes, with hover text for how far apart they are
  for (my $i = 0; $i < scalar(@$genes) - 1; $i++) {
    my $gene1 = $genes[$i];
    my $gene2 = $genes[$i+1];
    my $x1 = posToX($gene1->{end}, $xLeft, $begin, $end, $ntWidth, $invert);
    my $x2 = posToX($gene2->{begin}, $xLeft, $begin, $end, $ntWidth, $invert);
    ($x1, $x2) = ($x2, $x1) unless $x2 > $x1;
    next unless $x1 >= $xLeft && $x2 <= $xRight;
    # Ensure it is wide enough to hover on
    $x1 = max($xLeft, $x1 - 3);
    $x2 = min($xRight, $x2 + 3);
    my $distText;
    if ($invert) {
      $distText = "$gene2->{label} and $gene1->{label}";
    } else {
      $distText = "$gene1->{label} and $gene2->{label}";
    }
    $distText =~ s/[<>'"]//g;
    if ($gene1->{end} >= $gene2->{begin}) {
      $distText .= " overlap by " . ($gene1->{end} - $gene2->{begin} + 1 ) . " nucleotides";
    } else {
      $distText .= " are separated by " . ($gene2->{begin} - $gene1->{end} - 1 ) . " nucleotides";
    }
    my $dx = $x2 - $x1;
    my $y1 = $geneYMid - 3;
    my $dy = 6;
    push @svgLines,
      qq{<rect x="$x1" width="$dx" y="$y1" height="$dy" stroke="none" fill-opacity="0"><title>$distText</title></rect>};
  }
  return (svg => join("\n", @svgLines) . "\n",
          xMax => int(0.99 + $xRight), yMax => int(0.5 + $yTop + $geneHeight + 8));
}

# input should be a hash that includes yTop
# and optionally xLeft, kbWidth
# returns a hash with svg, xMax, and yMax
sub scaleBarSvg {
  my (%param) = @_;
  my $yTop = $param{yTop};
  die "yTop must be specified" unless defined $yTop;
  my $xLeft = $param{xLeft} || 0;
  my $kbWidth = $param{kbWidth} || $defaultKbWidth;

  my $barLeft = $xLeft + 1;
  my $barRight = $xLeft + $kbWidth;
  my $barHeight = 10;
  my $barY = $yTop + 20;
  my $barTop = $barY - $barHeight/2;
  my $barBottom = $barY + $barHeight/2;

  my @svgLines = ();
  push @svgLines, qq[<line x1="$barLeft" x2="$barRight" y1="$barY" y2="$barY" stroke="black" />];
  push @svgLines, qq[<line x1="$barLeft" x2="$barLeft" y1="$barTop" y2="$barBottom" stroke="black" />];
  push @svgLines, qq[<line x1="$barRight" x2="$barRight" y1="$barTop" y2="$barBottom" stroke="black" />];
  my $textLeft = $barRight + 5;
  push @svgLines, qq[<text x="$textLeft" y="$barBottom" text-anchor="left" dominant-baseline="bottom">1 kb</text>];

  return (svg => join("\n", @svgLines) . "\n",
          xMax => $textLeft + 4*20,
          yMax => $barBottom + 5);
}
