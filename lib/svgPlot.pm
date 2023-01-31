#!/usr/bin/perl -w
# Simple SVG scatterplots
package svgPlot;
require Exporter;
use strict;
use List::Util qw{min max};
use HTML::Entities;

# The input hash must include xRange and yRange (both reference to lists with 2 values)
# Options:
#   width, height,
#   xLeft, yTop,
#   margin # bottom left right bottom (lines)
# Returns an svgPlot object
sub new {
  my $class = shift;
  my %options = @_;
  die "Must provide xlim and ylim"
    unless exists $options{xRange}[1] && exists $options{yRange}[1];
  my $xRange = $options{xRange};
  my $yRange = $options{yRange};
  my $width = $options{width} || 400;
  my $height = $options{height} || 400;
  my $margin = $options{margin} || [3,3,2,1];
  my $xLeft = $options{xLeft} || 0;
  my $xRight = $xLeft + $width;
  my $yTop = $options{yTop} || 0;
  my $yBottom = $yTop + $height;
  my $lineSize = $options{lineSize} || 18;

  my $drawBottom = $yBottom - $margin->[0] * $lineSize;
  my $drawLeft = $xLeft + $margin->[1] * $lineSize;
  my $drawTop = $yTop + $margin->[2] * $lineSize;
  my $drawRight = $xRight - $margin->[3] * $lineSize;
  die "Not enough horizontal space for plot" unless $drawLeft < $drawRight;
  die "Not enough vertical space for plot" unless $drawBottom > $drawTop;

  return bless
    {
     xRange => $xRange,
     yRange => $yRange,
     drawLeft => $drawLeft,
     drawRight => $drawRight,
     drawTop => $drawTop,
     drawBottom => $drawBottom,
     xRight => $xRight,
     yBottom => $yBottom,
     lineSize => $lineSize,
     margin => $margin
    }, $class;
}

sub start {
  my ($self, $style) = @_;
  $style = "" if !defined $style;
  return qq{<SVG width="$self->{xRight}" height="$self->{yBottom}" style="$style"}
    . qq{<g transform="scale(1.0)">};
}

sub end($) {
  return "</g></svg>";
}

# $svgPlot->marginText(text, side, options)
# where side is one of bottom left top right.
# option may include
# lineAt -- defaults to max(objects' margin - 1.25, 0)
# style, title (hovertext), URL
# Returns a string
sub marginText {
  my $self = shift;
  my $text = shift;
  my $side = shift || die "no side given in svgPlot::marginText()";
  my (%options) = @_;
  my $lineAt = $options{lineAt};

  my $xMid = ($self->{drawLeft} + $self->{drawRight})/2;
  my $yMid = ($self->{drawTop} + $self->{drawBottom})/2;

  my ($x, $y);
  if ($side eq "bottom") {
    $x = $xMid;
    $lineAt = max(0, $self->{margin}[0] - 1.25) if !defined $lineAt;
    $y = $self->{drawBottom} + $self->{lineSize} * ($lineAt+1);
  } elsif ($side eq "top") {
    $x = $xMid;
    $lineAt = max(0, $self->{margin}[2] - 1.25) if !defined $lineAt;
    $y = $self->{drawTop} - $self->{lineSize} * $lineAt;
  } elsif ($side eq "left") {
    $lineAt = max(0, $self->{margin}[1] - 1.0) if !defined $lineAt;
    $x = $self->{drawLeft} - $self->{lineSize} * $lineAt;
    $y = $yMid;
  } elsif ($side eq "right") {
    $lineAt = max(0, $self->{margin}[3] - 1.0) if !defined $lineAt;
    $x = $self->{drawRight} + $self->{lineSize} * ($lineAt+1);
    $y = $yMid;
  } else {
    die "Unknown side $side";
  }


  my $out = qq{<TEXT text-anchor="middle" dominant-baseline="bottom" x="0" y="0"};
  my $transform;
  if ($side eq "left" || $side eq "right") {
    $transform = "translate($x,$y) rotate(-90)";
  } else {
    $transform = "translate($x,$y)";
  }
  $out .= qq{ transform="$transform"};
  my $style = $options{style};
  $out .= qq{style="$style"} if defined $style && $style ne "";
  $out .= ">";

  my $title = $options{title};
  if (defined $title && $title ne "") {
    $title = encode_entities($title);
    $out .= qq{<TITLE>$title</TITLE>};
  }
  $out .= encode_entities("$text") . "</TEXT>";

  my $URL = $options{URL};
  if (defined $URL && $URL ne "") {
    $URL = encode_entities($URL);
    $out = qq{<A xlink:href="$URL">$out</A>};
  }
  return $out;
}

sub convertX($$) {
  my ($self, $x) = @_;
  my ($rangeMin, $rangeMax) = @{ $self->{xRange} };
  return $self->{drawLeft}
    + ($self->{drawRight} - $self->{drawLeft}) * ($x - $rangeMin) / ($rangeMax - $rangeMin);
}

sub convertY($$) {
  my ($self, $y) = @_;
  my ($rangeMin, $rangeMax) = @{ $self->{yRange} };
  return $self->{drawBottom}
    - ($self->{drawBottom} - $self->{drawTop}) * ($y - $rangeMin) / ($rangeMax - $rangeMin);
}

sub axes {
  my $self = shift;
  return qq{<line x1="$self->{drawLeft}" x2="$self->{drawRight}" }
    . qq{ y1="$self->{drawBottom}" y2="$self->{drawBottom}" stroke="black" />}
    . qq{<line x1="$self->{drawLeft}" x2="$self->{drawLeft}" }
    . qq{ y1="$self->{drawBottom}" y2="$self->{drawTop}" stroke="black" />};
}

# arguments are axis ("x" or "y"), a list of positions to mark ticks at, and options
# Options may include labels
sub axisTicks {
  my $self = shift;
  my $axis = shift;
  my $ats = shift;
  my %options = (@_);
  my $labels = $options{labels} || $ats;
  die "Unknown axis $axis" unless $axis eq "x" || $axis eq "y";
  my @out = ();
  foreach my $i (0..(scalar(@$ats)-1)) {
    my $at = $ats->[$i];
    my $label = encode_entities($labels->[$i]);
    if ($axis eq "x") {
      my $x = $self->convertX($at);
      my $tickBottom = $self->{drawBottom} + 0.5 * $self->{lineSize};
      my $textAt = $self->{drawBottom} + 1.5 * $self->{lineSize};
      push @out, qq{<line x1="$x" x2="$x" y1="$self->{drawBottom}" y2="$tickBottom" stroke="black" />};
      push @out, qq{<text text-anchor="middle" dominant-baseline="bottom" x="$x" y="$textAt">$label</text>} if $label ne "";
    } else {
      my $y = $self->convertY($at);
      my $tickLeft = $self->{drawLeft} - 0.5 * $self->{lineSize};
      my $textAt = $self->{drawLeft} - 0.6 * $self->{lineSize};
      push @out, qq{<line x1="$self->{drawLeft}" x2="$tickLeft" y1="$y" y2="$y" stroke="black" />};
      push @out, qq{<text dominant-baseline="bottom" text-anchor="middle" x="0" y="0" }
        . qq{ transform = "translate($textAt,$y) rotate(-90)">$label</text>} if $label ne "";
    }
  }
  return join("\n", @out);
}

# arguments are self, x, y, and optionally a hash with URL, title, color, and size
# (By default, points are black empty circles, and the size is the radius)
sub point {
  my $self = shift;
  my $xIn = shift;
  my $yIn = shift;
  my %options = (@_);
  my $x = $self->convertX($xIn);
  my $y = $self->convertY($yIn);
  my $size = $options{size} || 2;
  my $color = $options{color} || "black";
  my $fill = $options{fill} || "none";
  my $out = qq{<circle cx="$x" cy="$y" r="$size" stroke="$color" fill="$fill"};
  my $style = $options{style};
  $out .= qq{ style="$style"} if defined $style && $style ne "";
  $out .= ">";
  my $title = $options{title};
  if (defined $title && $title ne "") {
    $out .= "<TITLE>" . encode_entities($title) . "</TITLE>";
  }
  $out .= "</circle>";
  my $URL = $options{URL};
  if (defined $URL && $URL ne "") {
    $URL = encode_entities($URL);
    $out = qq{<A xlink:href="$URL">$out</A>};
  }
  return $out;
}


