#!/usr/bin/perl -w
use strict;

use CGI qw(:standard Vars start_ul);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use lib "../lib";
use lib "../../PaperBLAST/lib";
# neighborWeb.pm relies on various PaperBLAST libraries
use neighborWeb;
use svgPlot;

start_page(title => "test");
my $plot = svgPlot->new(xRange => [0,100], yRange => [0,100],
                        margin => [3,3,2,2]);
my @ticks = (0,20,40,60,80,100);
my @tickLabels = map $_."%", @ticks;
print
  $plot->start(),
  $plot->axes(),
  $plot->axisTicks("x", \@ticks, labels => \@tickLabels),
  $plot->axisTicks("y", \@ticks, labels => \@tickLabels),
  $plot->marginText("x label", "bottom"),
  $plot->marginText("y label", "left", title => "Hover text"),
  $plot->point(30,100, title => "title", URL => "search.cgi", color => "red"),
  $plot->point(0,0),
  $plot->point(100,100);

finish_page();
