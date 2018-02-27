#!/usr/bin/perl
use strict;
use warnings;
use autodie;

open my $out, '>', 'read_counts.tsv';
while(<>){
  chomp;
  my @parts = split("\t");
  next if $parts[0] eq '"chromosome"';
  print $out join("\t", @parts ) . "\n";
}
