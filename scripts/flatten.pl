#!/usr/bin/perl
use strict;
use warnings;
use autodie;

my %seen;
open my $out, '>', 'flat_vars.tsv';
while(<>){
  chomp;
  my @parts = split("\t");
  next if $seen{$parts[0]}{$parts[2]}{$parts[3]}++;
  print $out join("\t", @parts ) . "\n";
}
