#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use feature qw/ say /;
use Data::Printer;

my $control = $ARGV[0];
my $tumour = $ARGV[1];

open my $control_counts, '<', $control;
open my $tumour_counts, '<', $tumour;

my %c_counts;

my @chroms = qw/2L 2R 3L 3R 4 X Y/;
my %chrom_filt;
@chrom_filt{@chroms} = ();

my $previous_pos = 0;
my $window_size;
while(<$control_counts>){
    chomp;
    my ($chrom, $start, $count) = split;
    next if not exists $chrom_filt{$chrom};

    if ($previous_pos>0){
      $window_size = $start - $previous_pos;
    }

    $previous_pos = $start;

    $c_counts{$chrom}{$start} = $count;
}

my %t_counts;

while(<$tumour_counts>){
    chomp;
    my ($chrom, $start, $count) = split;
    next if not exists $chrom_filt{$chrom};
    $t_counts{$chrom}{$start} = $count;
}

open my $out, '>', 'normalised_read_counts.txt';


print $out "chromosome\tstart\tend\ttest\tref\n";

for my $chrom (sort keys %chrom_filt){
    for my $pos (sort { $a <=> $b } keys %{$c_counts{$chrom}}){
        print $out join("\t", $chrom, $pos, $pos+$window_size, $t_counts{$chrom}{$pos}, $c_counts{$chrom}{$pos}) . "\n";
    }
}
