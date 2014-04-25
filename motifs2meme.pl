#!/usr/bin/env genome-perl

use strict;
use warnings;
use Data::Dumper;

my @motif_lines;

print "MEME version 4\n\n";
print "ALPHABET= ACGT\n\n";

while(my $line = <>) {
    chomp $line;
    if ($line =~ />/) {
        print_motif(@motif_lines);
        @motif_lines = ();
    }
    push @motif_lines, $line;
}
print_motif(@motif_lines);

sub print_motif {
    my $header = shift;
    my @motif_lines = @_;

    unless (@motif_lines) {
        return;
    }
    my $name = substr($header, 1);
    print "MOTIF $name\n";
    print "letter-probability matrix: alength= 4 w= ".scalar @motif_lines."\n";
    for my $line (@motif_lines) {
        print substr($line, 2)."\n";
    }
    print "\n";
}
