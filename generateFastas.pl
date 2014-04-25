#!/usr/bin/env genome-perl

use strict;
use warnings;
use Genome::File::Vcf::Reader;
use Genome::File::Vcf::Genotype;
use Data::Dumper;
use Genome;

my $variant_file = $ARGV[0];
my $fasta_file = $ARGV[1];
my $sample = $ARGV[2];
my $output_germline = $ARGV[3];
my $output_somatic = $ARGV[4];

my $in = Genome::File::Vcf::Reader->new($variant_file);
my $germline_out = Genome::Sys->open_file_for_writing($output_germline);
my $somatic_out = Genome::Sys->open_file_for_writing($output_somatic);

my $FLANK_SIZE = 25;

my $counter = 0;

my $sample_index = $in->header->index_for_sample_name($sample);

while (my $variant = $in->next) {
    next if $variant->is_filtered;
    my $ft = $variant->sample_field($sample_index,"FT");
    next if(defined $ft and ($ft ne '.' and $ft ne 'PASS'));

    my $left_sequence = get_left_flanking_region($variant);

    my $right_sequence = get_right_flanking_region($variant);

    my $variant_sequence = get_variant_sequence($variant);

    print $germline_out ">$counter\n";
    print $germline_out $left_sequence.$variant->{reference_allele}.$right_sequence."\n";

    print $somatic_out ">$counter\n";
    print $somatic_out $left_sequence.$variant_sequence.$right_sequence."\n";
    $counter++;
}
$germline_out->close;
$somatic_out->close;

sub get_left_flanking_region {
    my $variant = shift;
    my $chr = $variant->{chrom};
    my $start = $variant->{position} - $FLANK_SIZE;
    if ($start < 0) {
        $start = 0;
    }
    my $end = $variant->{position} - 1;
    return extract_sequence($chr, $start, $end);
}

sub get_right_flanking_region {
    my $variant = shift;
    my $chr = $variant->{chrom};
    my $length = length($variant->{reference_allele});
    my $start = $variant->{position} + $length;
    if ($start < 0) {
        $start = 0;
    }
    my $end = $start + $FLANK_SIZE - 1;
    return extract_sequence($chr,$start,$end);
}

sub extract_sequence {
    my ($chr, $start, $end) = @_;
    my $sequence =  `samtools faidx $fasta_file $chr:$start-$end | grep -v \">\"`;
    chomp $sequence;
    return $sequence;
}

sub  get_variant_sequence {
    my $variant = shift;
    my @entry_alleles = $variant->alleles;
    my @alleles = map{$entry_alleles[$_]} $variant->genotype_for_sample($sample_index)->get_alleles;
    my $variant_sequence = $alleles[1];
    unless (defined $variant_sequence) {
        $variant_sequence = "";
    }
    return $variant_sequence;
}
