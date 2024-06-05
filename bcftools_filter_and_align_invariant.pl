#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "
bcftools_filter_and_align_invariant.pl <reference.fasta> <table.txt>

Given the reference fasta sequence file and the 'table.txt' file
output by bcftools_filter_and_align.pl will output a breakdown
of invariant site counts by base in the format 'A,C,G,T'. This
string can be used as input to iqtree using the '-fconst' option.

";

die $usage unless scalar @ARGV >= 2;

my ($ref, $table) = @ARGV;

open (my $tin, "<$table") or die "ERROR: Can't open $table: $!\n";
my %skip;
my $preskip = 0;
while (my $line = <$tin>){
    chomp $line;
    next if $line =~ m/^id\spos/;
    my @tmp = split("\t", $line);
    $skip{$tmp[0]}{$tmp[1]} = 1;
    $preskip++;
}
close ($tin);

my @counts = (0) x 4;
my %hash = (A => 0,
            C => 1,
            G => 2,
            T => 3);
open (my $rin, "<$ref") or die "ERROR: Can't open $ref: $!\n";
my ($id, $seq);
my $skipcount = 0;
while (my $line = <$rin>){
    chomp $line;
    if ($line =~ m/^>/){
        if ($id){
            my $pos = 0;
            $seq = uc($seq);
            foreach my $base (split(//, $seq)){
                $pos++;
                if ($skip{$id}{$pos}){
                    $skipcount++;
                    next;
                }
                next unless $base =~ m/^[ACGT]$/;
                $counts[$hash{$base}]++;
            }
        }
        $seq = "";
        $id = substr($line, 1);
        $id =~ s/\s+.*//;
        next;
    }
    $line =~ s/\s//g;
    $seq .= $line;
}
close ($rin);
if ($id){
    my $pos = 0;
    $seq = uc($seq);
    foreach my $base (split(//, $seq)){
        $pos++;
        if ($skip{$id}{$pos}){
            $skipcount++;
            next;
        }
        next unless $base =~ m/^[ACGT]$/;
        $counts[$hash{$base}]++;
    }
}
print "Found $skipcount of $preskip snp sites in the reference.\n\n";

print join(",", @counts), "\tA,C,G,T\n\n";
