#!/usr/bin/env perl

my $version = "0.4";

use strict;
use warnings;

my $usage = "
bcftools_filter_consensus_to_table.pl
version: $version

Generates a table of variant positions in the same format as the one produced by
bcftools_filter_and_align.pl as well as a thinned multiple genome alignment.

Required:
  -r    fasta sequence of the reference genome
  -f    file of consensus genome seuquences ouptut by bcftools_filter.pl from 
        alignments against the reference sequence
        File format:
        /path/to/consensus.fasta <tab> genome_name

Optional:
  -i    include reference sequence in the output
        (default: reference sequence will not be included)
  -p    output file prefix
        (default: 'output')

";

use Getopt::Std;
use vars qw( $opt_f $opt_r $opt_i $opt_p );
getopts('r:f:p:i');

die $usage unless ($opt_r and $opt_f);

my $reffile = $opt_r;
my $fof     = $opt_f;
my $pref    = $opt_p ? $opt_p : "output";

## Read in the reference sequence;
my @order;
my @refseq;
open (my $rin, "<$reffile") or die "ERROR: Can't open $reffile: $!\n";
my ($id, $seq);
while (my $line = <$rin>){
    chomp $line;
    if ($line =~ m/^>/){
        if ($id){
            $seq = uc($seq);
            $seq =~ s/\s//g;
            utf8::downgrade($seq);
            push @refseq, $seq;
        }
        $seq = "";
        $id = substr($line, 1);
        $id =~ s/\s.+//;
        push @order, $id;
        next;
    }
    $seq .= $line;
}
close($rin);
if ($id){
    $seq = uc($seq);
    $seq =~ s/\s//g;
    $seq = utf8::downgrade($seq);
    push @refseq, $seq;
} else {
    die "ERROR: No fasta-formatted sequence found in $reffile";
}

## Read in the consensus sequences and compare 
open(my $fin, "<$fof") or die "ERROR: Can't open $fof: $!\n";
my %variants;
my %pos_bases;
my @gen_order;
my $first_gen = 1;
@refseq = () unless ($opt_i);
while (my $fline = <$fin>){
    chomp $fline;
    next if $fline =~ m/^\s*#/;
    next if $fline =~ m/^\s*$/;
    my ($path, $gen) = split("\t", $fline);
    print STDERR "Reading $gen\n";
    push @gen_order, $gen;
    open (my $in, "<$path") or die "ERROR: Can't open $path: $!\n";
    my $index = 0;
    ($id, $seq) = ("") x 2;
    while (my $line = <$in>){
        chomp $line;
        if ($line =~ m/^>/){
            if ($id){
                $seq = uc($seq);
                $seq =~ s/\s//g;
                utf8::downgrade($seq);
                if ($first_gen and !$opt_i){
                    push @refseq, $seq;
                } else {
                    my $rseq = $refseq[$index];
                    die "ERROR: lengths of sequence record $id and reference sequence record $order[$index] do not match\n" unless (length($seq) == length($rseq));
                    # xor the two strings
                    my $mask = $rseq ^ $seq;
                    while ($mask =~ m/[^\0]/g){
                        my $pos = $-[0];
                        my $qbase = substr($seq, $pos, 1);
                        my $rbase = substr($rseq, $pos, 1);
                        $variants{$index}{$pos}{$gen} = $qbase;
                        $pos_bases{$index}{$pos}{$qbase} = 1;
                        $pos_bases{$index}{$pos}{$rbase} = 1;
                    }
                }
                $index++;
            }
            $seq="";
            $id = substr($line, 1);
            next;
        }
        $seq .= $line;
    }
    close ($in);
    if ($id){
        $seq = uc($seq);
        $seq =~ s/\s//g;
        $seq = utf8::downgrade($seq);
        if ($first_gen and !$opt_i){
            push @refseq, $seq;
        } else {
            my $rseq = $refseq[$index];
            die "ERROR: lengths of sequence record $id and reference sequence record $order[$index] do not match\n" unless (length($seq) == length($rseq));
            # xor the two strings
            my $mask = $rseq ^ $seq;
            while ($mask =~ m/[^\0]/g){
                my $pos = $-[0];
                my $qbase = substr($seq, $pos, 1);
                my $rbase = substr($rseq, $pos, 1);
                $variants{$index}{$pos}{$gen} = $qbase;
                $pos_bases{$index}{$pos}{$qbase} = 1;
                $pos_bases{$index}{$pos}{$rbase} = 1;
            }
        }
    } else {
        die "ERROR: No No fasta-formatted sequence found in $path\n";
    }
    $first_gen = 0;
}

## Output the table and consensus sequence
my @align;
print STDERR "Outputting table...\n";
my $varcount = 0;
print STDERR "\r$varcount variant positions";
open (my $tout, ">$pref.table.txt");
print $tout "id\tpos";
print $tout "\tref" if $opt_i;
foreach my $gen (@gen_order){
    print $tout "\t$gen";
}
print $tout "\n";
for my $i (0 .. $#order){
    my $chrom = $order[$i];
    my $rseq = $refseq[$i];
    foreach my $pos_i (sort {$a <=> $b} keys %{$variants{$i}}){
        my @bases = keys %{$pos_bases{$i}{$pos_i}};
        my $is_gap;
        foreach (@bases){
            $is_gap = 1 if $_ =~ m/-/;
        }
        next if scalar @bases < 2 or (scalar @bases == 2 and $is_gap);
        my $rbase = substr($rseq, $pos_i, 1);
        my $pos = $pos_i + 1;
        print $tout "$chrom\t$pos";
        if ($opt_i){
            print $tout "\t$rbase,?,?,?";
            $align[0] .= $rbase; 
        }
        for my $j (0 .. $#gen_order){
            my $gen = $gen_order[$j];
            my $val = $rbase;
            if (exists $variants{$i}{$pos_i}{$gen}){
                $val = $variants{$i}{$pos_i}{$gen};
            }
            print $tout "\t$val,?,?,?";
            my $aj = $j;
            $aj += 1 if $opt_i;
            $align[$aj] .= $val;
        }
        $varcount++;
        print STDERR "\r$varcount variant positions" if $varcount % 100 == 0;
        print $tout "\n";
    }
}
close ($tout);
print STDERR "\r$varcount variant positions\n";


## Output the alignment
print STDERR "Outputting alignment...\n";
open (my $aout, ">$pref.alignment.fasta");
for my $i (0 .. $#align){
    my $seq = $align[$i];
    if ($i == 0){
        if ($opt_i){
            print $aout ">ref\n$seq\n";
        } else {
            my $gen = $gen_order[$i];
            print $aout ">$gen\n$seq\n";
        }
        next;
    }
    my $ai = $i;
    $ai-- if $opt_i;
    my $gen = $gen_order[$ai];
    print $aout ">$gen\n$seq\n";
}
close ($aout);
