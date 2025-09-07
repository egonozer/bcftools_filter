#!/usr/bin/env perl

my $version = "0.4";

use strict;
use warnings;

my $usage = "

bcftools_filter_and_align_filter.pl [options]
version: $version

Filters a table produced by bcftools_filter_and_align.pl to remove
non-core sites and generate a new table and new fasta alignment.

Required:
  -t    table.txt file output by bcftools_filter_and_align.pl

Options:
  -m    minimum fraction of genomes in which a SNP position must be found in
        order to be counted (i.e. \"1\" for SNP position found in all genomes,
        \"0.7\" for SNP position found in at least 70% of the genomes)
        (default: any SNP position found in at least 2 of the input genomes
        will be included, i.e. \"0\")
  -i    file with list of genomes to include in the ouput. List should be genome
        names, one per line.
        Fraction limits given by -m will be based on these genomes (i.e. if
        -m is 0.5 and you enter list of 4 genomes for -i, output will be loci
        present in at least 2 of these 4 genomes)
        (default: all genomes will be included)
  -o    output file prefix. Will produce files <prefix>.table_filtered.txt and
        <prefix>.alignment_filtered.fasta
        (default: 'output')

";

use Getopt::Std;
use vars qw( $opt_t $opt_m $opt_i $opt_o );
getopts('t:m:i:o:');

die $usage unless ($opt_t);

my $intable = $opt_t;
my $minfrac = $opt_m ? defined $opt_m : 0;
my $incfile = $opt_i if $opt_i;
my $pref    = $opt_o ? $opt_o : "output";

my %incl;
my $inclcount;
if ($incfile){
    open (my $ifile, "<$incfile") or die "ERROR: Can't open $incfile: $!\n";
    while (my $line = <$ifile>){
        chomp $line;
        $line =~ s/^\s*//;
        $line =~ s/\s*$//;
        $incl{$line} = 1 if $line;
    }
    close ($ifile);
    $inclcount = scalar(keys(%incl));
}

open (my $tin, "<$intable") or die "ERROR: Can't open $intable: $!\n";
open (my $tout, ">$pref.table_filtered.txt");
my @seq;
my @header;
my @keepcols;
my ($incount, $outcount) = (0) x 2;
my $gencount;
while (my $line = <$tin>){
    chomp $line;
    next if $line =~ m/^\s*$/;
    my @tmp = split("\t", $line);
    unless (@header){
        if ($incfile){
            @header = @tmp[0..1];
            my $hcount = 0;
            for my $i (2 .. $#tmp){
                my $id = $tmp[$i];
                if ($incl{$id}){
                    push @header, $id;
                    push @keepcols, $i;
                    $hcount++;
                }
            }
            if ($hcount != $inclcount){
                print STDERR "WARNING: $inclcount gennomes given for inclusion (-i), but only $hcount ids found in table\n";
            }
        } else {
            @header = @tmp;
            @keepcols = (2 .. $#tmp);
        }
        $gencount = scalar(@header) - 2;
        print $tout join("\t", @header), "\n";
        next;
    }
    $incount++;
    if ($incount % 1000 == 0){
        print STDERR "\rSNPs in: $incount, SNPs out: $outcount";

    }
    my @array = map{substr($_,0,1)} @tmp[@keepcols];
    my %hash;
    $hash{"F"} = 0;
    $hash{"-"} = 0;
    $hash{$_}++ for @array;
    my $gapcount = $hash{"F"} + $hash{"-"};

    unless (1 - ($gapcount/$gencount) < $minfrac){
        $outcount++;
        print $tout join("\t", @tmp[0,1,@keepcols]), "\n";
        for my $i (0 .. $#array){
            $seq[$i] .= $array[$i];
        }
    }
}
close ($tin);
close ($tout);
print STDERR "\rSNPs in: $incount, SNPs out: $outcount\n";
open (my $aout, ">$pref.alignment_filtered.fasta") or die "ERROR: Can't open alignment file for writing: $!\n";
for my $i (0 .. $#seq){
    my $head = $header[$i+2];
    print $aout ">$head\n$seq[$i]\n";
}
close ($aout);
