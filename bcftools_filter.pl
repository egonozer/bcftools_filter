#!/usr/bin/perl

use strict;
use warnings;

my $usage = "

bcftools_filter.pl [options] <results.vcf or piped from bcftools>

Filters bcftools variant results produced from alignments against a reference.

Generate the input file by running:
samtools mpileup -E -M 0 -Q 25 -q 30 -m 2 -D -S -g -f <reference.fasta> <alignment.bam> | bcftools view -Icg - 
(for samtools and bcftools versions 0.1.19)

Will also accept gzipped or bgzipped vcf files

Output will be sequence file of reference strain where SNVs passing filters
are substituted into the sequence.  Missing data or filtered SNVs will be
replaced with gap characters '-'

Note: As a reminder, the default of samtools mpileup (without the '-B' flag) is
to perform BAQ base quality calculation. Though this can avoid calling false
SNPs around INDELs, it may result in some true bases matching the reference to
be filtered out of the output. Hence there may less false SNPs at the cost of
more false gaps. Caveat emptor.

Required:
  -f    fasta file for reference strain

Options:
  -q    minimum SNV quality score
        [default: 200]
  -c    minimum read consensus, in percent
        [default: 75]
  -d    minimum read depth
          * will be calculated from DP4 to only count reads above quality filter
        [default: 5]
  -D    maximum read depth, in fold of median for the isolate
        [default: 3]
  -r    minimum number of reads in each direction
        [default: 1]
  -h    DO NOT require homozygositiy
        [default: only keep SNVs homozygous under diploid model, i.e. GT = 1/1]
  -m    reference genome masking file, in interval format output by NCBI dustmaker,
        i.e.
        >chromosome_name
        209 - 215
        415 - 421
        487 - 494
        1104 - 1110
        etc.
  -o    Output sequence ID. If there is only one sequence in the reference, then
        its name will be replaced in the output. If there are two or more
        sequences in the reference, their IDs will be prefixed.
        (default: output will have reference sequence IDs)
  -x    Filename of difference file to output. Will output a file with differences
        from the reference sequence and what type of filter was applied.
        Output file columns will be:
        - reference position
        - base change
        - filters applied
            bq = below minimum SNV quality score
            bc = below minimum read consensus
            bd = below minimum read depth (at snp site)
            bn = below minimum read depth (at non-snp site)
            ad = above maximum read depth
            br = below minimum number of reads in each direction
            nh = not homozygous
            ma = masked
            mi = missing
        
";

use Getopt::Std;
use vars qw( $opt_f $opt_q $opt_c $opt_d $opt_D $opt_r $opt_h $opt_m $opt_o $opt_x );
getopts('f:q:c:d:D:r:hm:o:x:');

die $usage unless @ARGV and $opt_f;

my $reffile     = $opt_f;
my $minqual     = $opt_q ? $opt_q : 200;
my $mincons     = $opt_c ? $opt_c : 75;
my $mindep      = $opt_d ? $opt_d : 5;
my $maxfold     = $opt_D ? $opt_D : 3;
my $mindir      = $opt_r ? $opt_r : 1;
my $maskfile    = $opt_m if $opt_m;
my $outid       = $opt_o if $opt_o;
my $difffile    = $opt_x if $opt_x;

my %mask;
if ($maskfile){
    my $totmask = 0;
    my $masklines = 0;
    open (my $in, "<$maskfile") or die "ERROR: Can't open $maskfile: $!\n";
    my $id;
    while (my $line = <$in>){
        chomp $line;
        if ($line =~ m/^>/){
            $id = substr($line, 1);
            $id =~ s/\s.*$//;
            next;
        }
        if ($line =~ m/^(\d+) - (\d+)/){
            $masklines++;
            my ($start, $stop) = ($1, $2);
            for my $i ($start .. $stop){
                $mask{$id}{$i} = 1;
                $totmask++;
            }
        }
    }
    close ($in);
    print STDERR "Total of $totmask masked positions in $masklines intervals.\n";
}

## Read in reference sequence
open (my $fin, "<$reffile") or die "ERROR: Can't open $reffile: $!\n";
my @seqs;
my %seqlengs;
my $totseqleng = 0;
my ($id, $seq);
while (my $line = <$fin>){
    chomp $line;
    next if $line =~ m/^\s*$/;
    if ($line =~ m/^>/){
        if ($id){
            push @seqs, ([$id, $seq]);
            my $leng = length($seq);
            $seqlengs{$id} = $leng;
            $totseqleng += $leng;
            print STDERR "$id: $leng bp\n";
            $seq = "";
        }
        $id = substr($line, 1);
        $id =~ s/\s.*$//;
        next;
    }
    $line =~ s/\s//g;
    $seq .= $line;
}
close ($fin);
if ($id){
    push @seqs, ([$id, $seq]);
    my $leng = length($seq);
    $seqlengs{$id} = $leng;
    $totseqleng += $leng;
    print STDERR "$id: $leng bp\n";
    $seq = "";
}
print STDERR "Total sequence length: $totseqleng\n";

#Read and filter the bcftools results
my $file = $ARGV[0];
my $in;
if ($file =~ m/\.bgz$/){
    open ($in, "bgzip -cd $file | ")
} elsif ($file =~ m/\.gz$/){
    open ($in, "gzip -cd $file | ")
} else {
    open ($in, "<file") or die "ERROR: Can't open $file: $!\n";
}
my %snvs; #contig, position, base, depth
my %missing_or_filtered;
#stats: below_min_qual(0), below_min_consensus(1), below_min_depth(2), above_max_depth(3), unidirectional(4), non-homozygous(5), masked(6), missing(7), filtered(8)
my @stats = (0) x 9;
my @depths;
my $last_id = " ";
my $last_pos = 0;
while (my $line = <$in>){
    chomp $line;
    next if $line =~ m/^\s*#/;
    my @tmp = split("\t", $line);
    my ($id, $pos, $alt, $qual, $stuff, $format, $format_val) = @tmp[0,1,4,5,7,8,9];
    if ($id ne $last_id){
        if ($last_pos != 0){
            my $last_leng = $seqlengs{$last_id};
            if ($last_pos < $last_leng){
                for my $i ($last_pos + 1 .. $last_leng){
                    push @{$missing_or_filtered{$last_id}{$i}}, "mi";
                    $stats[7]++;
                    #push @depths, 0;
                }
            }
        }
        if ($pos > 1){
            for my $i (1 .. $pos - 1){
                push @{$missing_or_filtered{$id}{$i}}, "mi";
                $stats[7]++;
                #push @depths, 0;
            }
        }
    }
    if ($last_pos > 0 and $last_pos + 1 != $pos){
        for my $i ($last_pos + 1 .. $pos - 1){
            push @{$missing_or_filtered{$id}{$i}}, "mi";
            $stats[7]++;
            #push @depths, 0;
        }
    }
    my @formats = split(":", $format);
    my @format_vals = split(":", $format_val);
    my %fmt;
    for my $i (0 .. $#formats){
        $fmt{$formats[$i]} = $format_vals[$i];
    }
    $stuff =~ m/DP4=(\d+),(\d+),(\d+),(\d+)/;
    my ($rf, $rr, $af, $ar) = ($1, $2, $3, $4);
    my $total_depth;
    if ($rf or $rr or $af or $ar){
        $total_depth = $rf + $rr + $af + $ar;
        print STDERR "WARNING: DP4 does not match DP\n$line\n" if !$fmt{"DP"} or $fmt{"DP"} != $total_depth;
    } else {
        print STDERR "WARNING: No DP value in FORMAT section:\n$line\n" unless exists $fmt{"DP"};
        $total_depth = $fmt{"DP"};
    }
    push @depths, $total_depth unless $total_depth == 0;
    #push @depths, $total_depth;
    if ($alt ne "."){
        my $filt = 0;
        my @filt_types;
        if ($qual < $minqual){
            $filt++;
            $stats[0]++;
            push @filt_types, "bq";
        }
        my $cons = 100*(($af + $ar) / $total_depth);
        if ($cons < $mincons){
            $filt++;
            $stats[1]++;
            push @filt_types, "bc";
        }
        if ($total_depth < $mindep){
            $filt++;
            $stats[2]++;
            push @filt_types, "bd";
            $stats[7]++; #will also count these as missing / uncovered
        }
        if ($af < $mindir or $ar < $mindir){
            $filt++;
            $stats[4]++;
            push @filt_types, "br";
        }
        my $gtval;
        if ($fmt{"GT"}){
            $gtval = $fmt{"GT"};
        } else {
            print STDERR "WARNING: Variant position without GT value:\n$line\n";
        }
        if ($gtval ne "1/1" and $gtval ne "1" and !$opt_h){
            $filt++;
            $stats[5]++;
            push @filt_types, "nh";
        }
        if ($mask{$id}{$pos}){
            $filt++;
            $stats[6]++;
            push @filt_types, "ma";
        }
        if ($filt > 0){
            push @{$missing_or_filtered{$id}{$pos}}, @filt_types;
            $stats[8]++;
        } elsif ($filt == 0) {
            @{$snvs{$id}{$pos}} = ($alt, $total_depth);
            if ($alt =~ m/,/){
                my $first = (split",", $alt)[0];
                @{$snvs{$id}{$pos}} = ($first, $total_depth);
            }
        }
    } else {
        #Eliminate non-SNP positions with high-qual coverage less than the minimum depth
        if ($total_depth < $mindep){
            push @{$missing_or_filtered{$id}{$pos}}, "bn";
            $stats[7]++;
        }
    }
    $last_id = $id;
    $last_pos = $pos;
}
close ($in);
my $last_leng = $seqlengs{$last_id};
die "ERROR: No file given\n" unless $last_leng;
if ($last_pos < $last_leng){
    for my $i ($last_pos + 1 .. $last_leng){
        push @{$missing_or_filtered{$last_id}{$i}}, "mi";
        $stats[7]++;
        #push @depths, 0;
    }
}

#calculate median depth
my $num_depths = scalar @depths;
print STDERR "num_depths = $num_depths\n";
@depths = sort{$a <=> $b} @depths;
my $mid = ($num_depths / 2) - 0.5;
my $median;
if ($mid == int($mid)){
    print STDERR "mid = $mid\n";
    $median = $depths[$mid];
} else {
    my $two = $num_depths / 2;
    my $one = $two - 1;
    $median = ($depths[$one] + $depths[$two]) / 2;
}
print STDERR "Median depth: $median\n";
print STDERR "Maximum depth: $depths[$#depths]\n";
my $maxdep = $median * $maxfold;
print STDERR "Maximum depth threshold: $maxdep\n";
#Filter SNVs above maximum depth
my $totalsnps = 0;
foreach my $id (keys %snvs){
    #print STDERR "$id\n";
    foreach my $pos (keys %{$snvs{$id}}){
        #print STDERR "^$pos\n";
        my ($alt, $depth) = @{$snvs{$id}{$pos}};
        #print STDERR "$id $pos $alt $depth\n";
        if ($depth > $maxdep){
            $stats[3]++;
            $stats[8]++;
            push @{$missing_or_filtered{$id}{$pos}}, "ad";
            delete $snvs{$id}{$pos};
        } else {
            $totalsnps++;
        }
    }
}

print STDERR "Total missing or below minimum depth ($mindep): $stats[7]\n";
my $pct_covered = 100 * (($totseqleng - $stats[7]) / $totseqleng);
print STDERR "\tPercent aligned: $pct_covered\n";
print STDERR "Total filtered SNVs: $stats[8]\n";
#below_min_qual(0), below_min_consensus(1), below_min_depth(2), above_max_depth(3), unidirectional(4), non-homozygous(5), masked(6),
print STDERR "\tBelow minimum quality ($minqual): $stats[0]\n";
print STDERR "\tBelow minimum consensus ($mincons%): $stats[1]\n";
print STDERR "\tBelow minimum depth ($mindep): $stats[2]\n";
print STDERR "\tAbove maximum depth ($maxfold x $median): $stats[3]\n";
print STDERR "\tUnidirectional ($mindir): $stats[4]\n";
print STDERR "\tNon-homozygous: $stats[5]\n" if !$opt_h;
print STDERR "\tMasked: $stats[6]\n";
print STDERR "Total SNVs: $totalsnps\n";

print STDERR "\nOutputting replaced sequence\n";
my @sdiff;
my @fdiff;
foreach my $slice (@seqs){
    my ($id, $seq) = @{$slice};
    if ($outid){
        my $newid = "$outid\_$id";
        $newid = $outid if scalar @seqs == 1;
        print ">$newid\n";
    } else {
        print ">$id\n";
    }
    my @tmp = split(//, $seq);
    for my $i (0 .. $#tmp){
        my $pos = $i + 1;
        my $base = $tmp[$i];
        if (my $snp = $snvs{$id}{$pos}){
            $base = @{$snp}[0];
            if ($opt_x){
                unless ($missing_or_filtered{$id}{$pos}){
                    my $string = "$pos\t$base";
                    push @sdiff, $string;
                }
            }
            #print STDERR "pos:". ($i + 1) . " SNV\n";
        }
        if ($missing_or_filtered{$id}{$pos}){
            $base = "-";
            if ($opt_x){
                my @tmp = @{$missing_or_filtered{$id}{$pos}};
                @tmp = sort{$a cmp $b}@tmp;
                my $string = "$pos\t$base\t";
                $string .= join(",", @tmp);
                push @fdiff, $string;
            }
            #print STDERR "pos:". ($i + 1) . " missing\n";
            #print STDERR "*** BOTH ***\n" if $snvs{$id}{$i+1};
        }
        print "$base";
    }
    print "\n";
}
if ($opt_x){
    open (my $diffout, ">$difffile") or die "ERROR: Can't open $difffile for writing: $!\n";
    print $diffout join("\n", @sdiff), "\n";
    print $diffout join("\n", @fdiff), "\n";
    close ($diffout);
}

