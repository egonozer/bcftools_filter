#!/usr/bin/perl

my $version = "0.2";

## Will take multiple vcf files as input, filter SNPs based on criteria given,
## then output a fasta alignment of just variant positions. Will also output
## a table of the contig and position of each variant locus.

use strict;
use warnings;

my $usage = "

bcftools_filter_and_align.pl [options]
version: $version

Filters vcf files produced from alignments against a reference and generates
an alignment

Generate the input files by running:
samtools mpileup \
  -E -M 0 -Q 25 -q 30 -m 2 -D -S -g \
  -f <reference.fasta> \
  <alignment.bam> | \
  bcftools view -Icg - \
  > variants.vcf
(for samtools and bcftools versions 0.1.19)

OR

bcftools mpileup \
  -E -Q 25 -q 30 -m 2 -a DP,SP,AD,ADF,ADR \
  -f <reference.fasta> \
  -Ou <alignment.bam> | \
  bcftools call -m -V indels -Ov --ploidy 1 - \
  > variants.vcf
(for bcftools version version 1.9)

Note: As a reminder, the default of samtools mpileup (without the '-B' flag) is
to perform BAQ base quality calculation. Though this can avoid calling false
SNPs around INDELs, it may result in some true bases matching the reference to
be filtered out of the output. Hence there may less false SNPs at the cost of
more false gaps. Caveat emptor.

In the output tab-separated table, the order of the records in each column
is:
1. base call
2. total number of raw reads aligned
3. total number of reads passing filter
4. total number of reads passing filter representing base call
Read depths will be X's in cases where the reference base was called.

Required:
  -v    file of VCF files. Must contain at least 1. VCF files can be gzipped.
        Format:
        /path/to/file.vcf<tab>ID
  -f    fasta file of reference sequence used in alignment

Options:
  -i    include reference sequence in the alignment
        [default: exclude reference from the alignment]
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
        [default: only keep SNVs homozygous under diploid model, i.e. GT = 1/1 or haploid model, i.e. GT = 1]
  -m    reference genome masking file, in interval format output by NCBI dustmaker,
        i.e.
        >chromosome_name
        209 - 215
        415 - 421
        487 - 494
        1104 - 1110
        etc.
  -o    output prefix
        [default: 'output']
  -t    threads
        [default: 2]

";

use Getopt::Std;
use vars qw( $opt_v $opt_f $opt_q $opt_c $opt_d $opt_D $opt_r $opt_h $opt_m $opt_o $opt_x $opt_i $opt_t );
getopts('v:f:q:c:d:D:r:hm:o:x:it:');

die $usage unless $opt_v and $opt_f;

my $vcffile     = $opt_v;
my $reffile     = $opt_f;
my $minqual     = $opt_q ? $opt_q : 200;
my $mincons     = $opt_c ? $opt_c : 75;
my $mindep      = $opt_d ? $opt_d : 5;
my $maxfold     = $opt_D ? $opt_D : 3;
my $mindir      = defined $opt_r ? $opt_r : 1;
my $maskfile    = $opt_m if $opt_m;
my $outid       = $opt_o ? $opt_o: "output";
my $threads     = $opt_t ? $opt_t : 2;

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
print STDERR "Total reference sequences: ", scalar @seqs, "\n";
print STDERR "Total reference sequence length: $totseqleng\n";

#Read and filter the vcf files
my @order;
my %snvs; #contig, position, base, depth
my $num_threads_running = 0;
#stats:
  #below_min_qual(0)
  #below_min_consensus(1)
  #below_min_depth(2)
  #above_max_depth(3)
  #unidirectional(4)
  #non-homozygous(5)
  #masked(6)
  #missing(7)
  #filtered(8)
open (my $vin, "<$vcffile") or die "ERROR: Can't open $vcffile: $!\n";
while (my $vline = <$vin>){
    if ($num_threads_running < $threads){
        chomp $vline;
        my ($path, $gen) = split("\t", $vline);
        my $pid = fork;
        if (0 == $pid){
            my $tpref = "$$";
            open (my $tsout, ">$tpref.stats.txt") or die "ERROR: Can't open $tpref.stats.txt\n";
            open (my $thout, ">$tpref.hash.txt") or die "ERROR: Can't open $tpref.hash.txt\n";
            my $in;
            if ($path =~ m/\.gz$/){
                open ($in, "gzip -cd $path | ");
            } else {
                open ($in, "<$path") or die "ERROR: Can't open $path: $!\n";
            }
            my @stats = (0) x 9;
            my @depths;
            my $last_id = " ";
            my $last_pos = 0;
            my $totalsnps = 0;
            my $progress = 0;
            my $progress_pct = 0;
            while (my $line = <$in>){
                chomp $line;
                next if $line =~ m/^\s*#/;
                my @tmp = split("\t", $line);
                my ($id, $pos, $ref, $alt, $qual, $stuff, $format, $format_val) = @tmp[0,1,3,4,5,7,8,9];
                next if (length($ref) > 1 or length($alt) > 1); #Skip indels
                if ($id ne $last_id){
                    if ($last_pos != 0){
                        my $last_leng = $seqlengs{$last_id};
                        if ($last_pos < $last_leng){
                            for my $i ($last_pos + 1 .. $last_leng){
                                print $thout "$last_id\t$i\t$gen\t-\t0\t0\t0\n";
                                #@{$snvs{$last_id}{$i}{$gen}} = ("-", 0, 0, 0);
                                $stats[7]++;
                            }
                        }
                        for my $i (0 .. $#seqs){
                            my ($ctg, $seq) = @{$seqs[$i]};
                            if ($ctg eq $last_id){
                                my $start = $i + 1;
                                my $nextid = $seqs[$start][0];
                                while ($nextid ne $id){
                                    $start++;
                                    $nextid = $seqs[$start][0];
                                }
                                last;
                            }
                        }
                    }
                    if ($pos > 1){
                        for my $i (1 .. $pos - 1){
                            print $thout "$id\t$i\t$gen\t-\t0\t0\t0\n";
                            #@{$snvs{$id}{$i}{$gen}} = ("-", 0, 0, 0);
                            $stats[7]++;
                        }
                    }
                }
                if ($last_pos > 0 and $last_pos + 1 != $pos){
                    for my $i ($last_pos + 1 .. $pos - 1){
                        print $thout "$id\t$i\t$gen\t-\t0\t0\t0\n";
                        #@{$snvs{$id}{$i}{$gen}} = ("-", 0, 0, 0);
                        $stats[7]++;
                    }
                }
                my @formats = split(":", $format);
                my @format_vals = split(":", $format_val);
                my %fmt;
                for my $i (0 .. $#formats){
                    $fmt{$formats[$i]} = $format_vals[$i];
                }
                (my $fulldepth) = $stuff =~ m/DP=(\d+)/;
                $stuff =~ m/DP4=(\d+),(\d+),(\d+),(\d+)/;
                my ($rf, $rr, $af, $ar) = ($1, $2, $3, $4);
                my $adep = $af+$ar;
                my $total_depth;
                if ($rf or $rr or $af or $ar){
                    $total_depth = $rf + $rr + $af + $ar;
                    print STDERR "WARNING: DP4 does not match DP\n$line\n" if !$fmt{"DP"} or $fmt{"DP"} != $total_depth;
                } else {
                    print STDERR "WARNING: No DP value in FORMAT section:\n$line\n" unless exists $fmt{"DP"};
                    $total_depth = $fmt{"DP"};
                }
                push @depths, $total_depth unless $total_depth == 0;

                ## Count all positions (snp or no snp) with depths below the minimum
                my $filt = 0;
                my @filt_types;
                if ($total_depth < $mindep){
                    print $thout "$id\t$pos\t$gen\tF\t$total_depth\t$adep\t$fulldepth\n";
                    #@{$snvs{$id}{$pos}{$gen}} = ("F", $total_depth, $adep, $fulldepth);
                    $stats[2]++;
                    $filt++;
                    push @filt_types, "bn";
                }
                if ($alt ne "."){

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
                    if (($gtval ne "1/1" and $gtval ne "1") and !$opt_h){
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
                        $stats[8]++;
                        print $thout "$id\t$pos\t$gen\tF\t$total_depth\t$adep\t$fulldepth\n";
                        #@{$snvs{$id}{$pos}{$gen}} = ("F", $total_depth, $adep, $fulldepth);
                    } elsif ($filt == 0) {
                        print $thout "$id\t$pos\t$gen\t$alt\t$total_depth\t$adep\t$fulldepth\n";
                        #@{$snvs{$id}{$pos}{$gen}} = ($alt, $total_depth, $adep, $fulldepth);
                        if ($alt =~ m/,/){
                            my $first = (split",", $alt)[0];
                            print $thout "$id\t$pos\t$gen\t$first\t$total_depth\t$adep\t$fulldepth\n";
                            #@{$snvs{$id}{$pos}{$gen}} = ($first, $total_depth, $adep, $fulldepth);
                        }
                        $totalsnps++;
                    }
                } else {
                    #Eliminate non-SNP positions with high-qual coverage less than the minimum depth

                }
                $last_id = $id;
                $last_pos = $pos;
            }
            close ($in);
            my $last_leng = $seqlengs{$last_id};
            die "\nERROR: No file given\n" unless $last_leng;
            if ($last_pos < $last_leng){
                for my $i ($last_pos + 1 .. $last_leng){
                    print $thout "$last_id\t$i\t$gen\t-\t0\t0\t0\n";
                    #@{$snvs{$last_id}{$i}{$gen}} = ("-", 0, 0, 0);
                    $stats[7]++;
                }
            }

            #calculate median depth
            my $num_depths = scalar @depths;
            print STDERR "num_depths = $num_depths\n";
            @depths = sort{$a <=> $b} @depths;
            my $mid = ($num_depths / 2) - 0.5;
            my $median;
            if ($mid == int($mid)){
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
            foreach my $id (keys %snvs){
                foreach my $pos (keys %{$snvs{$id}}){
                    next unless $snvs{$id}{$pos}{$gen};
                    my ($alt, $total_depth, $adep, $fulldepth) = @{$snvs{$id}{$pos}{$gen}};
                    if ($total_depth > $maxdep){
                        $stats[3]++;
                        $stats[8]++;
                        print $thout "$id\t$pos\t$gen\tF\t$total_depth\t$adep\t$fulldepth\n";
                        #@{$snvs{$id}{$pos}{$gen}} = ("F", $total_depth, $adep, $fulldepth);
                        $totalsnps--;
                    }
                }
            }
            my $totloci = $stats[8] + $totalsnps;
            print $tsout "$gen\t$median\t$stats[7]\t$stats[2]\t$totloci\t$stats[8]\t$stats[0]\t$stats[1]\t$stats[3]\t$stats[4]\t$stats[5]\t$stats[6]\t$totalsnps\n";
            close($thout);
            close($tsout);
            exit($?)
        }
        push @order, ([$gen, $path, $pid]);
        $num_threads_running++;
    }
    if ($num_threads_running == $threads){
        my $pid = wait;
        my $status = $?;
        die "ERROR: Process $pid returned a status of $status\n" if ($status != 0);
        open (my $tin, "<$pid.hash.txt") or die "ERROR: Can't open $pid.hash.txt: $!\n";
        while (my $tline = <$tin>){
            chomp $tline;
            my ($id, $pos, $gen, $alt, $tdep, $adep, $fdep) = split("\t", $tline);
            @{$snvs{$id}{$pos}{$gen}} = ($alt, $tdep, $adep, $fdep);
        }
        close ($tin);
        #unlink ("$pid.hash.txt");
        foreach my $i (0 .. $#order){
            if ($pid == $order[$i][2]){
                print STDERR "Finished processing $order[$i][0]\n";
                last;
            }
        }
        $num_threads_running--;
    }
}
close ($vin);
while (my $pid = wait){
    last if $pid < 0;
    my $status = $?;
    die "ERROR: Process $pid returned a status of $status\n" if ($status != 0);
    open (my $tin, "<$pid.hash.txt") or die "ERROR: Can't open $pid.hash.txt: $!\n";
    while (my $tline = <$tin>){
        chomp $tline;
        my ($id, $pos, $gen, $alt, $tdep, $adep, $fdep) = split("\t", $tline);
        @{$snvs{$id}{$pos}{$gen}} = ($alt, $tdep, $adep, $fdep);
    }
    close ($tin);
    #unlink ("$pid.hash.txt");
    foreach my $i (0 .. $#order){
        if ($pid == $order[$i][2]){
            print STDERR "Finished processing $order[$i][0]\n";
            last;
        }
    }
    $num_threads_running--;
}

open (my $sout, ">$outid.stats.txt");
print $sout "id\tmedian_depth\tmissing\t<min_depth\tsnvs\tfilt_snvs\t<min_qual\t<min_consensus\t>max_depth\tunidir\tnon-homozyg\tmasked\tsnvs_out\n";
foreach my $slice (@order){
    my ($gen, $path, $pid) = @{$slice};
    open (my $in, "<$pid.stats.txt") or die "ERROR: Can't open $pid.stats.txt: $!\n";
    while (my $line = <$in>){
        chomp $line;
        print $sout "$line\n";
    }
    close ($in);
    unlink("$pid.stats.txt");
}
close $sout;

print STDERR "Processing variants\n";
my @results;
foreach (@seqs){
    my ($id, $seq) = @{$_};
    if ($snvs{$id}){
        my %hash = %{$snvs{$id}};
        delete $snvs{$id}; ## free space
        foreach my $pos (sort{$a <=> $b} keys %hash){
            my $refbase = substr($seq, ($pos - 1), 1);
            my @bases;
            my %bhash;
            if ($opt_i){
                push @bases, ([$refbase, "NA", "NA", "NA"]);
                $bhash{$refbase}++;
            }
            foreach my $slice (@order){
                my ($gen, $path) = @{$slice};
                if ($hash{$pos}{$gen}){
                    my ($base, $total_depth, $adep, $fulldepth) = @{$hash{$pos}{$gen}};
                    delete $hash{$pos}{$gen}; # free space
                    push @bases, ([$base, $total_depth, $adep, $fulldepth]);
                    $base = "-" if $base eq "F";
                    $bhash{$base}++;
                } else {
                    push @bases, ([$refbase, "?", "?", "?"]);
                    $bhash{$refbase}++;
                }
            }

            my $bcount = scalar keys %bhash;
            next if $bcount == 1;
            next if ($bcount == 2 and $bhash{"-"});
            my @tmp = ($id, $pos, @bases);
            push @results, [@tmp];
        }
    }
}

print STDERR "Total SNP loci:", scalar @results, "\n";

## Output the alignment sequence and table
print STDERR "Outputting\n";
my @fasta;
open (my $tout, ">$outid.table.txt");
print $tout "id\tpos";
if ($opt_i){
   print $tout "\tref";
}
foreach (@order){
    my ($gen, $path) = @{$_};
    print $tout "\t$gen";
}
print $tout "\n";
foreach my $slice (@results){
    my @tmp = @{$slice};
    my ($contig, $pos, @bases) = @{$slice};
    print $tout "$contig\t$pos";
    for my $i (0 .. $#bases){
        my ($base, $total_depth, $adep, $fulldepth) = @{$bases[$i]};
        print $tout "\t$base,$fulldepth,$total_depth,$adep";
        $base = "-" if $base eq "F";
        $fasta[$i] .= $base;
    }
    print $tout "\n";
}
close ($tout);
open (my $fout, ">$outid.alignment.fasta");
if ($opt_i){
    my $seq = shift @fasta;
    print $fout ">ref\n$seq\n";
}
foreach (@order){
    my ($gen, $path) = @{$_};
    my $seq = shift @fasta;
    print $fout ">$gen\n$seq\n";
}
close ($fout);
