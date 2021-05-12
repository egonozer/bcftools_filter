#!/usr/bin/perl


## Will take multiple vcf files as input, filter SNPs based on criteria given,
## then output a fasta alignment of just variant positions. Will also output
## a table of the contig and position of each variant locus.

use strict;
use warnings;

my $usage = "

bcftools_filter_and_align.pl [options]

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
use vars qw( $opt_v $opt_f $opt_q $opt_c $opt_d $opt_D $opt_r $opt_h $opt_m $opt_o $opt_x $opt_i );
getopts('v:f:q:c:d:D:r:hm:o:x:i');

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
print STDERR "Total reference sequences: ", scalar @seqs, "\n";
print STDERR "Total reference sequence length: $totseqleng\n";

#Read and filter the vcf files
my @order;
my %snvs; #contig, position, base, depth
my %missing_or_filtered;
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
open (my $sout, ">$outid.stats.txt");
print $sout "id\tmedian_depth\tmissing\t<min_depth\tsnvs\tfilt_snvs\t<min_qual\t<min_consensus\t>max_depth\tunidir\tnon-homozyg\tmasked\tsnvs_out\n";
open (my $vin, "<$vcffile") or die "ERROR: Can't open $vcffile: $!\n";
while (my $vline = <$vin>){
    chomp $vline;
    my ($path, $gen) = split("\t", $vline);
    push @order, ([$gen, $path]);
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
    print STDERR "Reading $gen: $progress_pct%";
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
                        #push @{$missing_or_filtered{$last_id}{$i}{$gen}}, "mi";
                        @{$snvs{$last_id}{$i}{$gen}} = ("-", 0, 0, 0, 0);
                        $stats[7]++;
                        #push @depths, 0;
                        $progress += 100*(1/$totseqleng);
                        if (int($progress) != $progress_pct){
                            $progress_pct = int($progress);
                            print STDERR "\rReading $gen: $progress_pct%";
                        }
                    }
                }
                for my $i (0 .. $#seqs){
                    my ($ctg, $seq) = @{$seqs[$i]};
                    if ($ctg eq $last_id){
                        my $start = $i + 1;
                        my $nextid = $seqs[$start][0];
                        while ($nextid ne $id){
                            $progress += 100*(length($seqs[$start][1])/$totseqleng);
                            if (int($progress) != $progress_pct){
                                $progress_pct = int($progress);
                                print STDERR "\rReading $gen: $progress_pct%";
                            }
                            $start++;
                            $nextid = $seqs[$start][0];
                        }
                        last;
                    }
                }
            }
            if ($pos > 1){
                for my $i (1 .. $pos - 1){
                    #push @{$missing_or_filtered{$id}{$i}{$gen}}, "mi";
                    @{$snvs{$id}{$i}{$gen}} = ("-", 0, 0, 0, 0);
                    $stats[7]++;
                    #push @depths, 0;
                    $progress += 100*(1/$totseqleng);
                    if (int($progress) != $progress_pct){
                        $progress_pct = int($progress);
                        print STDERR "\rReading $gen: $progress_pct%";
                    }
                }
            }
        }
        if ($last_pos > 0 and $last_pos + 1 != $pos){
            for my $i ($last_pos + 1 .. $pos - 1){
                #push @{$missing_or_filtered{$id}{$i}{$gen}}, "mi";
                @{$snvs{$id}{$i}{$gen}} = ("-", 0, 0, 0, 0);
                $stats[7]++;
                #push @depths, 0;
                $progress += 100*(1/$totseqleng);
                if (int($progress) != $progress_pct){
                    $progress_pct = int($progress);
                    print STDERR "\rReading $gen: $progress_pct%";
                }
            }
        }
        $progress += 100*(1/$totseqleng);
        if (int($progress) != $progress_pct){
            $progress_pct = int($progress);
            print STDERR "\rReading $gen: $progress_pct%";
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
        my $rdep = $rf+$rr;
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
            push @{$missing_or_filtered{$id}{$pos}{$gen}}, "bn";
            @{$snvs{$id}{$pos}{$gen}} = ("F", $total_depth, $rdep, $adep, $fulldepth);
            $stats[2]++;
            $filt++;
            push @filt_types, "bn";
        }

        #push @depths, $total_depth;
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
                push @{$missing_or_filtered{$id}{$pos}{$gen}}, @filt_types;
                $stats[8]++;
                @{$snvs{$id}{$pos}{$gen}} = ("F", $total_depth, $rdep, $adep, $fulldepth);
            } elsif ($filt == 0) {
                @{$snvs{$id}{$pos}{$gen}} = ($alt, $total_depth, $rdep, $adep, $fulldepth);
                if ($alt =~ m/,/){
                    my $first = (split",", $alt)[0];
                    @{$snvs{$id}{$pos}{$gen}} = ($first, $total_depth, $rdep, $adep, $fulldepth);
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
            #push @{$missing_or_filtered{$last_id}{$i}{$gen}}, "mi";
            @{$snvs{$last_id}{$i}{$gen}} = ("-", 0, 0, 0, 0);
            $stats[7]++;
            #push @depths, 0;
            $progress += 100*(1/$totseqleng);
            if (int($progress) != $progress_pct){
                $progress_pct = int($progress);
                print STDERR "\rReading $gen: $progress_pct%";
            }
        }
    }
    print STDERR "\rReading $gen: 100%\n";

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
        #print STDERR "$id\n";
        foreach my $pos (keys %{$snvs{$id}}){
            #print STDERR "^$pos\n";
            next unless $snvs{$id}{$pos}{$gen};
            my ($alt, $total_depth, $rdep, $adep, $fulldepth) = @{$snvs{$id}{$pos}{$gen}};
            #print STDERR "$id $pos $alt $depth\n";
            if ($total_depth > $maxdep){
                $stats[3]++;
                $stats[8]++;
                push @{$missing_or_filtered{$id}{$pos}{$gen}}, "ad";
                @{$snvs{$id}{$pos}{$gen}} = ("F", $total_depth, $rdep, $adep, $fulldepth);
                $totalsnps--;
            }
        }
    }
    my $mi_and_bn = $stats[7] + $stats[2];
    print STDERR "Total missing: $stats[7]\n";
    print STDERR "Below minimum depth ($mindep): $stats[2]\n";
    print STDERR "Total missing or below minimum depth ($mindep): $mi_and_bn\n";
    my $pct_covered = 100 * (($totseqleng - $mi_and_bn) / $totseqleng);
    print STDERR "\tPercent aligned: $pct_covered\n\n";

    my $totloci = $stats[8] + $totalsnps;
    print STDERR "Total SNV loci: $totloci\n";
    print STDERR "Total filtered SNVs: $stats[8]\n";
    #below_min_qual(0), below_min_consensus(1), below_min_depth(2), above_max_depth(3), unidirectional(4), non-homozygous(5), masked(6),
    print STDERR "\tBelow minimum quality ($minqual): $stats[0]\n";
    print STDERR "\tBelow minimum consensus ($mincons%): $stats[1]\n";
    print STDERR "\tAbove maximum depth ($maxfold x $median): $stats[3]\n";
    print STDERR "\tUnidirectional ($mindir): $stats[4]\n";
    print STDERR "\tNon-homozygous: $stats[5]\n" if !$opt_h;
    print STDERR "\tMasked: $stats[6]\n" if $maskfile;
    print STDERR "Total unfiltered SNV loci: $totalsnps\n\n";

    print $sout "$gen\t$median\t$stats[7]\t$stats[2]\t$totloci\t$stats[8]\t$stats[0]\t$stats[1]\t$stats[3]\t$stats[4]\t$stats[5]\t$stats[6]\t$totalsnps\n";

}
close $sout;

my @results;
foreach (@seqs){
    my ($id, $seq) = @{$_};
    if ($snvs{$id}){
        my %hash = %{$snvs{$id}};
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
                    my ($base, $total_depth, $rdep, $adep, $fulldepth) = @{$hash{$pos}{$gen}};
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
