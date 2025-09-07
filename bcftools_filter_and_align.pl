#!/usr/bin/perl

my $version = "0.3";

## Changes from v0.2
## Replace large snv hash with arrays (memory saving)
## Second forking step for post-processing followed by rejoining output files. Should save lots of memory for very large data sets
## Removed redundancies and improve efficiency

## Will take multiple vcf files as input, filter SNPs based on criteria given,
## then output a fasta alignment of just variant positions. Will also output
## a table of the contig and position of each variant locus.

use strict;
use warnings;
use List::Util qw/sum/;

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

Note: As a reminder, the default of samtools or bcftools mpileup (without the 
'-B' flag) is to perform BAQ base quality calculation. Though this can avoid 
calling false SNPs around INDELs, it may result in some true bases matching the 
reference to be filtered out of the output. Hence there may fewer false SNPs at 
the cost of more false gaps. Caveat emptor.

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
  -H    DO NOT require homozygositiy
        [default: only keep SNVs homozygous under diploid model, i.e. GT = 1/1 or haploid model, i.e. GT = 1]
  -m    reference genome masking file, in interval format output by blast_masker.pl or NCBI dustmaker,
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
use vars qw( $opt_v $opt_f $opt_q $opt_c $opt_d $opt_D $opt_r $opt_H $opt_m $opt_o $opt_x $opt_i $opt_t );
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

## Read in reference sequence and assign dummies
open (my $fin, "<$reffile") or die "ERROR: Can't open $reffile: $!\n";
my @seqs;
my %seqlengs;
my %seqdum;
my %seqdum_rev;
my ($seqleng, $totseqleng, $maxleng) = (0) x 3;
my ($id, $id_tell, $seq_tell);
while (my $line = <$fin>){
    my $tell = tell $fin;
    if ($line =~ m/^>/){
        if ($id){
            my $distance = $seq_tell - $id_tell;
            push @seqs, ([$id, $id_tell, $distance]);
            my $order = $#seqs;
            $seqdum{$id} = $order;
            $seqdum_rev{$order} = $id;
            $seqlengs{$order} = $seqleng;
            $maxleng = $seqleng if $seqleng >= $maxleng;
            $totseqleng += $seqleng;
            print STDERR "$id: $seqleng bp\n";
            $seqleng = 0;
        }
        $id_tell = $tell;
        chomp $line;
        $id = substr($line, 1);
        $id =~ s/\s.*$//;
        next;
    }
    $seq_tell = $tell;
    chomp $line;
    $line =~ s/\s//g;
    $seqleng += length($line);
}
close ($fin);
if ($id){
    my $distance = $seq_tell - $id_tell;
    push @seqs, ([$id, $id_tell, $distance]);
    my $order = $#seqs;
    $seqdum{$id} = $order;
    $seqdum_rev{$order} = $id;
    $seqlengs{$order} = $seqleng;
    $maxleng = $seqleng if $seqleng >= $maxleng;
    $totseqleng += $seqleng;
    print STDERR "$id: $seqleng bp\n";
    $seqleng = 0;
}
print STDERR "Total reference sequences: ", scalar @seqs, "\n";
print STDERR "Total reference sequence length: $totseqleng, maximum length: $maxleng\n";
## Set number of leading zeroes for position and sequence numbers
my $lzp = int((log($maxleng)/log(10))+1);
my $lzs = int((log(scalar @seqs)/log(10))+1);

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
                $mask{"$id-$i"} = 1;
                $totmask++;
            }
        }
    }
    close ($in);
    print STDERR "Total of $totmask masked positions in $masklines intervals.\n";
}

#Read and filter the vcf files
my $vtot = 0;
open (my $vin, "<$vcffile") or die "ERROR: Can't open $vcffile: $!\n";
while (my $vline = <$vin>) {
    chomp $vline;
    next if $vline =~ m/^\s*$/;
    $vtot++;
}
close ($vin);

## calculate number of leading zeros needed in temporary names based on total number of input vcf files to process.  
my $lzv = int((log($vtot)/log(10))+1);
## Unique temporary id for keeping temporary files from other potential parallel processes on in a cluster separate. 
my $tempid = $$;
print STDERR "Process ID: $tempid\n";

my %baseidx = ( A => 0, C => 1, G => 2, T => 3);
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
print STDERR "Starting vcf file processing\n";
my ($pcount, $vcount) = (0) x 2;
my %pid_hash;
open ($vin, "<$vcffile") or die "ERROR: Can't open $vcffile: $!\n";
while (my $vline = <$vin>){
    if ($num_threads_running < $threads){
        chomp $vline;
        next if $vline =~ m/^\s*$/;
        my ($path, $gen) = split("\t", $vline);
        $pcount++;
        my $tpref = sprintf("%0${lzv}d", $pcount);
        my $pid = fork;
        if (0 == $pid){
            #my $tpref = "$$";
            open (my $tsout, ">tstats.$tempid.$tpref.txt") or die "ERROR: Can't open tstats.$tempid.$tpref.txt for writing: $!\n";
            open (my $thout, ">thash.$tempid.$tpref.txt") or die "ERROR: Can't open thash.$tempid.$tpref.txt for writing: $!\n";
            my @tsnvs;
            my $in;
            if ($path =~ m/\.gz$/){
                open ($in, "gzip -cd $path | ");
            } else {
                open ($in, "<$path") or die "ERROR: Can't open $path ($gen): $!\n";
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
                my $id_dum = $seqdum{$id};
                next if (length($ref) > 1 or length($alt) > 1); #Skip indels
                if ($id_dum ne $last_id){
                    if ($last_pos != 0){
                        my $last_leng = $seqlengs{$last_id};
                        if ($last_pos < $last_leng){
                            for my $i ($last_pos + 1 .. $last_leng){
                                print $thout "$last_id\t$i\t$gen\t-\t0\t0\t0\n";
                                $stats[7]++;
                            }
                        }
                    }
                    if ($pos > 1){
                        for my $i (1 .. $pos - 1){
                            print $thout "$id_dum\t$i\t$gen\t-\t0\t0\t0\n";
                            $stats[7]++;
                        }
                    }
                }
                if ($last_pos > 0 and $last_pos + 1 != $pos){
                    for my $i ($last_pos + 1 .. $pos - 1){
                        print $thout "$id_dum\t$i\t$gen\t-\t0\t0\t0\n";
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
                my $printed;
                if ($total_depth < $mindep){
                    print $thout "$id_dum\t$pos\t$gen\tF\t$total_depth\t$adep\t$fulldepth\n";
                    $stats[2]++;
                    $filt++;
                    push @filt_types, "bn";
                    $printed = 1;
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
                    if (($gtval ne "1/1" and $gtval ne "1") and !$opt_H){
                        $filt++;
                        $stats[5]++;
                        push @filt_types, "nh";
                    }
                    if ($mask{"$id-$pos"}){
                        $filt++;
                        $stats[6]++;
                        push @filt_types, "ma";
                    }
                    if ($filt > 0){
                        $stats[8]++;
                        print $thout "$id_dum\t$pos\t$gen\tF\t$total_depth\t$adep\t$fulldepth\n" unless $printed;
                    } elsif ($filt == 0) {
                        if ($alt =~ m/,/){
                            $alt = (split",", $alt)[0];
                        }
                        push @tsnvs, "$id_dum,$pos,$alt,$total_depth,$adep,$fulldepth";
                        $totalsnps++;
                    }
                } else {
                    #Eliminate non-SNP positions with high-qual coverage less than the minimum depth

                }
                $last_id = $id_dum;
                $last_pos = $pos;
            }
            close ($in);
            my $last_leng = $seqlengs{$last_id};
            die "\nERROR: No file given (id:$gen, last_id:$last_id)\n" unless $last_leng;
            if ($last_pos < $last_leng){
                for my $i ($last_pos + 1 .. $last_leng){
                    print $thout "$last_id\t$i\t$gen\t-\t0\t0\t0\n";
                    $stats[7]++;
                }
            }

            #calculate median depth
            my $num_depths = scalar @depths;
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
            ##### Should it be median depth per chromosome / contig???? (think plasmids, organelles ...)
            my $maxdep = $median * $maxfold;
            #Filter SNVs above maximum depth
            foreach (@tsnvs){
                my ($id_dum, $pos, $alt, $total_depth, $adep, $fulldepth) = split(",", $_);
                if ($total_depth > $maxdep){
                    $stats[3]++;
                    $stats[8]++;
                    $alt = "F";
                    $totalsnps--;
                }
                print $thout "$id_dum\t$pos\t$gen\t$alt\t$total_depth\t$adep\t$fulldepth\n";
            }
            my $totloci = $stats[8] + $totalsnps;
            print $tsout "$gen\t$median\t$stats[7]\t$stats[2]\t$totloci\t$stats[8]\t$stats[0]\t$stats[1]\t$stats[3]\t$stats[4]\t$stats[5]\t$stats[6]\t$totalsnps\n";
            close($thout);
            close($tsout);
            exit($?)
        }
        @{$pid_hash{$pid}} = ($gen, $tpref);
        push @order, ([$gen, $path, $tpref]);
        $num_threads_running++;
    }
    if ($num_threads_running == $threads){
        my $pid = wait;
        my $status = $?;
        my ($gen, $tpref) = @{$pid_hash{$pid}};
        delete $pid_hash{$pid};
        die "ERROR: Process $pid ($gen prefix $tpref) returned a status of $status\n" if ($status != 0);
        open (my $tin, "<thash.$tempid.$tpref.txt") or die "ERROR: Can't open thash.$tempid.$tpref.txt: $!\n";
        while (my $tline = <$tin>){
            chomp $tline;
            my ($id_dum, $pos, $gen, $alt, $tdep, $adep, $fdep) = split("\t", $tline);
            my $idlz = sprintf("%0${lzs}d", $id_dum);
            my $polz = sprintf("%0${lzp}d", $pos);
            unless ($snvs{"$idlz-$polz"}){
                $snvs{"$idlz-$polz"} = "0,0,0,0,0";
            }
            my @basearray = split(",", $snvs{"$idlz-$polz"});
            my $bpos = exists $baseidx{uc($alt)} ? $baseidx{uc($alt)} : 4;
            $basearray[$bpos]++;
            $snvs{"$idlz-$polz"} = join(",", @basearray);
        }
        close ($tin);
        $vcount++;
        print STDERR "Finished pre-processing $gen ($vcount/$vtot)\n";
        $num_threads_running--;
    }
}
close ($vin);
while (my $pid = wait){
    last if $pid < 0;
    my $status = $?;
    my ($gen, $tpref) = @{$pid_hash{$pid}};
    delete $pid_hash{$pid};
    die "ERROR: Process $pid ($gen prefix $tpref) returned a status of $status\n" if ($status != 0);
    open (my $tin, "<thash.$tempid.$tpref.txt") or die "ERROR: Can't open thash.$tempid.$tpref.txt: $!\n";
    while (my $tline = <$tin>){
        chomp $tline;
        my ($id_dum, $pos, $gen, $alt, $tdep, $adep, $fdep) = split("\t", $tline);
        my $idlz = sprintf("%0${lzs}d", $id_dum);
        my $polz = sprintf("%0${lzp}d", $pos);
        unless ($snvs{"$idlz-$polz"}){
            $snvs{"$idlz-$polz"} = "0,0,0,0,0";
        }
        my @basearray = split(",", $snvs{"$idlz-$polz"});
        my $bpos = exists $baseidx{uc($alt)} ? $baseidx{uc($alt)} : 4;
        $basearray[$bpos]++;
        $snvs{"$idlz-$polz"} = join(",", @basearray);
    }
    close ($tin);
    $vcount++;
    print STDERR "Finished pre-processing $gen ($vcount/$vtot)\n";
    $num_threads_running--;
}
undef %mask;

print STDERR "Total unique sites: ", scalar keys %snvs, "\n";
print STDERR "Counting variant sites ... ";
my ($seqidx, $seq) = (-1) x 2;
foreach my $idpos (sort {$a cmp $b} keys %snvs){
    my @array = split(",", $snvs{$idpos});
    my $gapsum = pop @array;
    my $basecount = 0;
    foreach (@array){
        $basecount ++ if $_ > 0;
    }
    my $basesum = sum(@array);
    my $skip;
    $skip = 1 if $basecount == 0;
    if ($basesum + $gapsum == $vtot){
        if (!$opt_i){
            $skip = 1 if $basecount == 1;
        }
    }
    if ($skip){
        delete $snvs{$idpos};
        next;
    }
    my ($id_dum, $pos) = split("-", $idpos);
    my $id_dum_nz = 0 + $id_dum;
    while ($id_dum_nz != $seqidx){
        $seqidx++;
        if ($id_dum_nz == $seqidx){
            my ($seqid, $idx, $dist) = @{$seqs[$seqidx]};
            open (my $in, "$reffile");
            seek $in, $idx, 0;
            read $in, $seq, $dist;
            close ($in);
            $seq =~ s/\s//g;
        }
    }
    my $refbase = substr($seq, ($pos - 1), 1);
    $snvs{$idpos} = $refbase;
}
print STDERR "total variant sites:", scalar keys %snvs, "\n";

## process hashes in parallel and stats
$num_threads_running = 0;
$vcount = 0;
my %npid_hash;
open (my $sout, ">$outid.stats.txt") or die "ERROR: Can't open $outid.stats.txt for writing: $!\n";
print $sout "id\tmedian_depth\tmissing\t<min_depth\tsnvs\tfilt_snvs\t<min_qual\t<min_consensus\t>max_depth\tunidir\tnon-homozyg\tmasked\tsnvs_out\n";
foreach my $j (0 .. $#order){
    if ($num_threads_running < $threads){
        my ($gen, $path, $tpref) = @{$order[$j]};
        my $npid = fork;
        if (0 == $npid) {
            open (my $tin, "<thash.$tempid.$tpref.txt") or die "ERROR: Can't open thash.$tempid.$tpref.txt: $!\n";
            my @tsnvs;
            while (my $tline = <$tin>){
                chomp $tline;
                my ($id, $pos, $gen, $alt, $tdep, $adep, $fdep) = split("\t", $tline);
                my $idlz = sprintf("%0${lzs}d", $id);
                my $polz = sprintf("%0${lzp}d", $pos);
                push @tsnvs, "$idlz-$polz#$alt,$fdep,$tdep,$adep" if $snvs{"$idlz-$polz"};
            }
            close ($tin);
            @tsnvs = sort{$a cmp $b} @tsnvs;

            open (my $pout, ">tproc.$tempid.$tpref.txt") or die "ERROR: Can't open tproc.$tempid.$tpref.txt ($gen) for writing: $!\n";
            my ($tidpos, $tval) = ("X") x 2;
            ($tidpos, $tval) = split("#", shift @tsnvs) if @tsnvs;
            foreach my $idpos (sort{$a cmp $b} keys %snvs){
                my $refbase = $snvs{$idpos};
                my $val = "$refbase,?,?,?";
                if ($tidpos eq $idpos){
                    $val = $tval;
                    ($tidpos, $tval) = split("#", shift @tsnvs) if @tsnvs;
                }
                print $pout "$val\n";
            }
            close ($pout);
            unlink("thash.$tempid.$tpref.txt");
            exit;
        }
        $num_threads_running++;
        $npid_hash{$npid} = $gen;
        ## process stats
        open (my $in, "<tstats.$tempid.$tpref.txt") or die "ERROR: Can't open tstats.$tempid.$tpref.txt: $!\n";
        while (my $line = <$in>){
            chomp $line;
            print $sout "$line\n";
        }
        close ($in);
        unlink("tstats.$tempid.$tpref.txt");
    }
    if ($num_threads_running == $threads){
        my $npid = wait;
        my $status = $?;
        die "ERROR: Process $npid ($npid_hash{$npid}) returned a status of $status\n" if ($status != 0);
        $num_threads_running--;
        $vcount++;
        print STDERR "\rFinished post-processing $vcount/$vtot";
    }
}
close ($sout);
while (my $npid = wait){
    last if $npid < 0;
    my $status = $?;
    die "ERROR: Process $npid ($npid_hash{$npid}) returned a status of $status\n" if ($status != 0);
    $num_threads_running--;
    $vcount++;
    print STDERR "\rFinished post-processing $vcount/$vtot";
}

## Process and output table
my $start_sites = scalar keys %snvs;
print STDERR "\nOutputting table: 0/$start_sites";
open (my $tout, ">$outid.table.txt") or die "ERROR: Can't open $outid.table.txt for writing: $!\n";
print $tout "id\tpos";
if ($opt_i){
   print $tout "\tref";
}
my @handarray;
foreach (@order){
    my ($gen, $path, $tpref) = @{$_};
    print $tout "\t$gen";
    open (my $in, "<tproc.$tempid.$tpref.txt") or die "ERROR: Can't open tproc.$tempid.$tpref.txt ($gen) for reading: $!\n";
    push @handarray, $in;
}
print $tout "\n";
my $in_snvs = 0;
my $total_snvs = 0;
foreach my $idpos (sort{$a cmp $b} keys %snvs) {
    $in_snvs++;
    print STDERR "\rOutputting table: $in_snvs/$start_sites" if $in_snvs % 100 == 0;
    my ($id_dum, $pos) = split("-", $idpos);
    $id_dum+=0;
    $pos+=0;
    #my ($id_dum, $pos, $refbase) = @{$ref_results[$i]};
    my $id = $seqdum_rev{$id_dum};
    my $outline = "$id\t$pos";
    $outline .= "\t$snvs{$idpos}" if $opt_i;
    for my $j (0 .. $#order){
        chomp(my $val = readline $handarray[$j]);
        $outline .= "\t$val";
    }
    print $tout "$outline\n";
    $total_snvs++;
}
print STDERR "\rOutputting table: $in_snvs/$start_sites\n";

## Output alignment
print STDERR "Outputting alignment\n";
open (my $fout, ">$outid.alignment.fasta") or die "ERROR: Can't open $outid.alignment.fasta for writing: $!\n";
if ($opt_i){
    print $fout ">ref\n";
    foreach my $idpos (sort {$a cmp $b} keys %snvs) {
        print $fout "$snvs{$idpos}";
    }
    print $fout "\n";
}
my $outcount = 0;
for my $i (0 .. $#order){
    my ($gen, $path, $tpref) = @{$order[$i]};
    close ($handarray[$i]); ## close the open filehandles from table step above
    print $fout ">$gen\n";
    my $count = -1;
    open (my $in, "<tproc.$tempid.$tpref.txt") or die "ERROR: Can't open tproc.$tempid.$tpref.txt ($gen) for reading: $!\n";
    while (my $line = <$in>){
        $count++;
        my $base = substr($line, 0, 1);
        $base = "-" if $base eq "F";
        print $fout "$base";
    }
    print $fout "\n";
    $outcount++;
    print STDERR "$gen ($outcount/", scalar @order, ")\n";
    close ($in);
    unlink ("tproc.$tempid.$tpref.txt");
}
close ($fout);
print STDERR "\nDone!\n";