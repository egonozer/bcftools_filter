#!/usr/bin/perl

use strict;
use warnings;

$|++;

my $usage = "
blast_masker.pl

Performs blast of sequences against themselves, outputs regions of self-alignment
in the same interval format as NCBI's dustmaker.

Required:
  -f    Sequence file, in fasta format
  
Optional:
  -s    size of fragments to break the genome up into when blasting against
        itself, in bp. Enter '0' to not fragment the genome
        (default: 500)
  -i    Minimum sequence identity, in %
        (default: 75)
  -d    Dustmaker interval file, for comparison
  -t    threads
        (default: 15)

";


use Getopt::Std;
use vars qw( $opt_f $opt_s $opt_i $opt_d $opt_t);
getopts('f:s:i:d:t:');

die $usage unless $opt_f;

my $seqfile     = $opt_f;
my $fragsize    = defined $opt_s ? $opt_s : 500;
my $minident    = $opt_i ? $opt_i : 75;
my $dustfile    = $opt_d if $opt_d;
my $threads     = $opt_t ? $opt_t : 15;

#Read the sequence file to get the total size and generate fragments
open (my $fin, "<$seqfile") or die "ERROR: Can't open $seqfile: $!\n";
my ($num_recs, $num_bases) = (0) x 2;
my @order;
my @fragments;
my ($id, $seq);
my $filenum = 0;
while (my $line = <$fin>){
    chomp $line;
    if ($line =~ m/^>/){
        if ($id){
            my $totleng = length($seq);
            $num_bases += $totleng;
            if ($fragsize == 0 or $totleng <= $fragsize){
                $filenum++;
                open (my $out, ">temp_fragment.$filenum.fasta") or die "ERROR: Can't open file for writing: $!\n";
                print $out ">$id\n$seq\n";
                close ($out);
                push @fragments, ([$id, 1, $filenum]);
            } else {
                my $remain = $totleng;
                my ($start, $stop) = (1, $fragsize);
                while ($remain >= $fragsize){
                    my $sub = substr($seq, $start - 1, $fragsize);
                    $filenum++;
                    open (my $out, ">temp_fragment.$filenum.fasta") or die "ERROR: Can't open file for writing: $!\n";
                    print $out ">$id\n$sub\n";
                    close ($out);
                    push @fragments, ([$id, $start, $filenum]);
                    $remain -= $fragsize;
                    $start += $fragsize;
                    $stop += $fragsize;
                }
                if ($remain > 0){
                    ($start, $stop) = (($totleng - $fragsize) + 1, $totleng);
                    my $sub = substr($seq, $start - 1, $fragsize);
                    $filenum++;
                    open (my $out, ">temp_fragment.$filenum.fasta") or die "ERROR: Can't open file for writing: $!\n";
                    print $out ">$id\n$sub\n";
                    close ($out);
                    push @fragments, ([$id, $start, $filenum]);
                }
            }
        }
        $seq = "";
        $num_recs++;
        $id = substr($line, 1);
        $id =~ s/\s.*$//;
        push @order, $id;
        next;
    }
    $line =~ s/\s//g;
    next unless $line;
    $seq .= $line;
}
close ($fin);
if ($id){
    my $totleng = length($seq);
    $num_bases += $totleng;
    if ($fragsize == 0 or $totleng <= $fragsize){
        $filenum++;
        open (my $out, ">temp_fragment.$filenum.fasta") or die "ERROR: Can't open file for writing: $!\n";
        print $out ">$id\n$seq\n";
        close ($out);
        push @fragments, ([$id, 1, $filenum]);
    } else {
        my $remain = $totleng;
        my ($start, $stop) = (1, $fragsize);
        while ($remain >= $fragsize){
            my $sub = substr($seq, $start - 1, $fragsize);
            $filenum++;
            open (my $out, ">temp_fragment.$filenum.fasta") or die "ERROR: Can't open file for writing: $!\n";
            print $out ">$id\n$sub\n";
            close ($out);
            push @fragments, ([$id, $start, $filenum]);
            $remain -= $fragsize;
            $start += $fragsize;
            $stop += $fragsize;
        }
        if ($remain > 0){
            ($start, $stop) = (($totleng - $fragsize) + 1, $totleng);
            my $sub = substr($seq, $start - 1, $fragsize);
            $filenum++;
            open (my $out, ">temp_fragment.$filenum.fasta") or die "ERROR: Can't open file for writing: $!\n";
            print $out ">$id\n$sub\n";
            close ($out);
            push @fragments, ([$id, $start, $filenum]);
        }
    }
}
if ($fragsize){
    print STDERR "Fragmented genome of length $num_bases into ".scalar @fragments." fragments of size $fragsize.\n";
}

#perform blast
my %hits;
my $fcount = 0;
my $num_threads_running = 0;
foreach my $slice (@fragments){
    if ($num_threads_running == $threads){
        my $pid = wait;
        open (my $pin, "<temp_out.$pid.txt") or die "ERROR: Can't open temp_out.$pid.txt: $!\n";
        while (my $line = <$pin>){
            chomp $line;
            my ($id, $qstart, $qstop) = split(",", $line);
            push @{$hits{$id}}, ([$qstart, $qstop]);
        }
        close ($pin);
        unlink("temp_out.$pid.txt");
        $num_threads_running--;
    }
    my ($id, $start, $fnum) = @{$slice};
    $fcount ++;
    print STDERR "\rBlasting fragment $fcount / ".scalar @fragments;
    my $pid = fork;
    if ($pid == 0){
        open (my $tout, ">temp_out.$$.txt") or exit(1);
        open (my $bin, "/usr/local/bin/blastn -query temp_fragment.$fnum.fasta -subject $seqfile -outfmt 6 | ");
        while (my $line = <$bin>){
            chomp $line;
            my @tmp = split("\t", $line);
            my ($qid, $sid, $pctid, $qstart, $qstop, $sstart, $sstop) = @tmp[0,1,2,6,7,8,9];
            $qstart += ($start - 1);
            $qstop += ($start - 1);
            next if ($qid eq $sid and $qstart == $sstart and $qstop == $sstop);
            next if $pctid < $minident;
            print $tout "$id,$qstart,$qstop\n";
            #push @{$hits{$id}}, ([$qstart, $qstop]);
        }
        close ($bin);
        close ($tout);
        unlink("temp_fragment.$fnum.fasta");
        exit(0);
    }
    $num_threads_running++;
}
print STDERR "\n";
while (my $pid = wait){
    last if $pid == -1;
    open (my $pin, "<temp_out.$pid.txt") or die "ERROR: Can't open temp_out.$pid.txt: $!\n";
    while (my $line = <$pin>){
        chomp $line;
        my ($id, $qstart, $qstop) = split(",", $line);
        push @{$hits{$id}}, ([$qstart, $qstop]);
    }
    close ($pin);
    unlink("temp_out.$pid.txt");
    $num_threads_running--;
}
my $masked = 0;
my %mcoords;
foreach my $id (@order){
    next unless $hits{$id};
    print ">$id\n";
    my @coords = @{$hits{$id}};
    @coords = sort{$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @coords;
    my ($seg_start, $seg_stop) = @{shift @coords};
    while (@coords){
        my ($start, $stop) = @{shift @coords};
        if ($start <= ($seg_stop + 1)){
            $seg_stop = $stop if ($stop >= $seg_stop);
            next;
        } else {
            $masked += ($seg_stop - $seg_start + 1);
            print "$seg_start - $seg_stop\n";
            push @{$mcoords{$id}}, ([$seg_start, $seg_stop]);
            ($seg_start, $seg_stop) = ($start, $stop);
        }
    }
    $masked += ($seg_stop - $seg_start + 1);
    print "$seg_start - $seg_stop\n";
    push @{$mcoords{$id}}, ([$seg_start, $seg_stop]);
}
print STDERR "Masked $masked positions out of $num_bases total.\n";

if ($dustfile){
    open (my $din, "<$dustfile") or die "ERROR: Can't open $dustfile: $!\n";
    my %dcoords;
    my $id;
    while (my $line = <$din>){
        chomp $line;
        if ($line =~ m/^>/){
            $id = substr($line, 1);
            $id =~ s/\s.*$//;
            next;
        }
        $line =~ m/(\d+)\s-\s(\d+)/;
        if ($1 and $2){
            push @{$dcoords{$id}}, ([$1, $2]);
        }
    }
    close ($din);
    
    my $dnotm = 0;
    foreach my $id (keys %dcoords){
        my @dc = @{$dcoords{$id}};
        my %mpos;
        if ($mcoords{$id}){
            my @mc = @{$mcoords{$id}};
            foreach my $slice (@mc){
                my ($start, $stop) = @{$slice};
                for my $i ($start .. $stop){
                    $mpos{$i} = 1;
                }
            }
            
        }
        foreach my $slice (@dc){
            my ($start, $stop) = @{$slice};
            for my $i ($start .. $stop){
                $dnotm++ unless $mpos{$i};
            }
        }
    }
    print STDERR "$dnotm positions found in dustmaker intervals not found in blast_masker intervals\n";
}
