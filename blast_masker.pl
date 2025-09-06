#!/usr/bin/env perl

my $version = "0.2";

## Changes from v0.1
## Output to file instead of STDOUT
## Reduced numbers of temporary files
## Fixed bug where fragments at the ends of contigs were being omitted from analysis

use strict;
use warnings;
use File::Basename;

$|++;

my $usage = "
blast_masker.pl

Performs blast of sequences against themselves, outputs regions of self-alignment
in the same interval format as NCBI's dustmaker.

Requires blastn to be in your PATH.

Required:
  -f    Sequence file, in fasta format

Optional:
  -s    size of fragments to break the genome up into when blasting against
        itself, in bp. Enter '0' to not fragment the genome
        (default: 500)
  -i    Minimum sequence identity, in %
        (default: 75)
  -o    output prefix. File will be output as: 
        '<prefix>_blastmasker<fragment size>.txt'
        (default: prefix will be derived from the input filename 
        up to the first '.' character)
  -d    Dustmaker interval file, for comparison
  -t    threads
        (default: 8)

";


use Getopt::Std;
use vars qw( $opt_f $opt_s $opt_i $opt_d $opt_t $opt_o);
getopts('f:s:i:d:t:o:');

die $usage unless $opt_f;

my $seqfile     = $opt_f;
my $fragsize    = defined $opt_s ? $opt_s : 500;
my $minident    = $opt_i ? $opt_i : 75;
my $dustfile    = $opt_d if $opt_d;
my $threads     = $opt_t ? $opt_t : 8;
my $prefix      = $opt_o if defined $opt_o;

unless($prefix){
    my ($dirs, $suffix);
    ($prefix, $dirs, $suffix) = fileparse($seqfile, qr/\..*/);
}
$prefix =~ s/\s/_/g;

## initialize temporary files
my @chars = ("A".."Z", "a".."z");
my $uid;
$uid .= $chars[rand @chars] for 1..8;
my @temp_handles;
for my $i (0 .. ($threads-1)){
    open (my $out, ">temp_fragment.$uid.$i.fasta") or die "ERROR: Can't open temp_fragment file for writing: $!\n";
    push @temp_handles, $out;
}
open (my $kout, ">temp_key.$uid.txt") or die "ERROR: Can't open temp_key file for writing: $!\n";

#Read the sequence file to get the total size and generate fragments
open (my $fin, "<$seqfile") or die "ERROR: Can't open $seqfile: $!\n";
print STDERR "Fragmenting genome...\n";
my ($num_recs, $num_bases) = (0) x 2;
my @order;
my ($id, $seq);
my $filenum = 0;
my $thand = 0;
while (my $line = <$fin>){
    chomp $line;
    if ($line =~ m/^>/){
        if ($id){
            my $totleng = length($seq);
            $num_bases += $totleng;
            if ($fragsize == 0 or $totleng <= $fragsize){
                $filenum++;
                my $out = $temp_handles[$thand];
                print $out ">q$filenum\n$seq\n";
                print $kout "q$filenum,1,$id\n";
                $thand++;
                $thand = 0 if $thand > $#temp_handles;
                #print STDERR "1:$id-$totleng thand$thand q$filenum 1 ($totleng)\n";
            } else {
                my $remain = $totleng;
                my ($start, $stop) = (1, $fragsize);
                while ($remain >= $fragsize){
                    my $sub = substr($seq, $start - 1, $fragsize);
                    $filenum++;
                    my $out = $temp_handles[$thand];
                    print $out ">q$filenum\n$sub\n";
                    print $kout "q$filenum,$start,$id\n";
                    $thand++;
                    $thand = 0 if $thand > $#temp_handles;
                    #print STDERR "2:$id-$remain thand$thand q$filenum $start ($fragsize)\n";
                    $remain -= $fragsize;
                    $start += $fragsize;
                    $stop += $fragsize;
                }
                if ($remain > 0){
                    my $sub = substr($seq, $start - 1, $remain);
                    $filenum++;
                    my $out = $temp_handles[$thand];
                    print $out ">q$filenum\n$sub\n";
                    print $kout "q$filenum,$start,$id\n";
                    $thand++;
                    $thand = 0 if $thand > $#temp_handles;
                    #print STDERR "3:$id-$remain thand$thand q$filenum $start ($remain)\n";
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
        my $out = $temp_handles[$thand];
        print $out ">q$filenum\n$seq\n";
        print $kout "q$filenum,1,$id\n";
        $thand++;
        $thand = 0 if $thand > $#temp_handles;
        #print STDERR "4:$id-$totleng thand$thand q$filenum 1 ($totleng)\n";
    } else {
        my $remain = $totleng;
        my ($start, $stop) = (1, $fragsize);
        while ($remain >= $fragsize){
            my $sub = substr($seq, $start - 1, $fragsize);
            $filenum++;
            my $out = $temp_handles[$thand];
            print $out ">q$filenum\n$sub\n";
            print $kout "q$filenum,$start,$id\n";
            $thand++;
            $thand = 0 if $thand > $#temp_handles;
            #print STDERR "5:$id-$remain thand$thand q$filenum $start ($fragsize) \n";
            $remain -= $fragsize;
            $start += $fragsize;
            $stop += $fragsize;
        }
        if ($remain > 0){
            my $sub = substr($seq, $start - 1, $remain);
            $filenum++;
            my $out = $temp_handles[$thand];
            print $out ">q$filenum\n$sub\n";
            print $kout "q$filenum,$start,$id\n";
            $thand++;
            $thand = 0 if $thand > $#temp_handles;
            #print STDERR "6:$id-$remain thand$thand q$filenum $start ($remain)\n";
        }
    }
}
if ($fragsize){
    print STDERR "Fragmented genome of length $num_bases into $filenum fragments of size $fragsize.\n";
}
foreach my $out (@temp_handles){
    close $out;
}
close ($kout);

#perform blast
my %hits;
my $fcount = 0;
my $num_threads_running = 0;
for my $i (0 .. ($threads-1)){
    my $pid = fork;
    if ($pid == 0){
        my $status = system("blastn -query temp_fragment.$uid.$i.fasta -subject $seqfile -outfmt '6 qaccver pident qstart qend sstart send' > temp_out.$uid.$$.txt");
        unlink "temp_fragment.$uid.$i.fasta";
        exit($status);
        #exit(0);
    }
    $num_threads_running++;
    print STDERR "\rStarting parallel blast job $num_threads_running / ", $threads;
}
print STDERR "\n";

my %hash;
open (my $kin, "<temp_key.$uid.txt") or die "ERROR: Can't open temp_key for reading: $!\n";
while (my $line = <$kin>){
    chomp $line;
    my ($qid, $pos, $id) = split(",", $line);
    $hash{$qid} = "$pos,$id";
}
close ($kin);
unlink "temp_key.$uid.txt";

my $readin = 0;
while (my $pid = wait){
    last if $pid == -1;
    $readin++;
    print STDERR "\rProcessing blast job $readin";
    open (my $pin, "<temp_out.$uid.$pid.txt") or die "ERROR: Can't open temp_out.$uid.$pid.txt: $!\n";
    while (my $line = <$pin>){
        chomp $line;
        my ($qid, $pident, $qstart, $qstop, $sstart, $sstop) = split("\t", $line);
        next if $pident < $minident;
        my ($pos,$id) = split(",", $hash{$qid});
        $qstart += ($pos - 1);
        $qstop  += ($pos - 1);
        next if ($qstart == $sstart and $qstop == $sstop);
        push @{$hits{$id}}, ([$qstart, $qstop]);
    }
    close ($pin);
    unlink("temp_out.$uid.$pid.txt");
    $num_threads_running--;
}
print STDERR "\n";

## Print intervals to output file
my $masked = 0;
my %mcoords;
open (my $maskout, ">", "${prefix}_blastmasker${fragsize}.txt") or die "ERROR: Can't open final interval file for writing: $!\n";
foreach my $id (@order){
    next unless $hits{$id};
    print $maskout ">$id\n";
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
            print $maskout "$seg_start - $seg_stop\n";
            push @{$mcoords{$id}}, ([$seg_start, $seg_stop]) if $dustfile;
            ($seg_start, $seg_stop) = ($start, $stop);
        }
    }
    $masked += ($seg_stop - $seg_start + 1);
    print $maskout "$seg_start - $seg_stop\n";
    push @{$mcoords{$id}}, ([$seg_start, $seg_stop]) if $dustfile;
}
close ($maskout);
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