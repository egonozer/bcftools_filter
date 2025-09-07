#!/usr/bin/perl

use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;
use File::Spec::Functions qw ( catfile path );

my $usage = "
bcftools_filter_and_align_annotate.pl

Required:
  -t    table.txt file output by bcftools_filter_and_align
  -g    genbank file of annotations
  
";

## command line processing.
use Getopt::Std;
use vars qw( $opt_t $opt_g );
getopts('t:g:');
die $usage unless $opt_t and $opt_g;

my ($tfile, $gbkfile) = ($opt_t, $opt_g);

## Read in the genbank file
my @c_seqs;
my %c_recs;
my $status = gbk_convert($gbkfile);
die "ERROR: gbk_convert failed" if $status;

## Read in SNP table
my $header;
my $last_contig;
my $cseq;
my @refrecs;
my %codontable = translate();
open (my $tin, "<$tfile") or die "ERROR: Can't open $tfile: $!\n";
while (my $line = <$tin>){
    chomp $line;
    my @tmp = split("\t", $line);
    my $cid = shift @tmp;
    my $pos = shift @tmp;
    unless ($header){
        print "contig_id\tpos\tref_base\ttype\tlocus_id\tproduct\tstrand\tgene_position\tcodon_position\tcodon\tAA";
        foreach my $i (0 .. $#tmp){
            my $id = $tmp[$i];
            #if ($i == 0){
            #    print "\tref_base\ttype\tlocus_id\tproduct\tstrand\tgene_position\tcoodn_position\tcodon\tAA"
            #} else {
            #    print "\t$id\_base\t$id\_aa"
            #}
            print "\t$id\_base\t$id\_aa"
        }
        print "\n";
        $header = 1;
        next;
    }
    print "$cid\t$pos";
    if (!$last_contig or $cid ne $last_contig){
        @refrecs = ();
        my $temp;
        for my $i (0 .. $#c_seqs){
            if ($c_seqs[$i][0] eq $cid){
                $temp = $i;
                last;
            }
        }
        if (defined $temp){
            $cseq = uc($c_seqs[$temp][1]);
            foreach my $start (sort{$a <=> $b} keys %{$c_recs{$cid}}){
                foreach my $stop (sort{$a <=> $b} keys %{$c_recs{$cid}{$start}}){
                    foreach my $dir (sort keys %{$c_recs{$cid}{$start}{$stop}}){
                        my ($type, $lid, $prod) = @{$c_recs{$cid}{$start}{$stop}{$dir}};
                        push @refrecs, ([$start, $stop, $dir, $type, $lid, $prod]);
                    }
                }
            }
        } else {
            print "\n";
            next;
        }
        $last_contig = $cid;
    }
    
    my $zerpos = $pos - 1;
    #my $refbase = (split(",", $tmp[0]))[0];
    my $refbase = substr($cseq, $zerpos, 1);
    my @vars;
    for my $i (0 .. $#tmp){
        my $base = (split(",", $tmp[$i]))[0];
        push @vars, $base;
    }
    
    my @hits;
    foreach my $slice (@refrecs){
        my ($start, $stop, $dir, $type, $lid, $prod) = @{$slice};
        last if ($pos > $stop and @hits);
        if ((sort{$a <=> $b}($start, $stop, $pos))[1] == $pos){
            ## is found within a record
            if ($type eq "CDS"){
                my $genseq = substr($cseq, ($start-1), ($stop - $start + 1));
                my $genpos = ($pos - $start) + 1;
                if ($dir eq "-"){
                    $genseq =~ tr/ACGT/TGCA/;
                    $genseq = reverse($genseq);
                    $genpos = ($stop - $pos) + 1;
                }
                my $cpos = 3;
                my $codon = substr($genseq, ($genpos - 3), 3);
                if ($genpos % 3 == 1){
                    $cpos = 1;
                    $codon = substr($genseq, ($genpos - 1), 3);
                } elsif ($genpos % 3 == 2){
                    $cpos = 2;
                    $codon = substr($genseq, ($genpos - 2), 3);
                }
                my $genaa = $codontable{$codon};
                push @hits, ([$start, $stop, $dir, $type, $lid, $prod, $genpos, $cpos, $codon, $genaa]);
            } else {
                push @hits, ([$start, $stop, $dir, $type, $lid, $prod]);
            }
        }
    }
    #position\tref_base\ttype\tlocus_id\tproduct\tgene_position\tcoodn_position\tcodon\tAA";
    unless (@hits){
        print "\t$refbase\tintergenic\t\t\t\t\t\t\t";
        foreach (@vars){
            print "\t$_\t";
        }
        print "\n";
        next;
    }
    foreach (@hits){
        my ($start, $stop, $dir, $type, $lid, $prod, $genpos, $cpos, $codon, $genaa) = @{$_};
        print "\t$refbase\t$type\t$lid\t$prod\t$dir";
        if ($type eq "CDS"){
            print "\t$genpos\t$cpos\t$codon\t$genaa";
        } else {
            print "\t\t\t\t";
        }
        foreach my $base (@vars){
            if ($dir eq "-"){
                $base =~ tr/ACGT/TGCA/;
            }
            print "\t$base";
            if ($type eq "CDS"){
                my $qcodon = $codon;
                substr($qcodon, ($cpos - 1), 1) = $base;
                my $aa = "X";
                $aa = $codontable{$qcodon} if $codontable{$qcodon};
                print "\t$aa";
            } else {
                print "\t";
            }
        }
        print "\n";
    }
    
}






#------------------------------------------------------------
sub translate {
    #table is gencode 11 taken from http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG11
    my $aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    my $b1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
    my $b2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
    my $b3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
    my @aa_a = split("", $aa);
    my @b1_a = split("", $b1);
    my @b2_a = split("", $b2);
    my @b3_a = split("", $b3);
    my %out;
    for my $i (0 .. $#aa_a){
        $out{"$b1_a[$i]$b2_a[$i]$b3_a[$i]"} = $aa_a[$i];
    }
    return (%out);
}

sub gbk_convert{
    my $file = shift;
    my $filename = shift;
    my $filenum = shift;
    my $return_status = 0;
    my $shortfile = basename($file);
    return(1) unless -e $file;
    my $gbkin;
    if ($file =~ m/\.gz$/){
        open ($gbkin, "gzip -cd $file | ") or return(1);
    } else {
        return(5) if -B $file;
        open ($gbkin, "<", $file) or return(1);
    }
    my $loccount = 0;
    my $seqcount = 0;
    my ($c_id, $c_seq);
    my $is_prod;
    my @tags;
    my %crecs;
    my @ctg_order;
    my $reading = 1; # 1 = front material, 2 = annotations, 3 = sequence
    
    while (my $fline = <$gbkin>){
        $fline =~ s/\R/\012/g; #converts to UNIX-style line endings
        my @lines = split("\n", $fline); #need to split lines by line-ending character in the case of Mac-formatted files which only have CR line terminators, not both CR and LF like DOS
        while (@lines){
            my $line = shift @lines;
            next if $line =~ m/^\s*$/;
            if ($line =~ m/^LOCUS\s+\S*\s+\d+\sbp/){
                if ($reading == 2){ #no ORIGIN sequence record was found between LOCUS records
                    return (2);
                }
                if ($reading == 3){
                    if ($c_seq and $c_id){
                        #print $seqout ">#$filename#$c_id\n$c_seq\n";
                        push @c_seqs, ([$c_id, $c_seq]);
                        $c_seq = "";
                        $reading = 1;
                    } else {
                        return (2);
                    }
                }
            }
            if ($line =~ m/^\/\//){ #reached the end of the file (or record)
                if ($c_seq and $c_id){
                    #print $seqout ">#$filename#$c_id\n$c_seq\n";
                    push @c_seqs, ([$c_id, $c_seq]);
                    $c_seq = "";
                    $reading = 1;
                } else {
                    return (2);
                }
            }
            if ($reading == 1){
                if ($line =~ m/^LOCUS\s+([^\s]+)/){
                    $seqcount++;
                    if ($line =~ m/^LOCUS\s+(\S+)\s+\d+ bp/){
                        $c_id = $1;
                    } else {
                        $c_id = "rec$seqcount";
                    }
                    push @ctg_order, $c_id;
                    next;
                }
                if ($line =~ m/^FEATURES\s+Location\/Qualifiers/){
                    $reading = 2;
                    next;
                }
            } elsif ($reading == 2){
                if ($line =~ m/^\s+(\S+)\s+(complement\()*[<>]*(\d+)<*\.\.[<>]*(\d+)>*\)*\s*$/){
                    $is_prod = "";
                    my ($type, $start, $stop) = ($1, $3, $4);
                    my $dir = "+";
                    $dir = "-" if $2;
                    unless ($type eq "source" or $type eq "gene" or $crecs{$c_id}{$start}{$stop}{$dir}){
                        @{$crecs{$c_id}{$start}{$stop}{$dir}} = ($type);
                    }
                    if ($type eq "CDS"){
                        #${$crecs{$c_id}{$start}{$stop}{$dir}}[0] = 1;
                        $loccount++;
                    }
                    if ($tags[0]){
                        my ($o_start, $o_stop, $o_dir) = @tags;
                        ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[1] = $tags[5] if $tags[5]; #use gene id if that's all there is
                        ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[1] = $tags[3] if $tags[3]; #prefer locus id
                        ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[2] = $tags[4] if $tags[4];
                        $loccount++;
                    }
                    @tags = ($start, $stop, $dir);
                    if ($type eq "source" or $type eq "gene"){
                        undef @tags;                        
                    }
                    next;
                }
                if ($line =~ m/^ORIGIN\s*$/){
                    $is_prod = "";
                    if ($tags[0]){
                        my ($o_start, $o_stop, $o_dir) = @tags;
                        ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[1] = $tags[5] if $tags[5]; #use gene id if that's all there is
                        ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[1] = $tags[3] if $tags[3]; #prefer locus id
                        ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[2] = $tags[4] if $tags[4];
                        $loccount++;
                    }
                    undef @tags;
                    $reading = 3;
                    next
                }
                if ($line =~ m/^\s+\/(\S+)=\"*([^"]*)\"*/){
                    $is_prod = "";
                    my ($key, $val) = ($1, $2);
                    if ($key eq "locus_tag"){
                        $tags[3] = $val;
                    }
                    if ($key eq "product"){
                        $tags[4] = $val;
                        $is_prod = 1;
                    }
                    if ($key eq "gene"){
                        $tags[5] = $val;
                    }
                    next;
                }
                if ($is_prod){
                    $line =~ s/^\s*//;
                    $line =~ s/"*\s*$//;
                    $tags[4] .= " $line";
                }
            } elsif ($reading == 3){
                $line =~ s/\d//g;
                $line =~ s/\s//g;
                $c_seq .= $line;
                next;
            }
        }
    }
    if ($c_seq and $c_id){
        #print $seqout ">#$filename#$c_id\n$c_seq\n";
        $c_seq = "";
        $reading = 1;
    }
    close ($gbkin);
    %c_recs = %crecs;
    return(0);
}
