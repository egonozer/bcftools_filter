# BCFTOOLS FILTER
***

## 1. INTRODUCTION

This repository contains a set of scripts for reading and filtering variant call format (vcf) files produced by bcftools from whole genome sequence alignments to generate consensus sequences or consensus alignments.  

## 2. REQUIREMENTS

* Perl
* [bcftools](https://anaconda.org/bioconda/bcftools)
* [blastn](https://anaconda.org/bioconda/blast) (only needed for [blast_masker.pl](#47-blast_maskerpl))

## 3. INSTALLATION

Download and save.

## 4. SCRIPTS

### 4.1 bcftools_filter.pl 

Filters bcftools variant result produced from sequencing read alignment against a reference and outputs a single consensus sequence. Output will be sequence file of reference strain where SNVs passing filters are substituted into the sequence.  Missing data or filtered SNVs will be replaced with gap characters '-'. 

Only single nucleotide variants will be identified. Indels are not included in the consensus sequence generated so all consensus sequences will have the same lengths as the reference seqeunce(s).

#### 4.1.1 Generating input files

Generate variant files for input into bcftools_filter.pl from SAM or BAM format read alignments by running (for samtools and bcftools versions 0.1.19):

```
samtools mpileup \
  -E -M 0 -Q 25 -q 30 -m 2 -D -S -g \
  -f <reference.fasta> \
  <alignment.bam> | \
  bcftools view -Icg - \
  > variants.vcf
```

OR

For bcftools version version 1.0 or higher:

```
bcftools mpileup \
  -E -Q 25 -q 30 -m 2 -a DP,SP,AD,ADF,ADR \
  -f <reference.fasta> \
  -Ou <alignment.bam> | \
  bcftools call -m -V indels -Ov --ploidy 1 - \
  > variants.vcf
```

In either case, "<reference.fasta>" should be replaced with the path to the reference sequence file in fasta format that was used go generate the SAM/BAM alignment file and "<alignment.bam>" should be replaced with the path to the alignment bam file.  

> [!NOTE]
> As a reminder, the default of samtools mpileup (without the '-B' flag) is to perform BAQ base quality calculation. Though this can avoid calling false SNPs around INDELs, it may result in some true bases matching the reference to be filtered out of the output. Hence there may fewer false SNPs at the cost of more false gaps. Caveat emptor.

> [!NOTE]
> `mpileup` output can also be piped directly to bcftools_filter.pl without the need to first create the vcf file. 

#### 4.1.2 Required Inputs:

> [!TIP]
> For list of options, call the script without any inputs: `perl bcftools_filter.pl`

* `-f`
Path to fasta file for reference strain used to generate the original SAM/BAM alignment

* `/path/to/<vcf file>`
The last argument in the command should either be the path to the vcf file created in step 4.1.1 above OR if piping directly from `mpileup` the command should end with the `-` character.

FOR EXAMPLE:

```
perl bcftools_filter.pl -f reference.fasta variants.vcf > consensus.fasta
```

Or if piping from `mpileup`:

```
bcftools mpileup \
  -E -Q 25 -q 30 -m 2 -a DP,SP,AD,ADF,ADR \
  -f reference.fasta \
  -Ou alignment.bam | \
  bcftools call -m -V indels -Ov --ploidy 1 - | \
  perl bcftools_filter.pl -f reference.fasta - > consensus.fasta
```

#### 4.1.3 Optional Inputs:

* `-q`: Minimum SNV quality score. Default: 200
* `-c`: Minimum read consensus, in percent. Default: 75
* `-d`: Minimum read depth. Will be calculated from the DP4 values in the VCF record to only count reads above quality filter. Default: 5
* `-D`: Maximum read depth, in fold of median read for the total alignment. Default: 3
* `-r`: Minimum number of reads in either direction required. Default: 1
* `-H`: DO NOT require homozygositiy. Default: only keep SNVs homozygous under diploid model, i.e. GT = 1/1 or haploid model, i.e. GT = 1. 
* `-m`: Reference genome masking file, in interval format output by blast_masker.pl ([see below](#47-blast_maskerpl)) or NCBI dustmaker.
Example format: 
```
        >chromosome_name
        209 - 215
        415 - 421
        487 - 494
        1104 - 1110
```
* `-o`: Output sequence ID. If there is only one sequence in the reference, then its name will be replaced with this value in the output. If there are two or more sequences in the reference, their IDs will be prefixed with this value. Default: output consensus sequence will have the reference sequence IDs.
* `-x`: Filename of difference file to output. Will output a file with differences from the reference sequence and what type of filter was applied.

#### 4.1.4 Outputs

Consensus sequence will be output to STDOUT.

Alignment and filtering statistics will be output to STDERR.

If you provide a filename to the `-x` option, a difference file will be output with details of each variant position and filtering applied.

Difference file columns will represent:

1. reference position

2. base change

3. filters applied:

            bq = below minimum SNV quality score
            bc = below minimum read consensus
            bd = below minimum read depth (at snp site)
            bn = below minimum read depth (at non-snp site)
            ad = above maximum read depth
            br = below minimum number of reads in each direction
            nh = not homozygous
            ma = masked
            mi = missing  

### 4.2 bcftools_filter_and_align.pl

Filters vcf files produced from alignments against a reference similar to bcftools_filter.pl above. The difference is that bcftools_filter_and_align.pl takes as input multiple one or more vcf files generated from alignments against the same reference and the output is a multiple sequence alignment of single nucleotide variant positions relative to the reference. 

#### 4.2.1 Generating input files

Refer to section 4.1.1 above for generating VCF files for each of the alignments to be filtered. 

#### 4.2.2 Required inputs

* `-f`: Path to fasta file for the reference strain used to generate the original SAM/BAM alignments.
* `-v`: Path to tab-separated file of VCF files and a unique ID for each alignment. Must contain at least 1 entry. VCF files can be gzipped.
Format:
```
/path/to/file.vcf<tab>ID
```

#### 4.2.3 Optional inputs

* `-i`: Include the reference sequence in the alignment. Will be labeled "ref" in the output. Default: exclude reference from the alignment
* `-q`: Minimum SNV quality score. Default: 200
* `-c`: Minimum read consensus, in percent. Default: 75
* `-d`: Minimum read depth. Will be calculated from the DP4 values in the VCF record to only count reads above quality filter. Default: 5
* `-D`: Maximum read depth, in fold of median read for the total alignment. Default: 3
* `-r`: Minimum number of reads in either direction required. Default: 1
* `-H`: DO NOT require homozygositiy. Default: only keep SNVs homozygous under diploid model, i.e. GT = 1/1 or haploid model, i.e. GT = 1.
* `-m`: Reference genome masking file, in interval format output by blast_masker.pl ([see below](#47-blast_maskerpl)) or NCBI dustmaker.
Example format: 
```
        >chromosome_name
        209 - 215
        415 - 421
        487 - 494
        1104 - 1110
```
* `-o`: Output files prefix. Dfault: 'output'
* `-t`: Parallel threads. Default: 2


#### 4.2.4 Outputs

Three files are output by bcftools_filer_and_align.pl: 

**\<prefix>.alignment.fasta**

Multiple sequence alignment of variant positions only. Sequence headers will be the IDs from the file of vcf paths and IDs given to `-v`. If `-i` was selected, the reference genome sequence will be included under the header ">ref". 

Filtered variant positions from each input file will be replaced with a "-" character in the output sequence alignment. 

**\<prefix>.stats.txt**

A tab-separated table of alignment and filtering results. Columns are:

1. id: Alignment ID
2. median_depth: Median read depth across all positions in the alignment
3. missing: Positions with zero read depth in the alignment
4. <min_depth: Positions with read depths below the threshold given to option `-d`
5. snvs: Total SNV positions relative to the reference that were identified in VCF file
6. filt_snvs: Number of SNV postitions that were flagged by one or more filter criteria
7. <min_qual: SNVs positions filtered for SNV quality scores below `-q` threshold.	
8. <min_consensus: SNV psositions filtered for being below the minimum consensus threshold `-c`.	
9. \>max_depth: SNV positions filtered for depths more than X-fold above the median where X was given to `-D`
10. unidir: SNV positions filtered for being covered by fewer than X forward and reverse reads where X was given by `-r`	
11. non-homozyg: SNV positions filtered for being identified as non-homozygous, i.e. without GT flags of '1/1' or '1', unless `-H` option was given. 
12. masked: SNV positions filtered for falling within an interval in a masking file given to option `-m`	
13. snvs_out: Total number of SNV positions relative to the reference sequence after applying filters 

**\<prefix>.table.txt**

Tab-separated column file of all variant positions output. Only positions with a variant in at least one of the input vcf files relative to the other input genomes (or relative to the reference sequence if `-i` is given) will be output.

For each position and isolate, a series of four values separated by commas will be given. The comma-separated values represent, in order:

1. base call. This will be "F" if it was filtered by one or more of the filter criteria.
2. total number of raw reads aligned to the postiion. 
3. number of reads passing filter
4. number of reads passing filter that represented the base call

Read depths will be ?'s in cases where the reference base was called.

### 4.3 bcftools_filter_and_align_annotate.pl

Generates table of annotations of SNVs output by bcftools_filter_and_align. 

#### 4.3.1 Required inputs

* `-t`: table.txt file output by bcftools_filter_and_align.pl
* `-g`: Genbank-formatted file (i.e. .gbk or .gbff) containing annotations of the sequence used as the original alignment reference. LOCUS record IDs in the genbank file must match the sequence IDs in the fasta file used to create the original SAM/BAM alignment files. 

#### 4.3.2 Outputs

A tab-separated table file.

The first 11 columns apply to the reference sequence:

1. contig_id: Reference sequence ID
2. pos: Position in the reference sequence
3. ref_base: Nucleotide base in the reference sequence
4. type: Coding sequence (CDS) or non-coding sequence.
5. locus_id: Locus ID of the genome feature
6. product: If the feature is a coding sequence, the annotated gene product
7. strand: Strand on which the feature is encoded, i.e. forward (+) or reverse (-)
8. gene_position: Nucleotide position within the gene 
9. codon_position: Position in the codon triplet, if applicable
10. codon: Codon sequence
11. AA: Translated amino acid sequence

All following columns apply to each of the aligned sequences:

12. \<id>_base: Nucleotide base for alignment \<id> at this position
13. \<id>_aa: Translated amino acid in \<id> 

etc. 

### 4.4 bcftools_filter_and_align_filter.pl

Filters a table produced by bcftools_filter_and_align.pl to remove
non-core sites, i.e. those without a variant position called in a given subset of the included alignments, and generate a new table and new fasta alignment.

#### 4.4.1 Required inputs

* `-t`: table.txt file output by bcftools_filter_and_align.pl

#### 4.4.2 Optional inputs

* `-m`: Minimum fraction of genomes in which a SNV position must have a base called in order to be counted (i.e. \"1\" for SNV position found in all genomes, \"0.7\" for SNV position found in at least 70% of the genomes). Default: any SNP position found in at least 2 of the genomes will be included, i.e. \"0\"
* `-i`: File with list of genomes to include in the ouput. List should be genome names, one per line. Fraction limits given by -m will be based on these genomes (i.e. if -m is 0.5 and you enter list of 4 genomes for -i, output will be loci present in at least 2 of these 4 genomes). Default: all genomes will be included
* `-o`: Output file prefix. Default: 'output'

#### 4.4.3 Outputs

Will produce files \<prefix>.table_filtered.txt and \<prefix>.alignment_filtered.fasta where "\<prefix>" is the value given to option `-o`. See section  [4.2.4](#424-outputs) above for details. 

### 4.5 bcftools_filter_and_align_invariant.pl

Calculates invariant site counts from the output of bcftools_filter_and_align.pl. When using iqtree2 to generate a phylogenetic tree from the consensus alignment fasta, the output of this script can be used with the '-fconst' option to simulate a complete genome alignment while saving time and memory. 

**USAGE:**

```
perl bcftools_filter_and_align_invariant.pl <reference.fasta> <table.txt>
```

#### 4.5.1 Required inputs

* \<reference.fasta>: Fasta-formatted file of the reference sequence that was used as input to bcftools_filter_and_align.pl. 

* \<table.txt>: table.txt file output by bcftools_filter_and_align.pl.

#### 4.5.2 Outputs

The program will output STDOUT a string of comma-separated numbers representing the number of invariant sites in the alignment corresponding to 'A,C,G,T'.

### 4.6 bcftools_filter_consensus_to_table.pl

Using the output(s) of bcftools_filter.pl from one or more separate alignments will Generate a table of variant positions in the same format as the table.txt file produced by bcftools_filter_and_align.pl as well as a thinned multiple genome alignment similar to alignment.fasta.

#### 4.6.1 Required inputs

* `-r`: Fasta sequence of the reference genome
* `-f`: Tab-separated file of consensus genome seuquences ouptut by bcftools_filter.pl from alignments against the reference sequence.
Format:
```
        /path/to/consensus.fasta <tab> genome_name
```

#### 4.6.2 Optional inputs

Optional:
* `-i`: Include reference sequence in the output. Default: reference sequence will not be included.
* `-p`: Output file prefix. Default: 'output'

#### 4.6.3 Outputs

Will produce two files \<prefix>.table.txt and \<prefix>.alignment.fasta. Format of these files will be similar to the output of bcftools_filter_and_align.pl described in section [4.2.4](#424-outputs). The only exception is that the read count values in the table file will all be "?" characters. 

### 4.7 blast_masker.pl

Performs blast of sequences against themselves and outputs regions of self-alignment
in the same interval format as NCBI's dustmaker. The output file can be used as input for the `-m` option of bcftools_filter.pl or bcftools_filter_and_align.pl.

Requires **blastn** to be in your PATH.

#### 4.7.1 Required inputs

* `-f`: Sequence file, in fasta format

#### 4.7.2 Optional inputs

* `-s`: Size of fragments to break the genome up into when blasting against
        itself, in bp. Enter '0' to not fragment the genome. Default: 500
* `-i`: Minimum percent sequence identity, in %. Default: 75
* `-o`: Output prefix. File will be output as: '\<prefix>_blastmasker\<fragment size>.txt'. Default: prefix will be derived from the input filename up to the first '.' character.
* `-d`  Dustmaker interval file, for comparison. Default: no comparison performed.
* `-t`  Threads for parallel processing. Default: 8.

#### 4.7.3 Outputs

Will produce a file of intervals for each sequence in the reference corresponding to regions with blast overlaps. Example format:

```
>sequence_ID
209 - 215
415 - 421
487 - 494
1104 - 1110
```

## 5. LICENSE:

**MIT License**

Copyright (c) 2021-2025 Egon A. Ozer

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## 6. CONTACT:

Contact [Egon Ozer](e-ozer@northwestern.edu) with questions or comments.

 