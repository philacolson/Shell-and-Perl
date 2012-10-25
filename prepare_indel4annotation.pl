#! /usr/bin/perl 


## this program merges the indels that have the same allele length, same type and overlapping location into a single site
## The input file will be list of file names that have the left-extreme and right-extreme of indels at the 3rd and 4th column
## The program will look into merge the indel of the current file and then across the file for indel mergeing
## Input file format looks like this:
## 17      65014141        65014138        65014143        NODE_14897      Good    deletion        155     152     157     48      CTT     ATGTCTCCTT(CTT>---)GTTTTTTTTT   HighQ:13;HQU:8;PassQ:5;LowQ:0;LQU:0
## 17      65014141        65014138        65014143        NODE_7814       Good    deletion        90      87      92      43      CTT     ATGTCTCCTT(CTT>---)GTTTTTTTTT   HighQ:7;HQU:5;PassQ:3;LowQ:0;LQU:0
## for replicate, the last field HighQ:7, the higher count of HighQ will be used as the line to represent the data

#=============================================
#NOTE NOTE NOTE NOTE NOTE
#strict prevents implicit data definition (eg by typo)
use strict;
#NOTE NOTE NOTE NOTE NOTE
#=============================================

use Getopt::Std;
use Carp;
use File::Basename;

#use vars necessary for getopt if you are using strict
use vars qw/$opt_i $opt_o $opt_s/;
getopt('i:o:s');

#multiline strings OK in perl
my $USAGE="
$0 takes the following arguement:
-i input file with the list of the files to be compared and their tag
-o Output file name for validate_somatic_loci.txt
-s Output file name for sub_sum.txt
";

croak "$USAGE" unless defined $opt_i;
croak "$USAGE" unless defined $opt_o;
croak "$USAGE" unless defined $opt_s;




my @combined_line_list;

open(FN, $opt_i) or die "Fail to open input file with the primer alignment:$opt_i\n";
open(FO, "> $opt_o") or die "Cannot open the output file:$opt_o\n";
open(FS, "> $opt_s") or die "Cannot open the output file:$opt_s\n";
printf FS "Gene	Sample	Repeat	Chr	Pos	#Unique	#Indel	#Total	#IndelN	#TotalN	RefAllele	MutantAllele	Flanking\n";
my $is_first_file=1;
while (<FN>)
{
    chomp;
#   break array into fields, using dummy vars for unused fields 
    if($is_first_file == 0)
    {
      my ($gene, $sample, $repeat, $chr, $pos, $unique, $indelCnt, $TotalCnt, $IndelNCnt, $TotalNCnt, $Flanking) = split "\t", $_;
      my @items = split /\(|>|\)/, $Flanking;
      my $allele1=$items[1];
      my $allele2=$items[2];
      my $allele_len=length($allele1);
      my $is_deletion=0;
      if($allele2 =~ "\-")
      {
        $is_deletion=1;
        $allele2="-";
      }
      else
      {
        $allele1="-";
      }
      my $from=$pos;
      my $to=$pos;
      if($is_deletion == 1)
      {
        $to=$pos + $allele_len -1;
      }
      else
      {
        $from -=1;
      }
      my $ref_allele=$allele1;
      my $other_allele=$allele2;
      printf FO "SJHQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $gene, $chr, $from, $to, $ref_allele, $ref_allele, $other_allele, $sample;
##      printf FO "SJHQ\t%s\tchr%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $gene, $chr, $from, $to, $ref_allele, $ref_allele, $other_allele, $sample;
      printf FS "SJHQ\t%s\t%s\t%s\tchr%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $gene, $sample, $repeat, $chr, $pos, $unique, $indelCnt, $TotalCnt, $IndelNCnt, $TotalNCnt, $ref_allele, $other_allele, $Flanking
    }
    $is_first_file=0;
}

close (FN);
close (FO);
close (FS);
