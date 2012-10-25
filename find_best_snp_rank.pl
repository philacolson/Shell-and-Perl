#! /usr/bin/perl
## This program goes through the output of SNP2AA.out
## for each SNP with multiple status, pick the highest (or the most damaging)
## protein change


use Getopt::Std;
use Carp;
use File::Basename;
getopt('i:o');

my $USAGE="$0 takes the following arguement:\n";
$USAGE .="-i input file with the SNP2AA.out\n";
$USAGE .="-o Output file name\n";




if(!defined $opt_i ){
  croak $USAGE;
}

if(!defined $opt_o){
  croak $USAGE;
}


undef my %SNP_rank_list;

$SNP_rank_list{"nonsense"} = 1;
$SNP_rank_list{"missense"} = 2;
$SNP_rank_list{"silent"} = 4;
$SNP_rank_list{"splice"} = 5;
$SNP_rank_list{"frameshift"} = 6;
$SNP_rank_list{"UTR_5"} = 7;
$SNP_rank_list{"UTR_3"} = 8;
$SNP_rank_list{"exon"} = 9;
$SNP_rank_list{"intron"} = 10;
$SNP_rank_list{"unknown"} = 11;
$SNP_rank_list{"proteinIns"} = 3;
$SNP_rank_list{"proteinDel"} = 3;

undef my %SNP_name_class;
open(FN, $opt_i) or die "Fail to open input file :$opt_i\n";
open(FO, "> $opt_o") or die "Cannot open the output file:$opt_o\n";

$SNP_name_index = 0;
$SNP_class_index = 4;
$line=0;
while (<FN>) {
  if($line >=1)
  {
  chomp;
  my @items = split "\t", $_;
  $SNP_name = $items[$SNP_name_index];
  $class = $items[$SNP_class_index];

  if (!defined $SNP_name_class{$SNP_name})
  {
    $SNP_name_class{$SNP_name} = $class;
  }
  else
  {
    $c_rank = $SNP_rank_list{$class};
    $old_class = $SNP_name_class{$SNP_name};
    $o_rank = $SNP_rank_list{$old_class};
    if($c_rank < $o_rank)
    {
      $SNP_name_class{$SNP_name} = $class;
    }
   }
  }
  ++$line;
}
close (FN);

foreach $SNP_name (keys %SNP_name_class)
{
   printf FO "%s\t%s\n", $SNP_name, $SNP_name_class{$SNP_name};
}
close (FO);

