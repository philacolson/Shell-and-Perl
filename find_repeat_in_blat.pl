#! /usr/bin/perl
## This program goes through the BLAT output and 
## for each query evaluate whether it is  mapped uniquely
## the difficult thing is that BLAT score was not provided here
## calculate tehe score by requiring match -mismatch -num_gap
## overlapping regions
## only check the sites that meet the target goal
## set the minimum distance to score 5
## presence of good.out file determines that it passes the check. bad.out is simply for debugging purpose
## some records may produce neither good nor bad output because they may not include any data that generates mismatch in the targeted region


## avoid penalizing the real gene alignment. Needs to add in an input file that extract the exon-boundaries of the targetted genes
## pseudo-gene is considered a better match only when the identity of the real gene is lower than the pseudogene
## 100     1       0       0       0       0       1       3500    +       HWUSI-EAS1680_8:1:41:4614:17208_5021    101     0       101     chr7    158821424       140124312       140127913       2       32,69,  0,32,   140124312,140127844,
## 98      3       0       0       0       0       0       0       -       HWUSI-EAS1680_8:1:41:4614:17208_5021    101     0       101     chrX    154913754       74720716        74720817        1       101,    0,      74720716,

##
## BRAF    1       NM_004333       chr7    140080282       140081039       758     UCSC
## BRAF    2       NM_004333       chr7    140086081       140086215       135     UCSC
## BRAF    3       NM_004333       chr7    140095556       140095687       132     UCSC
## BRAF    4       NM_004333       chr7    140099544       140099662       119     UCSC
## BRAF    5       NM_004333       chr7    140100456       140100502       47      UCSC
## BRAF    6       NM_004333       chr7    140123181       140123357       177     UCSC
## BRAF    7       NM_004333       chr7    140124260       140124344*       85      UCSC
## BRAF    8       NM_004333       chr7    *140127845       140127962       118     UCSC
## BRAF    9       NM_004333       chr7    140129290       140129426       137     UCSC
## BRAF    10      NM_004333       chr7    140133817       140133853       37      UCSC
## BRAF    11      NM_004333       chr7    140140577       140140736       160     UCSC
## BRAF    12      NM_004333       chr7    140146631       140146750       120     UCSC
## BRAF    13      NM_004333       chr7    140147681       140147829       149     UCSC
## BRAF    14      NM_004333       chr7    140154229       140154331       103     UCSC
## BRAF    15      NM_004333       chr7    140155161       140155264       104     UCSC
## BRAF    16      NM_004333       chr7    140180878       140181141       264     UCSC
## BRAF    17      NM_004333       chr7    140196380       140196481       102     UCSC
## BRAF    18      NM_004333       chr7    140270835       140271033       199     UCSC

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
use vars qw/$opt_i $opt_d $opt_c $opt_p $opt_b $opt_o $opt_m $opt_s $opt_g $opt_e/;
getopt('i:d:c:p:b:o:m:s:g:e');

my $USAGE="$0 takes the following arguement:\n";
$USAGE .="-i input data file produced by BLAT\n";
$USAGE .="-d minimum distance between the first and the second hit. default set to 5\n";
$USAGE .="-c targeted chromosome name. use chr1, etc\n";
$USAGE .="-p targeted chromosome position. eg 25271546, etc\n";
$USAGE .="-b input file with exon boundaries.\n";
$USAGE .="-o Output file name for recording non-repeat\n";
$USAGE .="-m Lower distance for reads with 1-bp mismatch: default=1\n";
$USAGE .="-s Lower distance for reads mapped to exon boudnary: default=1\n";
$USAGE .="-g Number of good reads required: default=2\n";
$USAGE .="-e Error file name for recording replicate hit\n";

if(!defined $opt_i ){
  croak $USAGE;
}

if(!defined $opt_o){
  croak $USAGE;
}

my $MIN_QUERY_DISTANCE=5;
if(defined $opt_d) {
  $MIN_QUERY_DISTANCE=$opt_d;
}

if(!defined $opt_c){
  croak $USAGE;
}
my $TARGET_CHR=$opt_c;


if(!defined $opt_p){
  croak $USAGE;
}
my $TARGET_POS=$opt_p;

if(!defined $opt_e) {
  croak $USAGE;
}

my $CK_SINGLE_MISMATCH=1;
if (defined $opt_m)
{
   $CK_SINGLE_MISMATCH=$opt_m;
}

my $CK_EXON_MATCH=1;
if (defined $opt_s)
{
  $CK_EXON_MATCH=$opt_s;
}

my $MIN_GOOD_READ=2;
if(defined $opt_g)
{
  $MIN_GOOD_READ=$opt_g;
}


undef my %exon_interval_list;

## BRAF    16      NM_004333       chr7    140180878       140181141       264     UCSC
## BRAF    17      NM_004333       chr7    140196380       140196481       102     UCSC
## BRAF    18      NM_004333       chr7    140270835       140271033       199     UCSC
my $EXON_INDEX_START=4;
my $EXON_INDEX_STOP=5;
my $BUF_SIZE=5;  ## buffer size of start and stop exons
my $MAX_QUERY_UNALIGNED_SIZE=30;   ##allow maximum of 30bp of unaligned to determine whether an unaligned segment is considered qualified for matching exon-boundary. this threshold is needed to differentiate unaligned segment due to poor alignment or due to shorter read bases across the adjacent exons

undef my @start_list;
undef my @stop_list;
my $exon_line = 0;

if (defined $opt_b)
{
  open(FN, $opt_b) or die "Fail to open the input data file from BLAT output:$opt_b\n";
  while(<FN>)
  {
    chomp;
    my @items=split(/\s/, $_);
    my $count = scalar(@items);
    if($count >= $EXON_INDEX_STOP+1)
    {
      $start_list[$exon_line] = $items[$EXON_INDEX_START];
      $stop_list[$exon_line] = $items[$EXON_INDEX_STOP];
      if($start_list[$exon_line] < 0)
      {
        $start_list[$exon_line] = 0;
      }
      $exon_line += 1;
    }
  }
  close(FN);
}





## psLayout version for BLAT output

## match|mis- |rep. |N's|Q gap|Q gap|T gap|T gap|strand|Q        |Q   |Q    |Q  |T        |T   |T    |T  |block|blockSizes |qStarts| tStarts
##      |match|match|   |count|bases|count|bases|      |name     |size|start|end|name     |size|start|end|count  

my $QUERY_INDEX=9;
my $MISMATCH_INDEX=1;
my $MATCH_INDEX=0;
my $QUERY_GAP_INDEX=4;
my $CHR_GAP_INDEX=6;
my $QUERY_SIZE_INDEX=10;
my $QUERY_START_INDEX=11;
my $QUERY_STOP_INDEX=12;
my $CHR_START_INDEX=15;
my $CHR_STOP_INDEX=16;
my $CHR_INDEX=13;
my $CHR_SEG_SIZE_INDEX=18;
my $CHR_SEG_START_INDEX=20;
my $QUERY_GAP_SIZE_INDEX=5;
my $CHR_GAP_SIZE_INDEX=7;


## only check the queries that fit into the categoriy
undef my %check_query_list;
my $index=0;
undef my @line_info_list;

open(FN, $opt_i) or die "Fail to open the input data file from BLAT output:$opt_i\n";
open(FO, "> $opt_o") or die "Cannot open the output file:$opt_o\n";
open(FE, "> $opt_e") or die "Cannot open error output file:$opt_e\n";
my $index = 0;
while (<FN>) {
  chomp;
  my $line = $_;
  my @items = split "\t", $_;
  undef my %line_info;
  $line_info{query} = $items[$QUERY_INDEX];
  $line_info{mismatch} = $items[$MISMATCH_INDEX];
  $line_info{match} = $items[$MATCH_INDEX];
  $line_info{query_gap} = $items[$QUERY_GAP_INDEX];
  $line_info{chr_gap}=$items[$CHR_GAP_INDEX];
  $line_info{match} = $items[$MATCH_INDEX];
  $line_info{qstart} = $items[$QUERY_START_INDEX];
  $line_info{qstop} = $items[$QUERY_STOP_INDEX];
  $line_info{size} = $line_info{qstop} - $line_info{qstart} + 1;
  $line_info{query_size} = $items[$QUERY_SIZE_INDEX];
  $line_info{chrstart} = $items[$CHR_START_INDEX];
  $line_info{chrstop} = $items[$CHR_STOP_INDEX];
  $line_info{chr} = $items[$CHR_INDEX];
  $line_info{identity} = $line_info{match}/($line_info{match} + $line_info{mismatch});	## ignore any difference <1/1000
  $line_info{takeit} = 0;  ## not taken yet
  $line_info{seg_size_list}=$items[$CHR_SEG_SIZE_INDEX];
  $line_info{seg_start_list}=$items[$CHR_SEG_START_INDEX];
  $line_info{query_gap_size}=$items[$QUERY_GAP_SIZE_INDEX];
  $line_info{chr_gap_size}=$items[$CHR_GAP_SIZE_INDEX];  ## this is for deletions in cDNA and genomic alignment
  $line_info{line} = $line;
  $line_info{replicate} = 0;
  $line_info_list[$index] = \%line_info;
  if($items[$CHR_INDEX] eq $TARGET_CHR)  ## some records may fail in this test and generate neither good nor bad.out file
  {
    if($TARGET_POS >=$items[$CHR_START_INDEX] && $TARGET_POS <=$items[$CHR_STOP_INDEX] && $items[$MISMATCH_INDEX] >=1)
    {
      my $query_name=$items[$QUERY_INDEX];
      $check_query_list{$query_name} = $index;
      $line_info{target} = 1;
      if($line_info{mismatch} == 1 && ($line_info{match}+$line_info{mismatch} == $line_info{query_size}) && $line_info{query_gap} == 0 && $line_info{chr_gap} == 0)
      {
        $line_info{single_mismatch}=1;
      }
    }
  }
  ++$index;
}
close (FN);


sub is_alignment_match_exon{
  my ($line_info_ref) = @_;
  my $seg_size_val=$line_info_ref->{seg_size_list};
  my $seg_start_list=$line_info_ref->{seg_start_list};
  my $match_exon_boundary = 0;

  if($exon_line < 1)
  {
    return 0;
  }
  my @chr_start_list=split "\,", $seg_start_list;
  my @chr_size_list=split "\,", $seg_size_val;
  my $t_num=scalar (@chr_size_list);

  for(my $i =0; $i <$t_num; ++$i)
  {
     my $seg_start=$chr_start_list[$i];
     my $seg_stop=$seg_start + $chr_size_list[$i] -1;
     for(my $j=0; $j <$exon_line; ++$j)
     {
        my $exon_start = $start_list[$j];
        my $exon_stop = $stop_list[$j];
        if($seg_start >=$exon_start -$BUF_SIZE && $seg_start <=$exon_start + $BUF_SIZE)
        {
          if($match_exon_boundary <3)  ## 3 means matching the start of an exon
          {
            $match_exon_boundary +=1;
          }
	}
        if($seg_stop >=$exon_stop -$BUF_SIZE && $seg_stop <=$exon_stop + $BUF_SIZE)
        {
          if($match_exon_boundary <3)
          {
            $match_exon_boundary+=2;
          }
	}
     }
  }
  return $match_exon_boundary;
}


## go through each of these queries to evaluate the 
## matching exon boundary means no need to penalize chr_gap
my $good_read=0;
my %good_query_list;
foreach my $query_name (keys %check_query_list)
{
  my $matching_index=$check_query_list{$query_name};
  my $matching_hit_info=$line_info_list[$matching_index];
  my $number_replicate=0;
  my $alignment_match_exon = 0;
  my $match_score = 0;
  my $check_score = 0;
  my $score_diff= $MIN_QUERY_DISTANCE;
  my $check_score_diff=0;

  $alignment_match_exon = is_alignment_match_exon ($matching_hit_info);
  if($alignment_match_exon == 1 || $alignment_match_exon == 2)  ## disqualify an un-aligned segment from being matching_exon if the unaligned size is expected to be alignable to the neighboring exon
  {
     my $unaligned_query_size=$matching_hit_info->{query_size} - $matching_hit_info->{size};
     if( $unaligned_query_size > $MAX_QUERY_UNALIGNED_SIZE || $matching_hit_info->{chr_gap} !=0)
     {
       $alignment_match_exon = 0;
     }
  }
  if($alignment_match_exon ==0)
  {
     $match_score = $matching_hit_info->{match} - $matching_hit_info->{mismatch} - $matching_hit_info->{chr_gap} - $matching_hit_info->{query_gap};
     if($matching_hit_info->{single_mismatch}==1 && $CK_SINGLE_MISMATCH ==1)  ## special case for those with only 1-bp mismatch
     {
       $score_diff=0;
     }
  }
  else  ## should consider add any bonus scores for matching??
  {
     if($CK_EXON_MATCH ==1)
     {
       $match_score = $matching_hit_info->{match} - $matching_hit_info->{mismatch};
       $score_diff=0;  ## for those that match the exon boudnaries, only take the alternative hit when the second location was higher
     }
  }
     
  for(my $i =0; $i <$index; ++$i)
  {
    if($i != $matching_index)
    {
      $check_score_diff=$score_diff;
      my $other_hit_info=$line_info_list[$i];
      if($other_hit_info->{query} eq $query_name)
      {
        if($query_name eq "HWUSI-EAS1680_8:2:49:1117:11506_186")
        {
          print "stop here\n";
        }
        $check_score = $other_hit_info->{match} - $other_hit_info->{mismatch} - $other_hit_info->{chr_gap} -$other_hit_info->{query_gap};
        $check_score_diff=$score_diff;
        if($alignment_match_exon != 0)
        {
          if($other_hit_info->{chr_gap_size} >100 || ($other_hit_info->{chr_gap_size} >20 && $matching_hit_info->{chr_gap_size}/$other_hit_info->{chr_gap_size} <10))
          {  ## hit may also conform to intron/exon boundary 
             $check_score = $other_hit_info->{match} - $other_hit_info->{mismatch};
             $check_score_diff=$MIN_QUERY_DISTANCE;
          }
        }
        if(($alignment_match_exon != 1 && $alignment_match_exon !=2) || $other_hit_info->{chr_gap_size} >100)  ## 1 or 2 = only one-end of the exons are matched
        {
          if($match_score - $check_score <= $check_score_diff)
          {
             $number_replicate = +1;
             $other_hit_info->{replicate} = 1;
          }
        }
        else  ## maybe too lenient
        {
           if($other_hit_info->{mismatch} + $other_hit_info->{query_gap} +$other_hit_info->{chr_gap}<=$matching_hit_info->{mismatch})
           {
              $number_replicate = +1;
              $other_hit_info->{replicate} = 1;
           }
        }
      }
    }
  }
  if($number_replicate ==0)
  {
    $good_read +=1;
    $good_query_list{$query_name} = $check_query_list{$query_name};
  }
  else
  {
     printf FE "%s\n", $matching_hit_info->{line};
     for(my $i =0; $i <$index; ++$i)
     {
       if($i != $matching_index)
       {
         my $other_hit_info=$line_info_list[$i];
         if($other_hit_info->{query} eq $query_name)
         {
           if($other_hit_info->{replicate} == 1)
           {
             printf FE "#%s\n", $other_hit_info->{line};
           }
           else
           {
             printf FE "*%s\n", $other_hit_info->{line};
           }
         }
       }
     }
   }
}
close (FE);

## only print out the good results when there are more than one good read
if ($good_read >= $MIN_GOOD_READ)
{
## may want to check the query start position to ensure that only multiple, non-overlapping reads are covered
  foreach my $query_name (keys %good_query_list)
  {
    my $matching_index=$good_query_list{$query_name};
    my $matching_hit_info=$line_info_list[$matching_index];
     printf FO "%s\n", $matching_hit_info->{line};
     for(my $i =0; $i <$index; ++$i)
     {
       if($i != $matching_index)
       {
         my $other_hit_info=$line_info_list[$i];
         if($other_hit_info->{query} eq $query_name)
         {
           printf FO "*%s\n", $other_hit_info->{line};
         }
       }
     }
   }
}
close (FO);
