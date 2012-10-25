#! /usr/bin/perl

## this program extract the intervals from the database file that intersects with 
## the input data file
## this file assumes that the first two columes of the database files are the 
## start and stop positions
## the database file is expected to be sorted
use Getopt::Std;
use Carp;
use File::Basename;
getopt('i:d:o:b:e:q:m:n:a');

my $USAGE="$0 takes the following arguement:\n";
$USAGE .="-i input data file \n";
$USAGE .="-d the database file\n";
$USAGE .="-o output file name\n";
$USAGE .="-b buffer_intervals default=0\n";
$USAGE .="-e error log file\n";
$USAGE .="-q print out the query line. 0=FALSE 1=TURE. default=0\n";
$USAGE .="-m the index for the start position of the -i file. default=0\n";
$USAGE .="-n the index for the start position of the -d file. default=0\n";
$USAGE .="-a record all the DB interval file. default=0\n";


if(!defined $opt_i ){
  croak $USAGE;
}

if(!defined $opt_d){
  croak $USAGE;
}


if(!defined $opt_o){
  croak $USAGE;
}

if(!defined $opt_e){
  croak $USAGE;
}

if(!defined $opt_q){
  $PRINT_QUERY=0;
}
else
{
  $PRINT_QUERY=$opt_q;
}

$BUF_SIZE = 0;
if(defined $opt_b)
{
  $BUF_SIZE = $opt_b;
}

## loading the database file
if(defined $opt_n)
{
  $INDEX_START=$opt_n;
  $INDEX_STOP=$INDEX_START+1;
}
else
{
  $INDEX_START=0;
  $INDEX_STOP=1;
}

undef my @start_list;  
undef my @stop_list;
undef my @all_lines;

open(FN, $opt_d) or die "Cannot open database file :$opt_d\n";
$line = 0;
while(<FN>) {
  chomp;
  $all_lines[$line] = $_;
  @items=split(/\s/, $_);
  $count = scalar(@items);
  if($count >= 2)
  {
    $start_list[$line] = $items[$INDEX_START];
    $stop_list[$line] = $items[$INDEX_STOP];
    if($BUF_SIZE > 0)
    {
      $start_list[$line] -= $BUF_SIZE;
      if($start_list[$line] < 0)
      {
        $start_list[$line] = 0;
      }
      $stop_list[$line] += $BUF_SIZE;
    }
    $line += 1;
  }
}
close(FN);

## printf "total_line = $line, start = %ld, stop = %ld\n", $start_list[$line-1], $stop_list[$line-1];
use integer;
## a function for doing binary search
sub bin_search_interval 
{
  my ($low, $high, $start, $stop) = @_;
  my $mid;
  my $i_start, $i_stop;

  if($low > $high)
  {
    return -1;
  }
  $mid = ($low + $high)/2;
  $i_start = $start_list[$mid];
  $i_stop = $stop_list[$mid];
  if(!($start > $i_stop || $stop < $i_start))
  {
    return $mid;
  }
  else
  {
     if($stop < $i_start)
     {
       return bin_search_interval ($low, $mid-1, $start, $stop);
     }
     if($start > $i_stop)
     {
       return  bin_search_interval ($mid+1, $high, $start, $stop);
     }
     else
     {
       return -1;
     }
  }
}

undef my @good_lines;
undef my @query_index;
$good_line_count = 0;

open (FO,"> $opt_e") || die "Cannot open error log file :$opt_e\n";
## process the query file
if(defined $opt_m)
{
  $INDEX_START=$opt_m;
  $INDEX_STOP=$INDEX_START+1;
}
else
{
  $INDEX_START=0;
  $INDEX_STOP=1;
}
open(FN, $opt_i) or die "Cannot open input file :$opt_i\n";
while(<FN>) {
  chomp;
  $line_val = $_;
  @items=split(/\s/, $_);
  $count = scalar(@items);
  if($count < 2)
  {
    printf STDERR "Invalid value at line $line\n";
    exit 1;
  }
  $start = $items[$INDEX_START];
  $stop = $items[$INDEX_STOP];
  $index = bin_search_interval (0, $line-1, $start, $stop);
  if($index == -1)
  {
    printf FO "Fail to find interval for %ld, %ld\n",  $start, $stop;
  }
  else
  {
     $good_lines[$good_line_count] = $index;
     $query_index[$good_line_count] = $line_val;
    
     $good_line_count += 1;

## this option allows recording of all database interval
     if(defined $opt_a && $opt_a == 1)
     {
       for($i = $index - 1; $i >=0 && $start <= $stop_list[$i]; --$i)
       {
          $good_lines[$good_line_count] = $i;
          $query_index[$good_line_count] = $line_val;
          $good_line_count += 1;
       }
       for($i = $index +1; $i <$line && $stop >= $start_list[$i]; ++$i)
       {
          $good_lines[$good_line_count] = $i;
          $query_index[$good_line_count] = $line_val;
          $good_line_count += 1;
       }
     }
  }
}
close (FN);
close (FO);



open (FO,"> $opt_o") || die "Cannot open output file :$opt_o\n";
for($i = 0; $i < $good_line_count; ++$i)
{
  $line = $good_lines[$i];
  if($PRINT_QUERY == 1)
  {
     printf FO "%s\t", $query_index[$i];
  }
  printf FO "%s\n", $all_lines[$line];
}
close (FO);
