#!/usr/bin/perl
# Name:   fisherexact
# Author: Gang Wu
# Date:   Feb, 21, 2010
# This script takes tab-delimited id and 4 integer values, calculates the p-value 
# of two-tailed fisher exact. 
# Example: grep ^rs test.txt |perl fisherexact

use strict;
use warnings;
use Text::NSP::Measures::2D::Fisher::twotailed;
use Log::Log4perl qw(:easy);
# http://search.cpan.org/~tpederse/Text-NSP-1.11/lib/Text/NSP/Measures/2D/Fisher/twotailed.pm
Log::Log4perl->easy_init($ERROR);
##The parameter here is the level of logging desired, going from more verbose to less messages. $debug, $info, $warn, $error, $fatal
## the default contingency table is 2X2
my ($row,$col) = (2,2);
my $logger = get_logger();
while(<>){
  chomp;
  my $line=$_;
  my @values = split(/\t/);
#          word2   ~word2
#  word1    n11      n12 | n1p
# ~word1    n21      n22 | n2p
#           --------------
#           np1      np2   npp
  my $npp = $values[1] + $values[2] + $values[3] + $values[4]; 
  my $n1p = $values[1] + $values[2];
  my $np1 = $values[1] + $values[3];
  my $n11 = $values[1];

  my $twotailed_value = calculateStatistic( n11=>$n11,
                                      n1p=>$n1p,
                                      np1=>$np1,
                                      npp=>$npp);

  if( (my $errorCode = getErrorCode()))
  {
		$logger->warn("$values[0] had the error: $errorCode - " . getErrorMessage());
  }
  else
  {
##    print "$_\t";
##    printf "%.1e",$twotailed_value;
##    print "\n";
$logger->debug("$values[0] had a value of $twotailed_value.");
      printf "%lf\t%s\n", $twotailed_value, $line;
  }
}


