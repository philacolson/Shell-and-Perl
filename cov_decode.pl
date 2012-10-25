#!/usr/bin/perl -w
# MNE 9/24/09

use strict;

## use lib "/usr/local/cgap/lib/perl";
use lib $ENV{HOME} . "/bin";

use Getopt::Long;
use POSIX qw(SEEK_SET);

use constant NONSPECIFIC_ROOT => "/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess/";

use constant MM_MASK => 0x01;
use constant PM_MASK => 0x02;

my %FLAGS;
GetOptions(\%FLAGS,
	   "-chr=s",
	   "-size=i",

	   "-start=i",
	   # starting BASE NUMBER (not 0-based index)

	   "-length=i",
	  );

my $chr = $FLAGS{chr} || die "specify -chr\n";
# chromosome
my $size = $FLAGS{size} || die "specify -size\n";
# sequence size (?)

my $file = sprintf '%s/chr%s.all.%d', NONSPECIFIC_ROOT, $chr, $size;

my $length = $FLAGS{length} || -s $file;
my $start = $FLAGS{start} || 0;

die "$file does not exist\n" unless -s $file;
open(IN, $file) || die;
binmode(IN);

my $wanted = $length;
my $pos = 1;
$start-- if $start;
# covert from base number to offset
if ($start) {
  seek(IN, $start, SEEK_SET) || die "can't seek to initial position";
  $pos += $start;
}

my $buf;
my $v;
my @what;

my $BUF_SIZE = 4096;
my $i;
while ($wanted) {
  my $read = read(IN, $buf, $BUF_SIZE);
  if ($read > 0) {
    for ($i=0; $i < $read; $i++) {
      $v = unpack 'c', substr($buf, $i, 1);
      die unless $v;
      @what = ();
      push @what, "MM" if $v & MM_MASK;
      push @what, "PM" if $v & PM_MASK;
      printf "%s\n", join ",", $pos, @what;
      #    printf "%s\n", join ",", @what;
      $pos++;
      last unless --$wanted;
    }
  } else {
    # EOF
    last;
  }
}
