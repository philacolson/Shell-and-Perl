#! /usr/bin/perl
## This program modifies genotype2trace.txt file by inserting a new 5th column
## to the output that records allele_label for samples that have the minor allele
## it takes three input file: *_snp_info.txt to get the final SNP_type; 
## *_allele_rpt.txt to get the allele label
## genotype2trace.txt
## -m option will change the genotype call on the input genotype2trace.txt if it does not fit into 
## the alleles defined in *_snp_info.txt. If major_peak = wt, it will be called homozygous major; otherwise
## it will be called N>N
## it outputs a modified genotype2trace.txt file






use Getopt::Std;
use Carp;
use File::Basename;
getopt('i:j:k:m:o');

my $USAGE="$0 takes the following arguement:\n";
$USAGE .="-i input file with the *_snp_info.txt\n";
$USAGE .="-j input file with the *_allele_rpt.txt\n";
$USAGE .="-k input file with the *_genotype2trace.txt\n";
$USAGE .="-m modify genotypes that do not fit alleles in *_snp_info.txt. default=0\n";
$USAGE .="-o Output file name\n";




if(!defined $opt_i ){
  croak $USAGE;
}

if(!defined $opt_j ){
  croak $USAGE;
}

if(!defined $opt_j ){
  croak $USAGE;
}

my $change_inconsistent_genotype = 0;
if(defined $opt_m)
{
   $change_inconsistent_genotype = 1;
}


if(!defined $opt_o){
  croak $USAGE;
}

undef my %SNPAllele2label;

open(FN, $opt_j) or die "Fail to open input *_allele_rpt.txt file :$opt_j\n";
$line=0;
while (<FN>) {
  if($line >=1)
  {
  	chomp;
  	my @items = split "\t", $_;
        my $t_count=scalar(@items);
        my $SNP_name = $items[0];
        my $SNP_class = $items[2];
        my $allele_label = $items[3];
        my $allele = $items[4];
        $allele_label = $allele_label.";";
        if($t_count >=6)
        {
          $allele_label = $allele_label.$items[5];
        }
        $allele_label = $allele_label.";";
        if($t_count >=7)
        {
           $allele_label = $allele_label.$items[6];  
        }
        $SNPAllele2label{$SNP_name}{$allele}{$SNP_class} = $allele_label;
  }
  $line +=1;
}
close (FN);

undef my %SNP2class;
open(FN, $opt_i) or die "Fail to open input *_snp_info.txt file :$opt_i\n";
$line=0;
undef my %SNPallele_list;
while (<FN>) {
  if($line >=1)
  {
  	chomp;
  	my @items = split "\t", $_;
        my $SNP_name = $items[0];
        my $SNP_class = $items[4];
        $SNP2class{$SNP_name} = $SNP_class;
        my @allele_list = split /\(|>|\)/, $items[5];  ## CAGCCACAGG(C>T>G)TCCCAGACAT
        $SNPallele_list{$SNP_name}=\@allele_list;
   
  }
  $line +=1;
}
close(FN);

sub is_genotype_match {
  my ($A1, $A2, $B1, $B2) = @_;

  if(($A1 eq $B1) && ($A2 eq $B2))
  {
    return 1;
  }
  if(($A1 eq $B2) && ($A2 eq $B1))
  {
    return 1;
  }
  return 0;
}

sub get_min {
  my ($first, $second) = @_;

  if($first > $second)
  {
      return $second;
  }
  else
  {
       return $first;
  }
}




open(FN, $opt_k) or die "Fail to open input *_genotype2trace.txt file :$opt_k\n";
open(FO, "> $opt_o") or die "Cannot open the output file:$opt_o\n";
$line=0;
while (<FN>) {
  chomp;
  my @items = split "\t", $_;
  my $count = scalar(@items);
  my $skip_line = 0;
  if($line == 0)
  {
     for($i = 0; $i < 4; ++$i)
     {
        printf FO "%s\t", $items[$i];
     }
     print FO  "GenotypeLabel\tProteinGI\tmRNA_acc";
     for($i=4; $i<$count; ++$i)
     {
        printf FO "\t%s", $items[$i];
     }
     print FO "\n";
  }
  else
  {
     my $SNP_name = $items[0];
     my $load_allele_label = 0;
     my $allele1 = $items[2];
     my $allele2 = $items[3];
     my $allele_label = "N";
     my $t_count = scalar(@items);
     if(defined $SNP2class{$SNP_name})
     {
        my $SNP_class = $SNP2class{$SNP_name};
	if(defined $SNPAllele2label{$SNP_name}{$allele1}{$SNP_class})
        {
          $allele_label = $SNPAllele2label{$SNP_name}{$allele1}{$SNP_class};
          $load_allele_label = 1;
          $allele_label =~ s/;/\t/g
        }
	elsif(defined $SNPAllele2label{$SNP_name}{$allele2}{$SNP_class})
        {
          $allele_label = $SNPAllele2label{$SNP_name}{$allele2}{$SNP_class};
          $load_allele_label = 1;
          $allele_label =~ s/;/\t/g
        }
     }
     
      

     if($change_inconsistent_genotype != 0 )
     {
        if(defined $SNPallele_list{$SNP_name})
        {
          my $allele_list_ref=$SNPallele_list{$SNP_name};
          my @t_allele_list = @$allele_list_ref;  ## 0=left_flank, 1=allele1, 2=allel2, ...n=right_flank
          if(($allele1 ne $t_allele_list[1]) || ($allele2 ne $t_allele_list[1]))  ## this is not wild-type genotype
          {
            my $allele_1_match = -1;
            my $allele_2_match = -1;
            my $allele_count=scalar (@t_allele_list) - 1;
            for (my $k=1; $k <$allele_count; ++$k)
            {
              if($allele1 eq $t_allele_list[$k])
              {
                $allele_1_match = $k;
              }
              if($allele2 eq $t_allele_list[$k])
              {
                $allele_2_match = $k;
              }
            }
            if($allele_1_match == -1 || $allele_2_match == -1)
            {
              if($allele_1_match == 1)   ## first allele matches the wild-type
              {
                if($t_count < 21)  ## only have one orientation
                {
                   $items[3] = $items[2];
                }
                elsif(($items[11] eq $items[12]) || ($items[19] eq $items[20]))  ## one of the allele is a homozygous
                {
                   $items[3] = $items[2];
                }
                elsif(is_genotype_match ($items[11],  $items[12], $items[19], $items[20]) == 0)  ## genotype do not match
                {
                  $items[3] = $items[2];
                }
              }
              elsif($t_count < 21 || (is_genotype_match ($items[10], $items[11], $items[18], $items[19]) == 0) || (($items[10] ne $items[11]) && (get_min($items[20], $items[12]) < 35)))
              {
                $skip_line = 1;
                ## $items[2] = "N";
                ## $items[3] = "N";
                ## $items[4] = "0";
              }
              else
              {
                 $items[2]=$items[10];
                 $items[3]=$items[11];
              }
            }
          }
        }
     }

     if($skip_line == 0)
     {
       for($i = 0; $i < 4; ++$i)
       {
          printf FO "%s\t", $items[$i];
       }
       if($load_allele_label ==1)
       {
     	  print FO  "$allele_label";
       }
       for($i=4; $i<$count; ++$i)
       {
          printf FO "\t%s", $items[$i];
       }
       print FO "\n";
    }
  }
  $line +=1;
}

close (FP);
close (FO);

