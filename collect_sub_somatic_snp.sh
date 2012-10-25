#!/bin/sh
INPUT_FILE_ORIGINAL=$1
PATIENT_ID=$2   #it will produce patientID_low_confidence_somatic_snp.txt and patientID_high_confidence_snp.txt
OUTPUT_RS_FALSE_FILE=$3  ## false positive rs counts
FILTER_RS_THRESHOLD=$4  ## how much percent of rs# is allowed to be there? 95 or 98 or 99
GERMLINE_SNP_FILE=$5
LOG_FILE=$6
SINGLE_STRAND_FILE=$7   ## information about the single-strand orientation

NEXTGEN_BIN_ROOT_DIR=/user/songliu/u2/group/Jianmin/snp_postprocess/snv_postprocess
PERLSRC_DIR=/user/songliu/u2/group/Jianmin/snp_postprocess/snv_postprocess
#BIN_DIR=${NEXTGEN_BIN_ROOT_DIR}/jzbin
FIND_SUB=${PERLSRC_DIR}/FindSub.pl
#PTR_FINDER=${BIN_DIR}/ptrfinder
PTR_FINDER=${NEXTGEN_BIN_ROOT_DIR}/ptrfinder
#SOMATIC_CRITER_LIST=${NEXTGEN_BIN_ROOT_DIR}/nextgensupport/SomaticCvgCritera.lst
SOMATIC_CRITER_LIST=${NEXTGEN_BIN_ROOT_DIR}/SomaticCvgCritera.lst
#CK_NON_SPECIFIC_PROGRAM=${NEXTGEN_BIN_ROOT_DIR}/scripts4all/check_replicate_region
CK_NON_SPECIFIC_PROGRAM=${NEXTGEN_BIN_ROOT_DIR}/check_replicate_region
FisherTest=${PERLSRC_DIR}/fisherexact


## need to do the export for running FisherTest
export PERL5LIB=${NEXTGEN_BIN_ROOT_DIR}/MikePerl/FishPerl/Text-NSP-1.11

## previous verison is collect_sub_somatic_snp_4_8_2010

if test -s ${INPUT_FILE_ORIGINAL}.gz
then
  gunzip ${INPUT_FILE_ORIGINAL}.gz
fi

if test ! -s ${INPUT_FILE_ORIGINAL}
then
  echo "Fail to find file ${INPUT_FILE_ORIGINAL}"
  exit 1
fi



if test -s ${LOG_FILE}
then
  rm ${LOG_FILE}
fi

## column definition in the original file
## col12      reference_normal_count
## col13      reference_tumor_count
## col14      alternative_normal_count 
## col15      alternative_tumor_count
## col16      dbSNP

#broad is perhaps a wider acceptance range? FIXME
## Mar. 29 2011. use the broad reference data which replaces the original one
## 20      broad_reference_normal_count
## 21      broad_reference_tumor_count FIXME
## 22      broad_alternative_normal_count
## 23      broad_alternative_tumor_count
## 24      dbSNP

INPUT_FILE=${INPUT_FILE_ORIGINAL}.switch  ## this file re-capitual the substitution variation format in the original file by using values in col 20-24 to replace those in the original col12-16

if test -s ${INPUT_FILE}
then
  rm ${INPUT_FILE}
fi

## Use broad count for tumor only
awk '{for(i=1; i<=11; ++i) printf("%s\t", $i); printf("%s\t%s\t%s\t%s", $20, $13, $22, $15); if(NF==24) printf("\t%s", $24); printf("\n")}' ${INPUT_FILE_ORIGINAL} >${INPUT_FILE}
## cut -f1-11,20-24 ${INPUT_FILE_ORIGINAL} >${INPUT_FILE}
## collect for the single-orientation variations:
## head -n1 SJINF001_bam_high_20.out |cut -f1,16-18,21
## Name    alternative_fwd_count   alternative_rev_count   alternative_bidirectional_confirmation  broad_reference_tumor_count
if test ! -s ${SINGLE_STRAND_FILE}
then
  grep SNP ${INPUT_FILE_ORIGINAL} |awk '{if($18 ==0) print $0}'  |cut -f1,16-18,21 >${INPUT_FILE_ORIGINAL}.single_strand.out
  cut -f1 ${INPUT_FILE_ORIGINAL}.single_strand.out |sort -u >${INPUT_FILE_ORIGINAL}.single_strand.lst
else
  cp ${SINGLE_STRAND_FILE} ${INPUT_FILE_ORIGINAL}.single_strand.lst
fi


## ${INPUT_FILE}.tmp
##Mutant alleles CAN occur in normal because they could be germline mutations
## This file generates the output in 6 colmn
## col1=snp_name (chr.loc)
## col2=% of mutant allele in TUMOR - % of mutant allele in NORMAL
## col3=% of mutant allele in NORMAL
## col4=% of mutant allele in TUMOR
## col5= # of all reads in NORMAL. This needs to be $12+$14 since it will combined the normal_reference+normal_alternative
## col6= # of mutant allele in TUMOR
##


if test ! -s ${SOMATIC_CRITER_LIST}
then
  echo "fail to find criteria file ${SOMATIC_CRITER_LIST}"
  exit 1
fi

#1) First print out the input file
#2) Filter and keep the bam reads that have mutations in either the reference OR the alternative because you cannot get the percentage of something with 0 occurences. 
#3) Print out the referenceID, % of alleles that are mutant in NORMAL, % of alleles that are mutant in TUMOR, # of all reads in NORMAL, and the number of mutant alleles in the tumor file.

cat ${INPUT_FILE} |awk '{if($14+$12!= 0 && $13+$15 != 0) print $0}' |awk '{printf("%s\t%ld\t%ld\t%s\t%s\n", $1, $14*100/($14+$12), $15*100/($13+$15), $12+$14, $15)}'  |awk '{printf("%s\t%ld\t%ld\t%ld\t%ld\t%ld\n", $1, $3-$2, $2, $3, $4, $5)}' >${INPUT_FILE}.tmp.original


## for now just filter out anything that has the normal in it. May consider something more flexible in the future that accomodate normal contamination with tumor tissue. Can run a few filter on the normal counts
## modify the file ${INPUT_FILE}.tmp.original by filtering out the germline SNPs found in normal
###PHRASING OF VARIABLES UNCLEAR:  original.clean has only the NON-GERMLINE snps.  Those are UNCLEAN.  They're mutations! FIXME
if test -s ${GERMLINE_SNP_FILE}
then
  grep SNP ${GERMLINE_SNP_FILE} |cut -f1 >${GERMLINE_SNP_FILE}.lst 
  if test -s ${GERMLINE_SNP_FILE}.lst
  then
    ${FIND_SUB} -i ${INPUT_FILE}.tmp.original -c ${GERMLINE_SNP_FILE}.lst -t 0 -d '\t' -n 0 -o  ${INPUT_FILE}.tmp.original.clean
    original_count=`wc ${INPUT_FILE}.tmp.original |awk '{printf("%s", $1)}'`
    filtered_count=`wc ${INPUT_FILE}.tmp.original.clean  |awk '{printf("%s", $1)}'`
    echo "original_count=$original_count	filtered_count=$filtered_count"
    echo "original_count=$original_count	germline_filtered_count=$filtered_count" >>${LOG_FILE}
#FILTERING AND RENAMING OCCURS HERE
    mv  ${INPUT_FILE}.tmp.original.clean  ${INPUT_FILE}.tmp.original
  fi
  rm  ${GERMLINE_SNP_FILE}.lst
fi

########################################
#NOW WE HAVE REMOVED ALL THE GERMLINE SNPS

if test ! -s ${INPUT_FILE}.tmp.original
then
  echo "Fail to generate ${INPUT_FILE}.tmp.original"
  exit 1
fi
##OPTIMIZE: Couldn't we strain the things before the FIND_SUB?  Carrying the extra snps that we aren't even using seems uselessly memory intensive.
## needs to eliminate cases where there are 0 reads in tumor normal normal
awk '{if($15+$13>0 && $14+$12 >0) printf("%s\t%ld\t%ld\t%ld\t%ld\n", $1, $15, $13, $14, $12)}' ${INPUT_FILE} |grep -v Name |${FisherTest} >fisher_output.txt

echo "high_confidence" >confidence.lst
echo "low_confidence" >>confidence.lst


HIGH_CCONFIDENCE_OUTPUT_FILE=${PATIENT_ID}_high_confidence_somatic_snp.txt
LOW_CCONFIDENCE_OUTPUT_FILE=${PATIENT_ID}_low_confidence_somatic_snp.txt
for confidence in `cat confidence.lst`; do
  echo $confidence >>${LOG_FILE}
  OUTPUT_FILE=${PATIENT_ID}_${confidence}_somatic_snp.txt
## doing both the high and low confidence analysis at the same time
if test -s ${OUTPUT_RS_FALSE_FILE}
then
  rm ${OUTPUT_RS_FALSE_FILE}
fi
touch ${OUTPUT_RS_FALSE_FILE}

if test -s ${OUTPUT_FILE}
then
  rm ${OUTPUT_FILE}
fi
touch ${OUTPUT_FILE}

REPEAT_LIST=${OUTPUT_FILE}.repeat.lst
if test -s ${REPEAT_LIST}
then
  rm ${REPEAT_LIST}
fi

if test -s ${REPEAT_LIST}.high
then
  rm ${REPEAT_LIST}.high
fi

  

if test -s snp_test.lst
then
  rm snp_test.lst
fi

if test -s fisher_snp_very_high.lst
then
  rm fisher_snp_very_high.lst
fi
## calculating p value using Fisher's exact test
awk '{if($1 <=0.001 && $3 >0 && $5 ==0) printf("%s\n", $2)}' fisher_output.txt  |sort -u >fisher_snp_very_high.lst

awk '{if($1 <=0.01) printf("%s\n", $2)}' fisher_output.txt  |sort -u >fisher_snp_high.lst
if test -s fisher_snp_high.lst
then
  echo "fisher_snp_high.lst" >snp_test.lst
fi
#Why are we ordering them and then putting them in the same file anyway? Is it because we want everything DOWN TO high quality?
if test ${confidence} = high_confidence
then
  awk '{if($1 <=0.05) printf("%s\n", $2)}' fisher_output.txt  |sort -u >fisher_snp_low.lst
  if test -s fisher_snp_low.lst
  then
    echo "fisher_snp_low.lst" >>snp_test.lst
  fi
else
  awk '{if($1 <=0.10) printf("%s\n", $2)}' fisher_output.txt  |sort -u >fisher_snp_low.lst
##  awk '{if($1 >0.10 && $1 <=0.25 && $3 >=3 && $6 >=10) printf("%s\n", $2)}' fisher_output.txt  |sort -u >fisher_snp_very_low.lst
  awk '{if($1 >0.10 && $1 <=0.25 && $3 >=3 && $6 >=5 && $5==0) printf("%s\n", $2)}' fisher_output.txt  |sort -u >fisher_snp_very_low.lst
  awk '{if($1 >0.25 && $1 <=0.50 && $3 >=3 && $5 == 0 && $6 >=8 && $3*100/($3+$4) >=20) printf("%s\n", $2)}' fisher_output.txt  |sort -u >>fisher_snp_very_low.lst
  sort -u fisher_snp_very_low.lst >fisher_snp_very_low.lst.sort 
  mv fisher_snp_very_low.lst.sort fisher_snp_very_low.lst #FIXME there has to be a more efficient way to do this
  if test -s fisher_snp_low.lst
  then
	##It has come to my attention that fisher snps differentiate themselves from other lists - fisher version merits its own file
    echo "fisher_snp_low.lst" >>snp_test.lst 
  fi

  if test -s fisher_snp_very_low.lst
  then
    cat fisher_snp_very_low.lst >>fisher_snp_low.lst
  fi
fi
#############################################
#We have now made a list of every snp that has a determined quality level or above.



#Whittles down the original input file for only the snps that we are concerned with.
for SIGNIFICANT_SNP_LIST in `cat snp_test.lst`; do
  echo ${SIGNIFICANT_SNP_LIST}
  ${FIND_SUB} -i ${INPUT_FILE}.tmp.original -c ${SIGNIFICANT_SNP_LIST} -t 1 -d '\t' -n 0 -o ${INPUT_FILE}.tmp

## check for tp53 mutations
  if test ${SIGNIFICANT_SNP_LIST} = fisher_snp_high.lst
  then
    awk '{if(index($2, "chr17") == 1) print $0}' ${INPUT_FILE} |awk '{if($3 >=7513394 && $3 <=7520826 && $14+$12!= 0 && $13+$15 != 0) print $0}' |awk '{printf("%s\t%ld\t%ld\t%s\t%s\n", $1, $14*100/($14+$12), $15*100/($13+$15), $12, $15)}'  |awk '{printf("%s\t%ld\t%ld\t%ld\t%ld\t%ld\n", $1, $3-$2, $2, $3, $4, $5)}' |awk '{if($2 > 20 && $3 <=3 && $6 >=3 && $5>=3) print $0}' |cut -f1 >TP53_snp.lst
    fgrep -f TP53_snp.lst ${INPUT_FILE} >>${OUTPUT_FILE}
    echo "TP53"
    cat ${OUTPUT_FILE}
  fi
##for the classification levels we care about, do
  for THRESHOLD in `cat ${SOMATIC_CRITER_LIST}`; do
    echo ${THRESHOLD}
## col2=% of mutant allele in TUMOR - % of mutant allele in NORMAL
	## col3=% of mutant allele in NORMAL
	## col4=% of mutant allele in TUMOR
	## col5= # of all reads in NORMAL. This needs to be $12+$14 since it will combined the normal_reference+normal_alternative
	## col6= # of mutant allele in TUMOR
## clean up the temp.out file which record those that meet the coverage criteria so we aren't dealing with data from an older run
    if test -s temp.out
    then
      rm temp.out
    fi

    if test ${THRESHOLD} = HIGH 
    then
	
	##GERMLINE COVERAGE IS HIGH: ($5 > 20)
	##TUMOR FREQUENCY IS BELOW 5%: ($3 <= 5) (tumor or normal?  Figure out parameters!)
	##MUTANT ALLELE COUNT IN TUMOR COMPARED TO GERMLINE: ($6/($5*$3/100)) >= [TIMES HIGHER]
	
## high-freqeuncy (frequency= high ratio of the snp in the tumor compared to the normal) criteria
## tumor-normal >=50%
## normal <=1%. 1% can only be achieved when you have a huge number of normal samples
##        cat ${INPUT_FILE}.tmp |awk '{if($2 >=50 && $3 <=1) print $0}' |awk '{if($5 >=4 && $6 >=5) print $0}' >temp.out
      if test ${confidence} = high_confidence
      then
        	if test ${SIGNIFICANT_SNP_LIST} = fisher_snp_high.lst
           	then
	#Further filter them based on stricter standards
	##If there is at least 50% more of a snp in the tumor cell than in the normal cell AND there is at most 3% in the normal - or 5% with a low total count.  Out of THOSE, just print out the ones that have a good number of reads
             cat ${INPUT_FILE}.tmp |awk '{if($2 >=50 && ($3 <=3 || ($3 <=5 && $5*$3/100 <=1))) print $0}' |awk '{if($5 >=10 && $6 >=4) print $0}' >temp.out
        	else
#fisher_snp_low: less tolerant of mutant allele in the normal
             cat ${INPUT_FILE}.tmp |awk '{if($2 >=50 && ($3 <=1 || ($3 <=5 && $5*$3/100 <=1))) print $0}' |awk '{if($5 >=10 && $6 >=4) print $0}' >temp.out
        	fi
	#initial filtering yes, but no filtering out the ones with low read counts or low mutant allele [tumor] counts
      else
        cat ${INPUT_FILE}.tmp |awk '{if($2 >=50 && ($3 <=3 || ($3 <=5 && $5*$3/100 <=1))) print $0}' >temp.out
      fi

    else
      if test ${THRESHOLD} = MED
      then
        if test ${SIGNIFICANT_SNP_LIST} = fisher_snp_high.lst
        then
          if test ${confidence} = high_confidence
          then
            cat ${INPUT_FILE}.tmp |awk '{if($2 >=25 && ($3 <=1 || ($3 <=5 && $5*$3/100 <=1)) && $5 >=10) print $0}' >temp.out
          else
### a bit more  lenient to tumor contamination in germline. If germline coverage is high ($5>20) and the normal (this used to say tumor, but that didn't make any sense PAC) frequency is below 5% ($3<=5) and mutant allele count in tumor is 10 times as high as germline $6/($5*$3/100) >=10

##So, essentially, if there's a big difference in occurence of a mutant gene in the cancer gene but not the germline gene, we're interested in it.
            cat ${INPUT_FILE}.tmp |awk '{if($2 >=25 && ($3 <=3 || ($3 <=5 && ($5*$3/100 <=1 || ($5 >=20 && $6/($5*$3/100) >=10))))) print $0}' >temp.out
##            cat ${INPUT_FILE}.tmp |awk '{if($2 >=25 && ($3 <=3 || ($3 <=5 && $5*$3/100 <=1))) print $0}' >temp.out
          fi
        else
## medium frequency info. 20%
          if test ${confidence} = high_confidence
          then
            cat ${INPUT_FILE}.tmp |awk '{if($2 >=40 && ($3 <=1 ||  ($3 <=5 && $5*$3/100 <=1))) print $0}' |awk '{if($5 >=15 && $6 >=4) print $0}' >temp.out
            cat ${INPUT_FILE}.tmp |awk '{if($2 >=30 && $2<40 && ($3 <=1 || ($3 <=5 && $5*$3/100 <=1))) print $0}' |awk '{if($5 >=18 && $6 >=4) print $0}' >>temp.out
            cat ${INPUT_FILE}.tmp |awk '{if($2 >=25 && $2<30 && ($3 <=1  ($3 <=5 && $5*$3/100 <=1))) print $0}' |awk '{if($5 >=20 && $6 >=4) print $0}' >>temp.out
##Are we ordering them and then putting them back in the same file?  Or are all these different qualities that make them high-confidence?
            cat ${INPUT_FILE}.tmp |awk '{if($2 >=25 && $3 == 0 && $5 >=15) print $0}' >>temp.out
          else
            cat ${INPUT_FILE}.tmp |awk '{if($2 >=25 && $3 <1) print $0}' >temp.out
            cat ${INPUT_FILE}.tmp |awk '{if($2 >=25 && $2<30 && ($3 <=1 ||  ($3 <=5 && $5*$3/100 <=1))) print $0}' |awk '{if($5 >=20 && $6 >=4) print $0}' >>temp.out
         fi
         sort -u temp.out >temp.out.sort
         mv temp.out.sort temp.out
##         cat ${INPUT_FILE}.tmp |awk '{if($2 >=25 && $2<50 && $3 <=1) print $0}' |awk '{if($5 >=20 && $6 >=4) print $0}' >>temp.out
        fi
      else
## low frequency info. >=10% <20%
        if test ${SIGNIFICANT_SNP_LIST} = fisher_snp_high.lst ##This if statement is useless- it does the same thing regardless.  See below.
        then
          if test ${confidence} = high_confidence
          then
            cat ${INPUT_FILE}.tmp |awk '{if($2 >=20 && $2<25 && $3 <1) print $0}' |awk '{if($5 >=20 && $6 >=4) print $0}' >temp.out##this code
            cat ${INPUT_FILE}.tmp |awk '{if($2 >=15 && $2<20 && $3 <1) print $0}' |awk '{if($5 >=25 && $6 >=4) print $0}' >>temp.out##this code
            cat ${INPUT_FILE}.tmp |awk '{if($2 >=10 && $2<15 && $3 <1) print $0}' |awk '{if($5 >=30 && $6 >=4) print $0}' >>temp.out##this code
          else
            cat ${INPUT_FILE}.tmp |awk '{if($2 >=15 && ($3 <=1 ||  ($3 <=5 && $5*$3/100 <=1))) print $0}' >temp.out
          fi
        else
          if test  ${confidence} = high_confidence
          then
            cat ${INPUT_FILE}.tmp |awk '{if($2 >=20 && $2<25 && $3 <1) print $0}' |awk '{if($5 >=20 && $6 >=4) print $0}' >temp.out ##is the same as this code
            cat ${INPUT_FILE}.tmp |awk '{if($2 >=15 && $2<20 && $3 <1) print $0}' |awk '{if($5 >=25 && $6 >=4) print $0}' >>temp.out ##is the same as this code #FIXME
            cat ${INPUT_FILE}.tmp |awk '{if($2 >=10 && $2<15 && $3 <1) print $0}' |awk '{if($5 >=30 && $6 >=4) print $0}' >>temp.out##is the same as this code
          else
            cat ${INPUT_FILE}.tmp |awk '{if($2 >=15 && ($3 <=1 ||  ($3 <=5 && $5*$3/100 <=1))) print $0}' >temp.out
          fi
        fi
      fi  ## end for low-frequency
    fi

## remove high-confidence markers from low-confidence list
    if test -s temp.out
    then
      sort -u temp.out >temp.out.sort
      mv temp.out.sort temp.out

## get rid of all single-strand coverage
      if test -s ${INPUT_FILE_ORIGINAL}.single_strand.lst
      then
        ${FIND_SUB} -i temp.out -c ${INPUT_FILE_ORIGINAL}.single_strand.lst -t 0 -d '\t' -n 0 -o temp.out.clean
        mv temp.out.clean temp.out
      fi

      if test ${confidence} = low_confidence
      then
        if test -s ${HIGH_CCONFIDENCE_OUTPUT_FILE}
        then
          cat ${HIGH_CCONFIDENCE_OUTPUT_FILE} ${HIGH_CCONFIDENCE_OUTPUT_FILE}.repeat |cut -f1 |sort -u>hc_snp.lst
          #cat ${HIGH_CCONFIDENCE_OUTPUT_FILE}  |cut -f1 |sort -u>hc_snp.lst
          ${FIND_SUB} -i temp.out -c hc_snp.lst -t 0 -d '\t' -n 0 -o temp.out.clean
          mv temp.out.clean temp.out
          rm hc_snp.lst
        fi
      fi
    fi

    if test -s temp.out
    then
      cut -f1 temp.out>x
      if test -s x
      then
## clean up the non-specific site first
        cat x |sed /\\./s//:/ >snp_list
##        ${CK_NON_SPECIFIC_PROGRAM} `pwd`/snp_list 50 1 `pwd` `pwd`/non_specific_site.txt
##        ${CK_NON_SPECIFIC_PROGRAM} `pwd`/snp_list 75 1 `pwd` `pwd`/non_specific_site.txt
        if test -s non_specific_site.txt
        then
          rm non_specific_site.txt
        fi

        if test -s non_specific_site.txt ##um? didnt you just do that?
        then
##          awk '{if($2 >0) print $0}' non_specific_site.txt  |cut -f1 |sed /\:/s//./ >non_specific_snp.lst
          awk '{if($2 >=60) print $0}' non_specific_site.txt  |cut -f1 |sed /\:/s//./ >non_specific_snp_high.lst
          if test -s non_specific_snp.lst
          then
           cat non_specific_snp.lst >>${REPEAT_LIST}
          fi
          if test -s non_specific_snp_high.lst
          then
            cat non_specific_snp_high.lst >>${REPEAT_LIST}.high
          fi

        fi
  
        if test -s x
        then
          ${FIND_SUB} -i  ${INPUT_FILE} -c x -t 1 -d '\t' -n 0 -o x.out
## run ptrfinder to clean up the polymers and simple repeats
          if test -s x.out
          then
            cut -f1,8,9,11 x.out |awk '{printf(">%s\n%s%s%s\n", $1, substr($4, 1, 20), $2, substr($4, length($4)-19, 20))}' >input4ptrfinder.seq
## screen for di, tri, etc repeats
            if test -s polymer.out
            then
              rm polymer.out
            fi
            ${PTR_FINDER} -seq input4ptrfinder.seq -repsize 2,5 -minrep 6  |awk '{printf("%s\t%ld\t%ld\n", $1, $10, $12)}' >polymer.out
            ${PTR_FINDER} -seq input4ptrfinder.seq -repsize 1 -minrep 8  |awk '{printf("%s\t%ld\t%ld\n", $1, $10, $12)}' >>polymer.out

## screen with the alternative allele
            cut -f1,8,9,11 x.out |awk '{printf(">%s\n%s%s%s\n", $1, substr($4, 1, 20), $3, substr($4, length($4)-19, 20))}' >input4ptrfinder.seq
            ${PTR_FINDER} -seq input4ptrfinder.seq -repsize 2,5 -minrep 6  |awk '{printf("%s\t%ld\t%ld\n", $1, $10, $12)}'>>polymer.out
#above would be like a CAGCAGCAGCAG repeat
            ${PTR_FINDER} -seq input4ptrfinder.seq -repsize 1 -minrep 8  |awk '{printf("%s\t%ld\t%ld\n", $1, $10, $12)}' >>polymer.out
#above would be like a AAAAAAAA repeat
            if test -s polymer.out
            then
	#Testing for repeat count or size maybe?  Only put the read in the output if 
              awk '{if(21 >=$2-1 && 21 <=$3+1) printf("%s\n", $1)}' polymer.out |cut -f2 -d">" |sort -u >polymer_snp.lst
              if test -s polymer_snp.lst
              then
                if test -s x.out.clean
                then
                  rm x.out.clean
                fi
                if test -s fisher_snp_very_high.lst
                then
                  cat fisher_snp_very_high.lst polymer_snp.lst |sort |uniq -d >d
                  cat d polymer_snp.lst |sort |uniq -u >u
                  mv u  polymer_snp.lst
                fi
                ${FIND_SUB} -i  x.out -c polymer_snp.lst -t 0 -d '\t' -n 0 -o x.out.clean
                mv x.out.clean x.out
                cut -f1 x.out |sort -u >x
              fi
            fi
          fi
        fi

        if test -s x
        then
          rs_count=`grep rs x.out |wc |awk '{printf("%s", $1)}'`
          total_count=`wc x |awk '{printf("%s", $1)}'`
          freq=`echo "$rs_count $total_count" |awk '{printf("%ld", $1*100/$2)}'`
          echo "SNP_list=${SIGNIFICANT_SNP_LIST} Threshold=$THRESHOLD: rs_count=$rs_count,  total_count=$total_count, freq=$freq" 
          echo " " >>${LOG_FILE}
          echo "SNP_list=${SIGNIFICANT_SNP_LIST} Threshold=$THRESHOLD: rs_count=$rs_count,  total_count=$total_count, freq=$freq" >>${LOG_FILE}

          if test -s ${REPEAT_LIST}
          then
            no_repeat_rs_count=`fgrep -f ${REPEAT_LIST} -v x.out |grep rs |wc |awk '{printf("%s", $1)}'`
            no_repeat_total_count=`fgrep -f ${REPEAT_LIST} -v x |wc |awk '{printf("%s", $1)}'`
            if test ${no_repeat_total_count} -gt 0
            then
              freq=`echo "$no_repeat_rs_count $no_repeat_total_count" |awk '{printf("%ld", $1*100/$2)}'`
              echo "No repeat: SNP_list=${SIGNIFICANT_SNP_LIST} Threshold=$THRESHOLD: no_repeat_rs_count=$no_repeat_rs_count, no_repeat_total_count=$no_repeat_total_count, freq=$freq" 
              echo "No repeat: SNP_list=${SIGNIFICANT_SNP_LIST} Threshold=$THRESHOLD: no_repeat_rs_count=$no_repeat_rs_count, no_repeat_total_count=$no_repeat_total_count, freq=$freq" >>${LOG_FILE}
            else
              echo "No repeat: SNP_list=${SIGNIFICANT_SNP_LIST} Threshold=${THRESHOLD}: 0 found"
              echo "No repeat: SNP_list=${SIGNIFICANT_SNP_LIST} Threshold=${THRESHOLD}: 0 found" >>${LOG_FILE}
            fi
          fi
        else
          echo "SNP_list=${SIGNIFICANT_SNP_LIST} Threshold=${THRESHOLD}: 0 found"
          echo " " >>${LOG_FILE}
          echo "SNP_list=${SIGNIFICANT_SNP_LIST} Threshold=${THRESHOLD}: 0 found" >>${LOG_FILE}
        fi

##   if test ${THRESHOLD} = HIGH  ## only run the revision analysis for high percentage
   
       freq=20
       if test $freq -gt 70  ## virtually disable this function
       then
  ## use 90% for determining a reasonable cut off for high-frequency mutation. This may run into trouble if bona-fide mutations ended up having low-coverage
        echo "awk 'BEGIN{percent=0;}{percent+=\$1*100/$rs_count; if(percent >=${FILTER_RS_THRESHOLD}) printf(\"%ld\t%ld\\n\", \$2, percent)}'" >awk.x
        chmod ugo+x awk.x
        grep rs x.out |awk '{printf("%ld\n", $12+$14)}' |sort -n |uniq -c |./awk.x >rs_cvg.out
  ##    grep rs x.out |cut -f13 |sort -n |uniq -c |awk 'BEGIN{percent=0;}{percent+=$1*100/2511; if(percent >=90) printf("%ld\t%ld\n", $2, percent)}' >rs_cvg.out
        if test -s rs_cvg.out
        then
          rm x.out
          NormalCvg=`head -n1 rs_cvg.out|cut -f1`
          echo "Normal coverage = $NormalCvg"
          echo "awk '{if(\$5 >=$NormalCvg) printf(\"%s\\n\", \$1)}' temp.out" >awk.x
          chmod uog+x awk.x
          ./awk.x >x
          if test -s x
          then
            cat x |sed /\\./s//:/ >snp_list
##            ${CK_NON_SPECIFIC_PROGRAM} `pwd`/snp_list 50 1 `pwd` `pwd`/non_specific_site.txt
##            ${CK_NON_SPECIFIC_PROGRAM} `pwd`/snp_list 75 1 `pwd` `pwd`/non_specific_site.txt
            if test -s non_specific_site.txt
            then
              rm non_specific_site.txt
            fi
            if test -s non_specific_site.txt
            then
              awk '{if($2 >0) print $0}' non_specific_site.txt  |cut -f1 |sed /\:/s//./ >non_specific_snp.lst
              ${FIND_SUB} -i x -c non_specific_snp.lst -t 0 -d '\t' -n 0 -o x.clean
              mv x.clean x
            fi
            if test -s x
            then
              ${FIND_SUB} -i  ${INPUT_FILE} -c x -t 1 -d '\t' -n 0 -o x.out
            fi
          fi
        fi
        rs_count=`grep rs x.out |wc |awk '{printf("%s", $1)}'`
        total_count=`wc x |awk '{printf("%s", $1)}'`
        echo "Revised: $THRESHOLD:frequency above 50% rs_count=$rs_count  total_count=$total_count"
      fi
  
      if test ${confidence} = high_confidence  ##keep high-confidence with rs number
      then
        cat x.out |fgrep -f TP53_snp.lst -v |grep -v rs>>${OUTPUT_FILE}
        cat x.out |fgrep -f TP53_snp.lst -v |grep rs >x_rs.out
        if test -s x_rs.out
        then
          fgrep -f fisher_snp_very_high.lst x_rs.out >>${OUTPUT_FILE}
        fi
        rm x_rs.out
      else
        cat x.out |fgrep -f TP53_snp.lst -v |grep -v rs >>${OUTPUT_FILE}
      fi
      cat x.out |fgrep -f TP53_snp.lst -v |grep rs |sed /\^/s//"$THRESHOLD	"/>> ${OUTPUT_RS_FALSE_FILE}
    fi
  fi
  done
done

if test -s ${OUTPUT_FILE}
then
  sed /chrMT/s//chrM/g  ${OUTPUT_FILE} |sort -u > ${OUTPUT_FILE}.sort
  mv  ${OUTPUT_FILE}.sort  ${OUTPUT_FILE}
fi

## Do a final despearate move on filtering repeats
if test -s ${OUTPUT_FILE}
then
  cut -f1 ${OUTPUT_FILE} |sed /\\./s//:/ >snp_list
##  ${CK_NON_SPECIFIC_PROGRAM} `pwd`/snp_list 75 1 `pwd` `pwd`/non_specific_site.txt
  if test -s non_specific_site.txt
  then
    awk '{if($2 >0) print $0}' non_specific_site.txt  |cut -f1 |sed /\:/s//./ >non_specific_snp.lst
    if test -s non_specific_snp.lst
    then
      cat non_specific_snp.lst >>${REPEAT_LIST}
    fi
    awk '{if($2 >=60) print $0}' non_specific_site.txt  |cut -f1 |sed /\:/s//./ >non_specific_snp.lst
    if test -s non_specific_snp.lst
    then
      cat non_specific_snp.lst >>${REPEAT_LIST}.high
    fi
  fi
fi

## non-specific but with <1/10000 fisher pvalue is not treated as non-specific
if test -s d
then
  rm d
fi

if test -s ${REPEAT_LIST}.high
then
## consider those in repeat but has very high-fisher value as something good
  echo "find repeat ${REPEAT_LIST}"
  sort -u ${REPEAT_LIST}.high >${REPEAT_LIST}.sort
  mv ${REPEAT_LIST}.sort ${REPEAT_LIST}.high
  cat ${REPEAT_LIST}.high fisher_snp_very_high.lst |sort |uniq -d >>d
fi

if test -s ${REPEAT_LIST}
then
  if test -s ${REPEAT_LIST}.high
  then
    fgrep -f ${REPEAT_LIST}.high -v ${REPEAT_LIST} |sort -u >>d
  else
    sort -u ${REPEAT_LIST} >>d
  fi
fi

if test -s d
then
    ${FIND_SUB} -i ${OUTPUT_FILE} -c d -t 1 -d '\t' -n 0 -o ${OUTPUT_FILE}.repeat
    if test -s ${OUTPUT_FILE}.repeat
    then
      grep -v rs ${OUTPUT_FILE}.repeat >${OUTPUT_FILE}.repeat.clean
      mv ${OUTPUT_FILE}.repeat.clean ${OUTPUT_FILE}.repeat
    fi
fi

if test -s ${REPEAT_LIST}
then
    ${FIND_SUB} -i ${OUTPUT_FILE} -c ${REPEAT_LIST} -t 0 -d '\t' -n 0 -o ${OUTPUT_FILE}.clean
    mv ${OUTPUT_FILE}.clean ${OUTPUT_FILE}
fi

## clean up
## rm ${INPUT_FILE}.tmp.original
## rm ${INPUT_FILE}.tmp
## rm temp.out
## rm rs_cvg.out
## rm x
## rm x.clean
## rm non_specific_site.txt
## rm non_specific_snp.lst
## rm awk.x
## rm fisher_snp_low.lst 
## rm fisher_snp_high.lst
## rm fisher_snp_very_high.lst
## rm ${REPEAT_LIST}
## rm ${INPUT_FILE}
## rm ${INPUT_FILE_ORIGINAL}.single_strand.out
## rm ${INPUT_FILE_ORIGINAL}.single_strand.lst
done
