#!/bin/sh

## this program is designed to merge the low_confidence_call and high_confidence call from two different groups
## low-confidence call will capture the variations that were missed for the following reasons
## a) quality check failed due to lack of the 5bp-flanking region or the 10bp-flanking threshold
## b) mismatch check failed due to the low-quality region mismatch
## False-positives
## a) Too-many high-quality mismatches
## b) indels
## c) dbsnp (no check on dbSNP)
## filtering criteria
## 1) remove everything that overlaps with dbSNP
## 2) remove everythig that overlaps with indel data
## 3) require a minimum count of minor alleles (>2 or a high frequency >=20%)

INPUT_DIR4_LOW_QUALITY_CALL=$1  ## e.g. temp3 for SJTALL007_high_confidence_somatic_snp.txt
INPUT_DIR4_HIGH_QUALITY_CALL=$2
SAMPLE_NAME=$3  ## e.g. SJTALL007
INPUT_BAM_FILE_4_INDEL_CHECK=$4  ## e.g SJTALL007_bam_high_20.out
OUTPUT_DIR=$5

HIGH_CONFIDENCE_TAG=high_confidence_somatic_snp
LOW_CONFIDENCE_TAG=low_confidence_somatic_snp
FIND_SUB=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess/perlsrc/FindSub.pl

HIGH_CONFIDENCE_FILE=${SAMPLE_NAME}_${HIGH_CONFIDENCE_TAG}.txt
LOW_CONFIDENCIDE_FILE=${SAMPLE_NAME}_${LOW_CONFIDENCE_TAG}.txt

DB_SNP_SNV_LIST=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess/dbsnp130_snv.lst

if test ! -d ${INPUT_DIR4_LOW_QUALITY_CALL}
then
if [ $DEBUG_LEVEL -gt 0 ]
  then
  echo "Fail to find directory ${INPUT_DIR4_LOW_QUALITY_CALL}"
  fi
  exit 1
fi

if test ! -d ${INPUT_DIR4_HIGH_QUALITY_CALL}
then
  if [ $DEBUG_LEVEL -gt 0 ]
  then
  echo "Fail to find directory ${INPUT_DIR4_HIGH_QUALITY_CALL}"
  fi
  exit 1
fi

if test ! -s ${INPUT_DIR4_LOW_QUALITY_CALL}/${HIGH_CONFIDENCE_FILE}
then
if [ $DEBUG_LEVEL -gt 0 ]
  then
  echo "Fail to find high-confidence call for low-quality data ${INPUT_DIR4_LOW_QUALITY_CALL}/${HIGH_CONFIDENCE_FILE}"
  fi
  exit 1
fi

if test ! -s ${INPUT_DIR4_LOW_QUALITY_CALL}/${LOW_CONFIDENCIDE_FILE}
then
  if [ $DEBUG_LEVEL -gt 0 ]
  then
  echo "Fail to find low-confidence call for low-quality data ${INPUT_DIR4_LOW_QUALITY_CALL}/${LOW_CONFIDENCIDE_FILE}"
  fi
  exit 1
fi

if test ! -s ${INPUT_DIR4_HIGH_QUALITY_CALL}/${HIGH_CONFIDENCE_FILE}
then
  if [ $DEBUG_LEVEL -gt 0 ]
  then
  echo "Fail to find high-confidence call for high-quality data ${INPUT_DIR4_HIGH_QUALITY_CALL}/${HIGH_CONFIDENCE_FILE}"
  fi
  exit 1
fi

if test ! -s ${INPUT_DIR4_HIGH_QUALITY_CALL}/${LOW_CONFIDENCIDE_FILE}
then
  if [ $DEBUG_LEVEL -gt 0 ]
  then
  echo "Fail to find low-confidence call for high-quality data ${INPUT_DIR4_HIGH_QUALITY_CALL}/${LOW_CONFIDENCIDE_FILE}"
  fi
  exit 1
fi

if test ! -s ${INPUT_BAM_FILE_4_INDEL_CHECK}
then
  if [ $DEBUG_LEVEL -gt 0 ]
  then
  echo "Fail to find the input indel file ${INPUT_BAM_FILE_4_INDEL_CHECK}"
  fi
  exit 1
fi

INDEL_LIST=${SAMPLE_NAME}_bam_indel.lst
if test -s ${INDEL_LIST}
then
  rm ${INDEL_LIST}
fi
if [ $DEBUG_LEVEL -gt 0 ]
  then
  echo "create indel list in merge_confidence_file"
fi
grep -v SNP ${INPUT_BAM_FILE_4_INDEL_CHECK}  |awk '{for(i=0-$5; i<=$5; ++i) printf("%s.%ld\n", $2, $3+i)}' |sort -u >${INDEL_LIST}



if test ! -d ${OUTPUT_DIR}
then
  mkdir ${OUTPUT_DIR}
fi

## collect all high-quality data (including both high and low confidence)
cat ${INPUT_DIR4_HIGH_QUALITY_CALL}/${HIGH_CONFIDENCE_FILE} ${INPUT_DIR4_HIGH_QUALITY_CALL}/${LOW_CONFIDENCIDE_FILE} |cut -f1|sort -u >hq_snp.lst

echo ${HIGH_CONFIDENCE_FILE} >input_file.lst
echo ${LOW_CONFIDENCIDE_FILE} >>input_file.lst
echo ${HIGH_CONFIDENCE_FILE}.repeat >>input_file.lst

for j in `cat input_file.lst`; do
  echo ${j}
  i=${INPUT_DIR4_LOW_QUALITY_CALL}/${j}
  if [ $DEBUG_LEVEL -gt 0 ]
  then
  echo $i
  echo "filter existing records"
  echo "${FIND_SUB} -i ${i} -c hq_snp.lst -t 0 -d '\t' -n 0 -o ${i}.new"
  fi
  ${FIND_SUB} -i ${i} -c hq_snp.lst -t 0 -d '\t' -n 0 -o ${i}.new
  if test -s ${i}.new
  then
## clean-up dbSNP
    if [ $DEBUG_LEVEL -gt 0 ]
  then
  echo "filter dbSNP"
  fi
    ${FIND_SUB} -i ${i}.new -c ${DB_SNP_SNV_LIST} -t 0 -d '\t' -n 0 -o ${i}.new.2
    rm ${i}.new
    mv ${i}.new.2 ${i}.new
## clean-up indels
    if test -s ${i}.new
    then
      if [ $DEBUG_LEVEL -gt 0 ]
  then
  echo "filter indel"
  fi
      ${FIND_SUB} -i ${i}.new -c ${INDEL_LIST} -t 0 -d '\t' -n 0 -o ${i}.new.2
      rm ${i}.new
      mv ${i}.new.2 ${i}.new
    fi
## evaluate the coverage for the minor allele. requires >=3 minor allele or >=2 alleles but with MAF >=20% and normal >=10x coverage
    if test -s ${i}.new
    then
      grep -v Name ${i}.new |grep SNP |awk '{if($15 >=3 || ($15 >=2 && $15*100/($13+$15) >=20 && $12 >=10)) print $0}' >${i}.new.2
      if test -s ${i}.new.2
      then
        new_count=`wc ${i}.new.2 |awk '{printf("%s", $1)}'`
        original_count=`wc ${INPUT_DIR4_HIGH_QUALITY_CALL}/${j} |awk '{printf("%s", $1)}'`
        if [ $DEBUG_LEVEL -gt 0 ]
  then
        echo "original_count=${original_count} new_count=${new_count}"
        echo "cat ${INPUT_DIR4_HIGH_QUALITY_CALL}/${j} ${i}.new.2"
        echo "OUTPUT file=${OUTPUT_DIR}/${j}"
        fi
        cat ${INPUT_DIR4_HIGH_QUALITY_CALL}/${j} ${i}.new.2 >${OUTPUT_DIR}/${j}
      fi
##      rm ${i}.new.2
    fi
##    rm ${i}.new
  fi
done

