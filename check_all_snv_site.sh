#!/bin/sh
## input file format
## 9|139038415|A|1
## 17|7190283|G|1
## col1=chr; col2=pos; col3=mutant allele; col4=check_single_strand
## some of the sites only have single-strand coverage. Those need to be excluded from coverage analysis

PROGRAM=/nfs_exports/apps/internal/scripts/check_snv_variation_call
## BAM_FILE=/nfs_exports/genomes/1/projects/RNASEQ/PCGP/BucketRaw/SJTALL002_D-TB-05-2043-reordered.bam
## INPUT_FILE=SJTALL012_test_input.txt
INPUT_FILE=$1
BAM_FILE=$2
GERMLINE_BAM_FILE=$3   ## this is optional. use no_germline_bam file as a possibility
OUTPUT_FILE=$4
## OUTPUT_FILE=/nfs_exports/genomes/1/projects/RNASEQ/PCGP/BucketIntermediate/SnpDetect/SJTALL002_SNV_var_report.txt

if test -s ${OUTPUT_FILE}
then 
  rm ${OUTPUT_FILE}
fi

if test ! -s ${BAM_FILE}
then
  echo "Fail to find bam file ${BAM_FILE}"
  exit 1
fi

if test ! -s ${INPUT_FILE}
then
  echo "Fail to find input file ${INPUT_FILE}"
  exit 1
fi

for i in `cat ${INPUT_FILE}`; do
##   BAM_FILE=`echo $i |cut -f1 -d"|"`
  CHR=`echo $i |cut -f1 -d"|"`
  POS=`echo $i |cut -f2 -d"|"`
  ALLELE=`echo $i |cut -f3 -d"|"`
  IGNORE_STRAND_CHECK=`echo $i|cut -f4 -d"|"`
  CK_GERMLINE_ALLELE_COUNT=`echo $i|cut -f5 -d"|"`
  USE_GERMLINE_BAM_FILE=${GERMLINE_BAM_FILE}

## CK_GERMLINE_ALLELE_COUNT means that presence of mutant allele in germline by itself does not disqualify. But
## if germline allele count is >10% of the mutant allele in tumor, then it is considered a poor hit
  echo "${PROGRAM} ${BAM_FILE} ${CHR} ${POS} ${ALLELE} ${OUTPUT_FILE} ${IGNORE_STRAND_CHECK} ${USE_GERMLINE_BAM_FILE} ${CK_GERMLINE_ALLELE_COUNT}"
  ${PROGRAM} ${BAM_FILE} ${CHR} ${POS} ${ALLELE} ${OUTPUT_FILE} ${IGNORE_STRAND_CHECK} ${USE_GERMLINE_BAM_FILE} ${CK_GERMLINE_ALLELE_COUNT}
done
