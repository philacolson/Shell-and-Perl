#!/bin/sh

BAM_FILE=$1  ## bam file for data analysis
CHR=$2
START=$3
END=$4
OUTPUT_FILE=$5

## output file will be named as ${SAMPLE_NAME}_chr[1-x]_coverage.txt
NEXTGEN_BIN_ROOT_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess

export CLASSPATH=${NEXTGEN_BIN_ROOT_DIR}/av.jar:${NEXTGEN_BIN_ROOT_DIR}/sam-1.13.jar:${NEXTGEN_BIN_ROOT_DIR}/mysql-connector-java-5.1.10-bin.jar

if test ! -s ${BAM_FILE}
then
  echo "Fail to open bam file ${BAM_FILE}"
  exit 1
fi

if test ! -s ${BAM_FILE}.bai
then
  echo "Fail to find the index file  ${BAM_FILE}.bai"
  exit 1
fi

java  -Xmx3072m Ace2.SAMCoverage -bam ${BAM_FILE} -tname chr${CHR} -tstart ${START} -tend ${END} -of ${OUTPUT_FILE}
