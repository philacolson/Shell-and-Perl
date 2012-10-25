#!/bin/sh
## this program creates a fasta *.fa and a *.qual file in a genomic region of interest from *.bam file

BAM_FILE=$1  ## bam file with complete path. e.g /nfs_exports/genomes/1/PCGP/BucketRaw/SJTALL/SJTALL011_D-TB-09-2789.bam
CHR=$2 ## chromosome number, with no chr 
CHR_FROM=$3  ## start_chromosome of interest
CHR_TO=$4    ## end position of a chromosome
ERR_LOG=$5
PRE_FIX4OUTPUT_FILE=$6  ## create the *.fa file or *.fa.qual file

NEXTGEN_BIN_ROOT_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess
SAMTOOL=samtools
FASTQ_TRIM=${NEXTGEN_BIN_ROOT_DIR}/perlsrc/fastq_trim.pl

if test -s ${ERR_LOG}
then
  rm ${ERR_LOG}
fi

if test ! -s ${BAM_FILE}
then
  echo "Fail to find bam file ${BAM_FILE}"
  echo "Fail to find bam file ${BAM_FILE}" >${ERR_LOG}
  exit 1
fi

## create a temporary bam file
## ${SAMTOOL} view -b ${BAM_FILE} 11:32374390-32374590> extract.bam
echo "${SAMTOOL} view -F 1024 -b ${BAM_FILE} ${CHR}:${CHR_FROM}-${CHR_TO} >extract.bam"
${SAMTOOL} view -F 1024 -b ${BAM_FILE} ${CHR}:${CHR_FROM}-${CHR_TO} >extract.bam
if test ! -s extract.bam
then
  echo "Fail to generate bam file. ${SAMTOOL} view -b ${BAM_FILE} $CHR}:${CHR_FROM}-${CHR_TO} extract.bam"
  echo "Fail to generate bam file. ${SAMTOOL} view -b ${BAM_FILE} $CHR}:${CHR_FROM}-${CHR_TO} extract.bam" >${ERR_LOG}
  exit 1
fi

export PERL5LIB=${NEXTGEN_BIN_ROOT_DIR}/MikePerl/SitePerl:${NEXTGEN_BIN_ROOT_DIR}/MikePerl/MultiThread1:${NEXTGEN_BIN_ROOT_DIR}/MikePerl/MultiThread

if test -s ${PRE_FIX4OUTPUT_FILE}
then
  rm ${PRE_FIX4OUTPUT_FILE}
fi
java -jar ${NEXTGEN_BIN_ROOT_DIR}/picard-tools-1.13/SamToFastq.jar INPUT=extract.bam FASTQ=${PRE_FIX4OUTPUT_FILE} VALIDATION_STRINGENCY=SILENT
if test ! -s ${PRE_FIX4OUTPUT_FILE}
then
  echo "Fail in SamtoFastq: java -jar ${NEXTGEN_BIN_ROOT_DIR}/picard-tools-1.13/SamToFastq.jar INPUT=extract.bam FASTQ=${PRE_FIX4OUTPUT_FILE} VALIDATION_STRINGENCY=SILENT" 
  echo "Fail in SamtoFastq: java -jar ${NEXTGEN_BIN_ROOT_DIR}/picard-tools-1.13/SamToFastq.jar INPUT=extract.bam FASTQ=${PRE_FIX4OUTPUT_FILE} VALIDATION_STRINGENCY=SILENT" >${ERR_LOG}
  exit 1
fi
awk 'BEGIN{i=0;}{++i; val=i%4; if(val==1) printf("%s|%s%ld\n", $1, $1, i)}' ${PRE_FIX4OUTPUT_FILE} >>old_new_id_map.txt
awk 'BEGIN{i=0;}{++i; val=i%4; if(val==1) printf("%s%ld\n", $1, i); else print $0}' ${PRE_FIX4OUTPUT_FILE} >${PRE_FIX4OUTPUT_FILE}.mod
mv ${PRE_FIX4OUTPUT_FILE}.mod ${PRE_FIX4OUTPUT_FILE}
${FASTQ_TRIM} -convert-fastq-to-fasta ${PRE_FIX4OUTPUT_FILE}
if test ! -s ${PRE_FIX4OUTPUT_FILE}.fa
then
  echo "Fail to build fasta file ${FASTQ_TRIM} -convert-fastq-to-fasta ${PRE_FIX4OUTPUT_FILE}" 
  echo "Fail to build fasta file ${FASTQ_TRIM} -convert-fastq-to-fasta ${PRE_FIX4OUTPUT_FILE}" >${ERR_LOG}
  exit 1
fi

if test ! -s ${PRE_FIX4OUTPUT_FILE}.fa.qual
then
  echo "Fail to build quality file ${FASTQ_TRIM} -convert-fastq-to-fasta ${PRE_FIX4OUTPUT_FILE}" 
  echo "Fail to build quality file ${FASTQ_TRIM} -convert-fastq-to-fasta ${PRE_FIX4OUTPUT_FILE}" >${ERR_LOG}
  exit 1
fi

rm ${PRE_FIX4OUTPUT_FILE}
rm extract.bam

