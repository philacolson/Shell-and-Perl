#!/bin/sh

SAMPLE_NAME=$1
PROCESS_DIR=$2
BAD_SNP_FILE=$3  ## use no_bad_snp_file for default

OUTPUT_SUMMARY_FILE=${PROCESS_DIR}/${SAMPLE_NAME}_sub_sum.txt
OUTPUT_VALIDATED_FILE=${PROCESS_DIR}/${SAMPLE_NAME}_validate_somatic_loci.txt

GENE_EXON_REGION_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess/GENE_EXON_REGION

GENOME_INFO_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess
FIND_SUB=${GENOME_INFO_DIR}/perlsrc/FindSub.pl

SCRIPT_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess
FIND_GENE_PROGRAM=${SCRIPT_DIR}/find_gene4sub



HIGH_CONFIDENCE_SOMATIC_SNP_FILE=${SAMPLE_NAME}_high_confidence_somatic_snp.txt
## high-confidence but has repeat sequence
HIGH_CONFIDENCE_SOMATIC_SNP_FILE_REPEAT=${HIGH_CONFIDENCE_SOMATIC_SNP_FILE}.repeat

LOW_CONFIDENCE_SOMATIC_SNP_FILE=${SAMPLE_NAME}_low_confidence_somatic_snp.txt
## low-confidence and has repetitive regions
LOW_CONFIDENCE_SOMATIC_SNP_FILE_REPEAT=${LOW_CONFIDENCE_SOMATIC_SNP_FILE}.repeat

if test ! -x ${FIND_GENE_PROGRAM}
then
 if [ $DEBUG_LEVEL -gt 0 ]
then
echo "Fail to find program ${FIND_GENE_PROGRAM}"
fi
exit 1
fi

if test ! -d ${GENE_EXON_REGION_DIR}
then
if [ $DEBUG_LEVEL -gt 0 ]
then
echo "Fail to find gene-exon direcotry ${GENE_EXON_REGION_DIR}"
fi
exit 1
fi

if test ! -d ${PROCESS_DIR}
then
if [ $DEBUG_LEVEL -gt 0 ]
then
echo "Fail to find diretory ${PROCESS_DIR}"
fi
exit 1
fi

cd ${PROCESS_DIR}
## for coding SNP comparisons, just get all of them for now
if test -d temp
then
  rm -rf temp
fi
mkdir temp

cat ${HIGH_CONFIDENCE_SOMATIC_SNP_FILE} ${HIGH_CONFIDENCE_SOMATIC_SNP_FILE_REPEAT} ${LOW_CONFIDENCE_SOMATIC_SNP_FILE} ${LOW_CONFIDENCE_SOMATIC_SNP_FILE_REPEAT} |sort -u >${SAMPLE_NAME}_somatic_snp.txt

if test -s ${SAMPLE_NAME}_somatic_snp.txt
then
  if test -s ${BAD_SNP_FILE}
  then
    ${FIND_SUB} -i ${SAMPLE_NAME}_somatic_snp.txt -c ${BAD_SNP_FILE} -t 0 -d '\t' -n 0 -o ${SAMPLE_NAME}_somatic_snp.txt.clean
    mv ${SAMPLE_NAME}_somatic_snp.txt.clean ${SAMPLE_NAME}_somatic_snp.txt
  fi
fi
if test ! -s ${SAMPLE_NAME}_somatic_snp.txt
then
 if [ $DEBUG_LEVEL -gt 0 ]
then echo "Fail to create ${SAMPLE_NAME}_somatic_snp.txt"
  echo "cat ${HIGH_CONFIDENCE_SOMATIC_SNP_FILE} ${HIGH_CONFIDENCE_SOMATIC_SNP_FILE_REPEAT} ${LOW_CONFIDENCE_SOMATIC_SNP_FILE} ${LOW_CONFIDENCE_SOMATIC_SNP_FILE_REPEAT}"
 fi
exit 1
fi

mv ${SAMPLE_NAME}_somatic_snp.txt temp/.
cd temp
mkdir SOMATIC_RESULTS

${FIND_GENE_PROGRAM} `pwd` ${BAD_SNP_FILE} `pwd`/SOMATIC_RESULTS/sub_sum.txt `pwd`/SOMATIC_RESULTS/validate_somatic_loci.txt ${GENE_EXON_REGION_DIR} somatic_snp.txt

if test ! -s SOMATIC_RESULTS/sub_sum.txt
then
  if [ $DEBUG_LEVEL -gt 0 ]
then
echo "Fail to build sub_sum.txt after running ${FIND_GENE_PROGRAM} `pwd` ${BAD_SNP_FILE} `pwd`/SOMATIC_RESULTS/sub_sum.txt `pwd`/SOMATIC_RESULTS/validate_somatic_loci.txt ${GENE_EXON_REGION_DIR} somatic_snp.txt"
fi
exit 1
fi

OUTPUT_SUMMARY_FILE=${PROCESS_DIR}/${SAMPLE_NAME}_sub_sum.txt
OUTPUT_VALIDATED_FILE=${PROCESS_DIR}/${SAMPLE_NAME}_validate_somatic_loci.txt
cd ${PROCESS_DIR}
if test -s ${OUTPUT_SUMMARY_FILE}
then
  rm ${OUTPUT_SUMMARY_FILE}
fi
head -n1 temp/SOMATIC_RESULTS/sub_sum.txt >${OUTPUT_SUMMARY_FILE}

if test -s ${OUTPUT_VALIDATED_FILE}
then
  rm ${OUTPUT_VALIDATED_FILE}
fi


awk '{printf("%s.%ld\t", $4, $5); print $0}' temp/SOMATIC_RESULTS/sub_sum.txt  >temp/SOMATIC_RESULTS/sub_sum.txt.mod
awk '{printf("chr%s.%ld\t", $2, $3); print $0}' temp/SOMATIC_RESULTS/validate_somatic_loci.txt >temp/SOMATIC_RESULTS/validate_somatic_loci.txt.mod


if test -s ck_file.lst
then
  rm ck_file.lst
fi
echo ${HIGH_CONFIDENCE_SOMATIC_SNP_FILE} >>ck_file.lst
echo ${HIGH_CONFIDENCE_SOMATIC_SNP_FILE_REPEAT} >>ck_file.lst
echo ${LOW_CONFIDENCE_SOMATIC_SNP_FILE} >>ck_file.lst
echo ${LOW_CONFIDENCE_SOMATIC_SNP_FILE_REPEAT} >>ck_file.lst

for i in `cat ck_file.lst`; do
  if test -s ${i}
  then
    if test $i = ${HIGH_CONFIDENCE_SOMATIC_SNP_FILE}
    then
      label_val=SJHQ
    else
      if test $i = ${HIGH_CONFIDENCE_SOMATIC_SNP_FILE_REPEAT}
      then
        label_val=SJHQR
      else
        if test $i = ${LOW_CONFIDENCE_SOMATIC_SNP_FILE}
        then
          label_val=SJLQ
        else
          label_val=SJLQR
        fi
      fi 
    fi

    cut -f1 ${i} |sort -u >c_snp.lst
    ${FIND_SUB} -i temp/SOMATIC_RESULTS/sub_sum.txt.mod -c c_snp.lst -t 1 -d '\t' -n 0 -o c_snp.lst.out
    if test -s c_snp.lst.out
    then
##      echo ">${label_val}" >>${OUTPUT_SUMMARY_FILE}
      cut -f2-100 c_snp.lst.out |sed /^/s//"${label_val}	"/>>${OUTPUT_SUMMARY_FILE}
    fi

    ${FIND_SUB} -i temp/SOMATIC_RESULTS/validate_somatic_loci.txt.mod -c c_snp.lst -t 1a -d '\t' -n 0 -o c_snp.lst.out
    if test -s c_snp.lst.out
    then
##      echo ">${label_val}" >>${OUTPUT_VALIDATED_FILE}
      cut -f2-100 c_snp.lst.out |sed /^/s//"${label_val}	"/>>${OUTPUT_VALIDATED_FILE}
    fi
  fi
done
