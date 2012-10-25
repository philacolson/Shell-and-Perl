#!/bin/sh

## go through the tempReview_SampleName.txt dir and check whether sites listed in *_putative_exon_somatic_indel_mutation.txt
## match dbSNP insertion/deletion event
## nees to be incoporated

SAMPLE_NAME=$1  ## e.g. an file with output of DONOR_tumor.txt
DIR=$2
OUTPUT_FILE=$3  ##e.g. SJRB002_indel_overlap_dbSNP.txt
MODIFY_MUTATION_FILE=$4
## dbSNP_DELETION_DIR=/nfs_exports/apps/gnu-apps/NextGen/dbSNP_DELETION_REGION
## dbSNP_INSERTION_DIR=/nfs_exports/apps/gnu-apps/NextGen/dbSNP_INSERTION_REGION
dbSNP_DELETION_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess/dbSNP_DELETION_REGION
dbSNP_INSERTION_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess/dbSNP_INSERTION_REGION

MUTATION_FILE=${SAMPLE_NAME}_putative_exon_somatic_indel_mutation.txt

if test ! -d ${DIR}
then
  echo "no indel analysis dir ${DIR}"
  exit 1
fi

cd ${DIR}

if test ! -s ${MUTATION_FILE}
then
  echo "Fail to find file ${MUTATION_FILE}"
  exit 1
fi

if test ! ${SAMPLE_NAME}_tumor.txt
then
  echo "Fail to find file ${SAMPLE_NAME}_tumor.txt"
  exit 1
fi

cut -f4,5 ${MUTATION_FILE} |sed /^chr/s/// |grep -v Chr >somatic_exon_indel.lst
if test ! -s somatic_exon_indel.lst
then
  echo "Fail to find somatic indel file somatic_exon_indel.lst for ${SAMPLE_NAME}"
  exit 1
fi


fgrep -f  somatic_exon_indel.lst ${SAMPLE_NAME}_tumor.txt |grep insertion |sort -u>${SAMPLE_NAME}_insertion.txt
fgrep -f  somatic_exon_indel.lst ${SAMPLE_NAME}_tumor.txt |grep deletion |sort -u>${SAMPLE_NAME}_deletion.txt

cut -f1 somatic_exon_indel.lst |sort -u>all_indel_chr.lst
INPUT_CHR_LIST=all_indel_chr.lst



NEXTGEN_BIN_ROOT_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess
PERLSRC_DIR=${NEXTGEN_BIN_ROOT_DIR}/perlsrc
FIND_OVERLAP=${PERLSRC_DIR}/find_overlap_interval.pl
FIND_SUBP=${PERLSRC_DIR}/FindSub.pl


if test -s ${OUTPUT_FILE}
then
  rm ${OUTPUT_FILE}
fi

for i in `cat ${INPUT_CHR_LIST}`; do
   echo $i
   echo ${i} >t_chr.lst
   ${FIND_SUBP} -i ${SAMPLE_NAME}_insertion.txt -c t_chr.lst -d '\t' -n 1 -t 1 -o t_chr.lst.out
   if test -s t_chr.lst.out
   then
      if test -s t_chr_dbsnp_overlap.out
      then
        rm t_chr_dbsnp_overlap.out
      fi
      ${FIND_OVERLAP} -i t_chr.lst.out -d ${dbSNP_INSERTION_DIR}/chr${i}_region.txt -o t_chr_dbsnp_overlap.out -e err.log -m 3 -n 4 -q 1
      if test -s t_chr_dbsnp_overlap.out
      then
        cat t_chr_dbsnp_overlap.out >>${OUTPUT_FILE}
      fi
    fi
       
   ${FIND_SUBP} -i ${SAMPLE_NAME}_deletion.txt -c t_chr.lst -d '\t' -n 1 -t1 -o t_chr.lst.out
   if test -s t_chr.lst.out
   then
      if test -s t_chr_dbsnp_overlap.out
      then
        rm t_chr_dbsnp_overlap.out
      fi
      ${FIND_OVERLAP} -i t_chr.lst.out -d ${dbSNP_DELETION_DIR}/chr${i}_region.txt -o t_chr_dbsnp_overlap.out -e err.log -m 3 -n 4 -q 1
      if test -s t_chr_dbsnp_overlap.out
      then
        cat t_chr_dbsnp_overlap.out >>${OUTPUT_FILE}
      fi
    fi
done

if test ${MODIFY_MUTATION_FILE} = 1
then
  cut -f2,3 ${OUTPUT_FILE} |sort -u >dbSNP_overlap.lst
  if test -s dbSNP_overlap.lst
  then
    head -n1 ${MUTATION_FILE} >${MUTATION_FILE}.original
    fgrep -f dbSNP_overlap.lst -v ${MUTATION_FILE} |grep -v GeneName |sed /\$/s//\\\tGood/ >>${MUTATION_FILE}.original
    fgrep -f dbSNP_overlap.lst ${MUTATION_FILE} |sed /\$/s//\\\tdbSNP/ >>${MUTATION_FILE}.original
    grep -v dbSNP ${MUTATION_FILE}.original >${MUTATION_FILE}
  else
    sed /\$/s//\\\tGood/ ${MUTATION_FILE} >${MUTATION_FILE}.original
    cp ${MUTATION_FILE}.original ${MUTATION_FILE}
  fi
fi
