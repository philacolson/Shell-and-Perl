#!/bin/sh
## this program was a derivative of manual review
## it was derived from run_missing_indel
## it is expected to be used for two purposes:
## a) perform analysis to rescue the indels calculated from Mike's *bam.out low and high. Sometimes the assembly results can miss indels due to the quality value used in extraction
## b) perform analysis to verify the somatic indels collected from either the assembly process or from the IndelRescue process. There were a lot of falsely identified indels due to issues related to ExtractUnmapped
## the program uses all reads in the bam file to run the evaluation
## this version should be used only for assessing the state of tumor sample


CHR_LIST=$1

SAMPLE_MATCHING_SAMPLE_BAM_LIST=$2  ## input looks like sample|sample_bam|match_sample_bam. eg. SJTALL006|/nfs_exports/genomes/1/PCGP/BucketRaw/SJTALL/SJTALL006_D-TB-06-0060.bam|/nfs_exports/genomes/1/PCGP/BucketRaw/SJTALL/SJTALL006_G-TB-06-0386.bam

DIR4RUN_verify_indel=$3   ## put the data into one of the subdirectories that have already indexes built for unmapped reads. 
## DIR4RUN_verify_indel=/tcga_next_gen/TCGA_WholeGenome/IndelAnalysis/Velvet_trim_unmapped_k21.  /nfs_exports/genomes/1/PCGP/BucketIntermediate/SJTALL/IndelAnalysis/Velvet_trim_unmapped_k21

INPUT_POS_LIST=$4   ## this can be either a 2-field (chr|pos or three field chr|pos|indel_size. the 3rd field is useful for dealing with large indel size

CLEAN_ALL_INTERMEDIATE_FILES=$5  ## set to 1 if all intermediate files need to be cleaned

OUTPUT_DIR=$6   ## directory for collecting the results in the format  *_tumor.txt *_normal.txt

NEXTGEN_BIN_ROOT_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess
BIN_DIR=${NEXTGEN_BIN_ROOT_DIR}
PERLSRC_DIR=${NEXTGEN_BIN_ROOT_DIR}/perlsrc
## AUTO_REVIEW_PROGRAM=${BIN_DIR}/verify_Nextgen_indel_by_bam_reads  ## this is the program that performs automated review of bam reads
AUTO_REVIEW_PROGRAM=${NEXTGEN_BIN_ROOT_DIR}/verify_Nextgen_indel_by_bam_reads
## AUTO_REVIEW_PROGRAM=/nfs_exports/genomes/1/PCGP/BucketIntermediate/SJNBL/IndelAnalysis/verify_Nextgen_indel_by_bam_reads

FIND_SUB=${PERLSRC_DIR}/FindSub.pl

if test ! -d ${DIR4RUN_verify_indel}
then
  echo "Fail to find directory for running verify indel ${DIR4RUN_verify_indel}"
  exit 1
fi

if test ! -s ${INPUT_POS_LIST}
then
  echo "Fail to find input file ${INPUT_POS_LIST}"
  exit 1
fi

TAG4INPUT=reviewed

for j in `cat ${SAMPLE_MATCHING_SAMPLE_BAM_LIST}`; do
  SAMPLE_NAME=`echo $j|cut -f1 -d"|"`
  SAMPLE_BAM=`echo $j |cut -f2 -d"|"`
  MATCH_SAMPLE_BAM=`echo $j |cut -f3 -d"|"`
  DONOR_NAME=`echo ${SAMPLE_NAME}|cut -f1 -d"_"`
  TISSUE_STATUS=`echo ${SAMPLE_NAME}|cut -f2 -d"_"`
  disease_tag=normal
  if test ${TISSUE_STATUS} = D
  then
    disease_tag=tumor
  else
    if test ${TISSUE_STATUS} = R
    then
      disease_tag=relapse
    else
      if test ${TISSUE_STATUS} = X
      then
        disease_tag=Xenograph
      fi
    fi
  fi
  echo "DONOR=${DONOR_NAME} disease=${disease_tag}"

  
## clean up possible existing files
  rm ${DIR4RUN_verify_indel}/${SAMPLE_NAME}_*_${TAG4INPUT}_chr*_pos.lst_rpt.txt

  INPUT_FILE=${INPUT_POS_LIST}
  if test ! -s ${INPUT_FILE}
  then
    echo "Fail to find input file with the indel location ${INPUT_FILE}"
  else
    echo $j
    for i in `cat ${CHR_LIST}`; do
      echo $i
##      if test -d ${DIR4RUN_verify_indel}/run_${SAMPLE_NAME}.${i}  ## this may create a problem if jz method fails to find something
        if test -d ${DIR4RUN_verify_indel}
      then
##        PROCESS_DIR=${DIR4RUN_verify_indel}/run_${SAMPLE_NAME}.${i}
        PROCESS_DIR=${DIR4RUN_verify_indel}
        cd ${PROCESS_DIR}
        grep ^"chr${i}|" ${INPUT_FILE} >${TAG4INPUT}_chr${i}_pos.lst
        if test -s ${TAG4INPUT}_chr${i}_pos.lst
        then
          echo "${AUTO_REVIEW_PROGRAM} ${DIR4RUN_verify_indel}  ${TAG4INPUT}_chr${i}_pos.lst ${SAMPLE_NAME} $i ${DIR_4_ASSEMBLY} ${BIN_DIR} 1 ${SAMPLE_BAM} ${MATCH_SAMPLE_BAM}"
          ${AUTO_REVIEW_PROGRAM} ${DIR4RUN_verify_indel}  ${TAG4INPUT}_chr${i}_pos.lst ${SAMPLE_NAME} $i ${DIR_4_ASSEMBLY} ${BIN_DIR} ${CLEAN_ALL_INTERMEDIATE_FILES} ${SAMPLE_BAM} ${MATCH_SAMPLE_BAM}
        fi
      fi
    done
    cd ${DIR4RUN_verify_indel}
    ls ${SAMPLE_NAME}_*_${TAG4INPUT}_chr*_pos.lst_rpt.txt >rescue.lst
    if test -s rescue.lst
    then
      cat ${SAMPLE_NAME}_*_${TAG4INPUT}_chr*_pos.lst_rpt.txt |awk '{if($NF>=0) print $0}' |awk '{printf("%ld\t", length($12)); print $0}' |sort +0 -1 -n -r |cut -f1-6,8,13,14,15,16,17 |sed /\$/s//"\\tIndelRescue"/ >${OUTPUT_DIR}/${DONOR_NAME}_${disease_tag}.txt
    fi
  fi
  rm ${DIR4RUN_verify_indel}/${SAMPLE_NAME}_*_${TAG4INPUT}_chr*_pos.lst_rpt.txt
done
