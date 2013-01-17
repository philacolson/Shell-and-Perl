#!/bin/sh
INPUT_SAMPLE=$1 #list of chromosomes to test
D_BAM_FILE=$2
G_BAM_FILE=$3
##BAM_DIR=/nfs_exports/genomes/1/PCGP/BucketRaw/${DISEASE_CODE}
#BAM_DIR=/user/songliu/u2/group/Qiang/Exome/bladder_output/${DISEASE_CODE}
BAM_DIR=$4
DEBUG_LEVEL=$5
#Debug level can be 1-3, increasing with verbosity.  As a general rule, the grades pertain to the following: 1 just logs fatal events and ones that otherwise prematurely end the code, 2 mentions anything that deviates from expected code/error handler code is run, or other noteworthly events.  3 logs variable names before and after programs are run, and upon initialization.  Higher numbers post their description + any number's description lower than their own.  For example, a debug level of 2 would post every log statement that is not a level 3 log event.  It would post 1's, and 2's.

##SCRIPT_DIR=/nfs_exports/apps/internal/scripts
SCRIPT_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess


FILTER_DBSNP=$SCRIPT_DIR/filter_low_quality_dbSNP_4_tier1 


##EXCEL_ROOT_DIR=/nfs_exports/linux-file1/home/naevegrp/jzhang2/NextGen/PCGP
EXCEL_ROOT_DIR=$6
HIGH_LIST=$7
LOW_LIST=$8
GERM_LIST=$9
LOG=${BAM_DIR}/${INPUT_SAMPLE}_LOG.txt
echo $DEBUG_LEVEL >> $LOG
`cat $LOG`
#reimplement below if we want a list of them to be entered
#for i in `cat ${INPUT_SAMPLE}`; do
  if test ! -d ${EXCEL_ROOT_DIR}/${INPUT_SAMPLE} #if there is no directory for the SNP
  then
	if [ $DEBUG_LEVEL -gt 1 ]
	then 
  echo "INFO: no dir for ${EXCEL_ROOT_DIR}/${INPUT_SAMPLE}, making directory." >> $LOG
	fi
    mkdir -p ${EXCEL_ROOT_DIR}/${INPUT_SAMPLE}
    #exit 1
  fi
done

if [ $DEBUG_LEVEL -gt 2 ]
then 
echo "DEBUG: created EXCEL directories." >> $LOG
fi

for i in `cat ${INPUT_SAMPLE}`; do
  ##echo $i
  cp ../${INPUT_SAMPLE}_G_bam.out .
  cp ../${INPUT_SAMPLE}_bam_low.out ${INPUT_SAMPLE}_bam_low.out
  cp ../${INPUT_SAMPLE}_bam_high_20.out ${INPUT_SAMPLE}_bam_high_20.out
if [ $DEBUG_LEVEL -gt 2 ]
	then echo "DEBUG:  ${SCRIPT_DIR}/run_somatic_mutation_analysis ${INPUT_SAMPLE} no_false_snp" >> $LOG
	fi
  #hq
  #D_BAM=`ls ${BAM_DIR}/${INPUT_SAMPLE}_D*.bam`
  #G_BAM=`ls ${BAM_DIR}/${INPUT_SAMPLE}_G*.bam`
  ${SCRIPT_DIR}/run_somatic_mutation_analysis ${i} no_false_snp ${D_BAM_FILE} ${G_BAM_FILE} ${DEBUG_LEVEL} ${LOG}
done
if [ $DEBUG_LEVEL -gt 2 ]
	then
echo "DEBUG:  ${SCRIPT_DIR}/run_check_snv_script ${INPUT_SAMPLE} ${DISEASE_CODE}" >> $LOG
fi

${SCRIPT_DIR}/run_check_snv_script ${INPUT_SAMPLE} ${DISEASE_CODE}

for i in `cat ${INPUT_SAMPLE}`; do
  echo $i
  SAMPLE_NAME=$i
if [ $DEBUG_LEVEL -gt 2 ]
	then
       echo "DEBUG:  ${FILTER_DBSNP} $i ${DISEASE_CODE} hg19" >> LOG_DIR/log.txt
   fi
  ${FILTER_DBSNP} $i ${DISEASE_CODE} hg19
  cp ${i}_tier1_putative_mutation.xls $EXCEL_ROOT_DIR/${i}/.
  
  D_BAM_FILE=`ls ${BAM_DIR}/${SAMPLE_NAME}_D*.bam`
  G_BAM_FILE=`ls ${BAM_DIR}/${SAMPLE_NAME}_G*.bam`

  echo "${SAMPLE_NAME}_D|${D_BAM_FILE}|${G_BAM_FILE}" >${SAMPLE_NAME}_bam_file_info.lst
if [ $DEBUG_LEVEL -gt 2 ]
	then
  echo "${SCRIPT_DIR}/run_coding_indel_analysis_mod ${SAMPLE_NAME}  ${DISEASE_CODE} `pwd`/${SAMPLE_NAME}_bam_file_info.lst `pwd`" >> LOG_DIR/log.txt
fi
  ${SCRIPT_DIR}/run_coding_indel_analysis_mod ${SAMPLE_NAME}  ${DISEASE_CODE} `pwd`/${SAMPLE_NAME}_bam_file_info.lst `pwd`
done
