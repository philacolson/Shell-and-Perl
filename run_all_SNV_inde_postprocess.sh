#!/bin/sh
INPUT_SAMPLE_LIST=$1 #list of chromosomes to test
DISEASE_CODE=$2
##BAM_DIR=/nfs_exports/genomes/1/PCGP/BucketRaw/${DISEASE_CODE}
BAM_DIR=/user/songliu/u2/group/Qiang/Exome/bladder_output/${DISEASE_CODE}
DEBUG_LEVEL=$3
#Debug level can be 1-3, increasing with verbosity.  As a general rule, the grades pertain to the following: 1 just logs fatal events and ones that otherwise prematurely end the code, 2 mentions anything that deviates from expected code/error handler code is run, or other noteworthly events.  3 logs variable names before and after programs are run, and upon initialization.  Higher numbers post their description + any number's description lower than their own.  For example, a debug level of 2 would post every log statement that is not a level 3 log event.  It would post 1's, and 2's.

##SCRIPT_DIR=/nfs_exports/apps/internal/scripts
SCRIPT_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess


FILTER_DBSNP=$SCRIPT_DIR/filter_low_quality_dbSNP_4_tier1 


##EXCEL_ROOT_DIR=/nfs_exports/linux-file1/home/naevegrp/jzhang2/NextGen/PCGP
EXCEL_ROOT_DIR=/user/songliu/u2/group/Qiang/Exome/bladder_output/EXCEL
LOG_DIR = EXCEL_ROOT_DIR

echo $DEBUG_LEVEL >> LOG_DIR/debug.lvl
`cat $LOG`
for i in `cat ${INPUT_SAMPLE_LIST}`; do
  if test ! -d ${EXCEL_ROOT_DIR}/${i} #if there is no directory for the SNP
  then
	if [ $DEBUG_LEVEL -gt 1 ]
	then 
  echo "INFO: no dir for ${EXCEL_ROOT_DIR}/${i}, making directory." >> LOG_DIR/log.txt
	fi
    mkdir -p ${EXCEL_ROOT_DIR}/${i}
    #exit 1
  fi
done

if [ $DEBUG_LEVEL -gt 2 ]
then 
echo "DEBUG: created EXCEL directories." >>LOG_DIR/log.txt
fi

for i in `cat ${INPUT_SAMPLE_LIST}`; do
  ##echo $i
  cp ../${i}_G_bam.out .
  cp ../${i}_bam_low.out ${i}_bam_low.out
  cp ../${i}_bam_high_20.out ${i}_bam_high_20.out
if [ $DEBUG_LEVEL -gt 2 ]
	then echo "DEBUG:  ${SCRIPT_DIR}/run_somatic_mutation_analysis ${i} no_false_snp" >> LOG_DIR/log.txt
	fi
  #hq
  D_BAM=`ls ${BAM_DIR}/${i}_D*.bam`
  G_BAM=`ls ${BAM_DIR}/${i}_G*.bam`
  ${SCRIPT_DIR}/run_somatic_mutation_analysis ${i} no_false_snp ${D_BAM} ${G_BAM} ${LOG_DIR}
done
if [ $DEBUG_LEVEL -gt 2 ]
	then
echo "DEBUG:  ${SCRIPT_DIR}/run_check_snv_script ${INPUT_SAMPLE_LIST} ${DISEASE_CODE}" >> LOG_DIR/log.txt
fi

${SCRIPT_DIR}/run_check_snv_script ${INPUT_SAMPLE_LIST} ${DISEASE_CODE}

for i in `cat ${INPUT_SAMPLE_LIST}`; do
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
