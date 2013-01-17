#!/bin/sh

## this program was designed to target tiers 2,3,4 mutations
## tier2=conserved-region; tier3=non-repetitive; tier4=others  ## repeatmasker
SAMPLE_NAME=$1
PROCESS_DIR=$2
BAD_SNP_FILE=$3  ## use no_bad_snp_file for default
DEBUG_LEVEL=$4
$LOG=$5

PUTATIVE_SOMATIC_MUTATION_FILE=${PROCESS_DIR}/${SAMPLE_NAME}_putative_somatic_mutation.txt  ## file that includes all the mutations mapped to the genes

PROCESSED_SNP_FILE=${PROCESS_DIR}/${SAMPLE_NAME}_processed_snp.lst

PERL_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess/perlsrc
FIND_SUB=${PERL_DIR}/FindSub.pl

SCRIPT_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess
FIND_GENE_PROGRAM=${SCRIPT_DIR}/find_gene4sub


## tier2 tracks
GENOME_INFO_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess
TFBS_REGION=${GENOME_INFO_DIR}/TFBS_ConsSites  ## transcription binding sites
CPG_REGION=${GENOME_INFO_DIR}/CpG_REGION
CONSERVED_REGION=${GENOME_INFO_DIR}/HUMAN_MOUSE_CONSERVATION
REPEAT_REGION=${GENOME_INFO_DIR}/REPEATMASKER


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
    echo "Fail to find program ${FIND_GENE_PROGRAM}" >> $LOG
  fi
exit 1
fi

if test ! -d ${TFBS_REGION}
then
  if [ $DEBUG_LEVEL -gt 0 ]
then
echo "Fail to find transcription binding sites direcotry ${TFBS_REGION} in Annotate_lowtier" >> $LOG
fi
exit 1
fi

if test ! -d ${CPG_REGION}
then
  if [ $DEBUG_LEVEL -gt 0 ]
then
echo "Fail to find CpG island direcotry ${CPG_REGION} in Annotate_lowtier" >> $LOG
fi
  exit 1
fi

if test ! -d ${CONSERVED_REGION}
then
if [ $DEBUG_LEVEL -gt 0 ]
then
  echo "Fail to find conserved direcotry ${CONSERVED_REGION} in Annotate_lowtier" >> $LOG
fi
  exit 1
fi

if test ! -d ${REPEAT_REGION}
then
if [ $DEBUG_LEVEL -gt 0 ]
then
  echo "Fail to find repeat direcotry ${REPEAT_REGION} in Annotate_lowtier" >> $LOG
  fi
exit 1
fi

if test ! -d ${PROCESS_DIR}
then
if [ $DEBUG_LEVEL -gt 0 ]
then
  echo "Fail to find diretory ${PROCESS_DIR} in Annotate_lowtier" >> $LOG
fi
  exit 1
fi



cat ${HIGH_CONFIDENCE_SOMATIC_SNP_FILE} ${HIGH_CONFIDENCE_SOMATIC_SNP_FILE_REPEAT} ${LOW_CONFIDENCE_SOMATIC_SNP_FILE} ${LOW_CONFIDENCE_SOMATIC_SNP_FILE_REPEAT} |sort -u >${SAMPLE_NAME}_somatic_snp.txt

if test -s ck_file.lst
then
  rm ck_file.lst
fi
echo ${HIGH_CONFIDENCE_SOMATIC_SNP_FILE} >>ck_file.lst
echo ${HIGH_CONFIDENCE_SOMATIC_SNP_FILE_REPEAT} >>ck_file.lst
echo ${LOW_CONFIDENCE_SOMATIC_SNP_FILE} >>ck_file.lst
echo ${LOW_CONFIDENCE_SOMATIC_SNP_FILE_REPEAT} >>ck_file.lst

## keep track of which snp has already been processed
if test -s ${PROCESSED_SNP_FILE}
then
  rm ${PROCESSED_SNP_FILE}
fi
if test -s ${BAD_SNP_FILE}
then
  cat ${BAD_SNP_FILE} >>${PROCESSED_SNP_FILE}
fi

if test -s ${PUTATIVE_SOMATIC_MUTATION_FILE}
then
  awk '{printf("%s.%s\n", $4, $5)}' ${PUTATIVE_SOMATIC_MUTATION_FILE} >>${PROCESSED_SNP_FILE}
fi


## processing the tier2 info
echo "2|${TFBS_REGION}" >dir.lst
echo "2|${CPG_REGION}" >>dir.lst
echo "2|${CONSERVED_REGION}" >>dir.lst
echo "4|${REPEAT_REGION}" >>dir.lst
echo "3|no_dir" >>dir.lst

TIER2_SUMMARY_FILE=${PROCESS_DIR}/${SAMPLE_NAME}_tier2_putative_mutation.txt
TIER3_SUMMARY_FILE=${PROCESS_DIR}/${SAMPLE_NAME}_tier3_putative_mutation.txt
TIER4_SUMMARY_FILE=${PROCESS_DIR}/${SAMPLE_NAME}_tier4_putative_mutation.txt

if test -s ${TIER2_SUMMARY_FILE}
then
  rm ${TIER2_SUMMARY_FILE}
fi

if test -s ${TIER3_SUMMARY_FILE}
then
  rm ${TIER3_SUMMARY_FILE}
fi

if test -s ${TIER4_SUMMARY_FILE}
then
  rm ${TIER4_SUMMARY_FILE}
fi

for k in `cat dir.lst`; do
  cd ${PROCESS_DIR}
  tier=`echo $k|cut -f1 -d"|"`
  ck_dir=`echo $k|cut -f2 -d"|"`
  if [ $DEBUG_LEVEL -gt 0 ]
then
echo "tier=$tier; and ck_dir=$ck_dir" >> $LOG
fi
  if test $tier = 2
  then
    OUTPUT_SUMMARY_FILE=${TIER2_SUMMARY_FILE}
  else
    if test $tier = 3
    then
      OUTPUT_SUMMARY_FILE=${TIER3_SUMMARY_FILE}
    else
      OUTPUT_SUMMARY_FILE=${TIER4_SUMMARY_FILE}
    fi
  fi


  if test -s ${SAMPLE_NAME}_somatic_snp.txt
  then
    if test -s ${PROCESSED_SNP_FILE}
    then
      ${FIND_SUB} -i ${SAMPLE_NAME}_somatic_snp.txt -c ${PROCESSED_SNP_FILE} -t 0 -d '\t' -n 0 -o ${SAMPLE_NAME}_somatic_snp.txt.clean
      mv ${SAMPLE_NAME}_somatic_snp.txt.clean ${SAMPLE_NAME}_somatic_snp.txt
    fi
  fi

  if test -s ${SAMPLE_NAME}_somatic_snp.txt
  then
    if test -d temp
    then
      rm -rf temp
    fi
    mkdir temp
    mkdir temp/SOMATIC_RESULTS

    if test -d ${ck_dir}
    then
      cp ${SAMPLE_NAME}_somatic_snp.txt temp/.
      cd temp
      ${FIND_GENE_PROGRAM} `pwd` ${BAD_SNP_FILE} `pwd`/SOMATIC_RESULTS/sub_sum.txt `pwd`/SOMATIC_RESULTS/validate_somatic_loci.txt ${ck_dir} somatic_snp.txt
    else
      cd temp
      cat ../${SAMPLE_NAME}_somatic_snp.txt |sed /^/s//"NoRepeat	${SAMPLE_NAME}	"/ |awk '{printf("%s\t%s\t0\t%s\t%s\t%ld\t%ld\t%ld\t%ld\t%ld\t%s\t%s\t%s\n", $1, $2, $4, $5, $17, $17, $15+$17, $16, $16+$14, $10, $11, $13)}' >SOMATIC_RESULTS/sub_sum.txt
    fi


    if test -s SOMATIC_RESULTS/sub_sum.txt
    then
      cd ${PROCESS_DIR}
      awk '{printf("%s.%ld\t", $4, $5); print $0}' temp/SOMATIC_RESULTS/sub_sum.txt  >temp/SOMATIC_RESULTS/sub_sum.txt.mod

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
            cut -f2-100 c_snp.lst.out |sed /^/s//"${label_val}	"/|awk '{printf("%s\t%s\t%s\t%s\t%s", $2, $1, $3, $5, $6); for(i=8; i<=NF; ++i) printf("\t%s", $i); printf("\n")}'>>${OUTPUT_SUMMARY_FILE}
            cut -f1 c_snp.lst.out >>${PROCESSED_SNP_FILE}
          fi
        fi
      done
    fi
  fi
done

if test -s ${TIER2_SUMMARY_FILE}
then
  sort +1 -2 ${TIER2_SUMMARY_FILE} >${TIER2_SUMMARY_FILE}.sort
  echo "Class	SJQuality	Sample	Chr	WU_HG18_Pos	#Mutant_In_Tumor	#Total_In_Tumor	#Mutant_In_Normal	#Total_In_Normal	MutantAllele	ReferenceAllele	Flanking" >${TIER2_SUMMARY_FILE}
  cat ${TIER2_SUMMARY_FILE}.sort >>${TIER2_SUMMARY_FILE}
  rm ${TIER2_SUMMARY_FILE}.sort
fi

if test -s ${TIER3_SUMMARY_FILE}
then
  sort +1 -2 ${TIER3_SUMMARY_FILE} >${TIER3_SUMMARY_FILE}.sort
  echo "Class	SJQuality	Sample	Chr	WU_HG18_Pos	#Mutant_In_Tumor	#Total_In_Tumor	#Mutant_In_Normal	#Total_In_Normal	MutantAllele	ReferenceAllele	Flanking" >${TIER3_SUMMARY_FILE}
  cat ${TIER3_SUMMARY_FILE}.sort >>${TIER3_SUMMARY_FILE}
  rm ${TIER3_SUMMARY_FILE}.sort
fi

if test -s ${TIER4_SUMMARY_FILE}
then
  sort +1 -2 ${TIER4_SUMMARY_FILE} >${TIER4_SUMMARY_FILE}.sort
  echo "Class	SJQuality	Sample	Chr	WU_HG18_Pos	#Mutant_In_Tumor	#Total_In_Tumor	#Mutant_In_Normal	#Total_In_Normal	MutantAllele	ReferenceAllele	Flanking" >${TIER4_SUMMARY_FILE}
  cat ${TIER4_SUMMARY_FILE}.sort >>${TIER4_SUMMARY_FILE}
  rm ${TIER4_SUMMARY_FILE}.sort
fi
