#!/bin/sh
## convert the format of TCGA-09-0365_somatic_snp.txt.clean in /tcga_next_gen/TCGA_WholeGenome/SubAnalysis/ME_WashU
## into a format similar to somatic_indel_sum.txt
## it has the following info
## Name    Chr     Pos     Type    Size    Coverage        Percent_alternative_allele      Chr_Allele      Alternative_Allele       P-value Text    reference_normal_count  reference_tumor_count   alternative_normal_count        alternative_tumor_count  dbSNP
## Gene    Sample  Repeat  Chr     Pos     #Unique #Indel  #Total  #IndelN #TotalN
## ADCY3   TCGA-13-0913    1       2       24907110        4       9       19      0       19      CATCACCAGG(ATGAAGATGGCGAGCATGGCCCAGG>-------------------------)TGTTCCTGGC

INPUT_SUB_DIR=$1  ##/tcga_next_gen/TCGA_WholeGenome/SubAnalysis/ME_WashU
BAD_SNP_FILE=$2  ## either the false positive ones to get rid of or the germline snps to be rid of in the somatic SNP analysis
OUTPUT_FILE=$3  ## this is the sub_sum.txt file that has all the frequency info
OUTPUT_4_Annotation=$4  ## this is the file of validate_somatic_loci.txt for running annotation
CHR_REGION=$5	 ##this region records chromosomal region for exons of interest. It has the following information: geneName, exon, acc, chr, from, to
ORIGINAL_SNP_FILE_EXT=$6  ## extension of somatic snp file prior to filtering. e.g. ${sample}_somatic_snp.txt. U
## ORIGINAL_SNP_FILE_EXT=somatic_snp.txt
SNP_FILE_EXT=${ORIGINAL_SNP_FILE_EXT}.clean  ## things appended to the end of each SNPs

NEXTGEN_BIN_ROOT_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess
PERLSRC_DIR=${NEXTGEN_BIN_ROOT_DIR}/perlsrc

FIND_SUB=${PERLSRC_DIR}/FindSub.pl
FIND_OVERLAP=${PERLSRC_DIR}/find_overlap_interval.pl

if [ $DEBUG_LEVEL -gt 0 ]
then
echo "region=${CHR_REGION}"
fi

if test -s ${OUTPUT_FILE}
then
  rm ${OUTPUT_FILE}
fi

if test -s ${OUTPUT_4_Annotation}
then
  rm ${OUTPUT_4_Annotation}
fi

cd ${INPUT_SUB_DIR}

cd ${INPUT_SUB_DIR}  ## directories where all the *_somatic_snp.txt reside
ls *_${ORIGINAL_SNP_FILE_EXT} >somatic_snp_file.lst
if test ! -s somatic_snp_file.lst
then
if [ $DEBUG_LEVEL -gt 0 ]
then
  echo "Fail to find files _${ORIGINAL_SNP_FILE_EXT} in subdirectory ${INPUT_SUB_DIR}"
  fi
exit 1
fi

## this step is not necessary as it may eliminate the possibility of having to avoid getting the consistent format## this is a legacy of previous effort to identify putative variations

for i in `cat somatic_snp_file.lst`; do
  cp ${i} ${i}.clean
done

ls *_${SNP_FILE_EXT}  |cut -f1 -d "_" >t_sample.lst
SAMPLE_LIST=t_sample.lst
echo "Gene	Sample	Repeat	Chr	Pos	#Unique	#Indel	#Total	#IndelN	#TotalN	RefAllele	MutantAllele	Flanking" >${OUTPUT_FILE}

for i in `cat ${SAMPLE_LIST}`; do
  SNP_FILE=${INPUT_SUB_DIR}/${i}_${SNP_FILE_EXT}
  if test ! -s ${SNP_FILE}
  then
if [ $DEBUG_LEVEL -gt 0 ]
then
    echo "Fail to find input file ${SNP_FILE}"
  fi
else
    cut -f2 ${SNP_FILE} |sort -u >t_chr.lst
    for j in `cat t_chr.lst`; do
      echo ${j} >t.lst
      ${FIND_SUB} -i ${SNP_FILE} -c t.lst -t 1 -d '\t' -n 1 -o t.lst.out
      if test -s t.lst.out
      then
        awk '{printf("%s\t%s\n", $3, $3)}' t.lst.out >t_interval.lst
        ${FIND_OVERLAP} -i t_interval.lst -d ${CHR_REGION}/${j}_region.txt -o k2 -e err.log -m 0 -n 4 -q1 -b10
        if test -s k2
        then
          cut -f1,3 k2 |sort +0 -1 >x
          sort +2 -3 t.lst.out >t.lst.out.sort
          join -1 1 -2 3 x t.lst.out.sort >k
          if test -s k
          then
            sed /^/s//${i}" "/ k|awk '{printf("%s\t%s\t0\t%s\t%s\t%ld\t%ld\t%ld\t%ld\t%ld\t%s\t%s\t%s\n", $3, $1, $5, $2, $17, $17, $15+$17, $16, $16+$14, $10, $11, $13)}' >>${OUTPUT_FILE}
          fi
        fi
      fi
    done
  fi
done


if test -s ${OUTPUT_FILE}
then
  cat ${OUTPUT_FILE}  |grep -v Gene |awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $4, $5, $5, $11, $11, $12, $2)}' |sort +0 -1 |sed /chr/s/// |sort +0 -1 >${OUTPUT_4_Annotation}
fi
