#!/bin/sh

## this program goes through one-chr at a time to add coverage to the indel analysis output
DONOR_NAME=$1  ## input file are expected as ${DONOR_NAME}_tumor.txt and ${DONOR_NAME}_normal.txt
DISEASE_CODE=$2
INPUT_TumorNormal_DIR=$3
BAM4TUMOR=$4   ## put a fake name if you would like it to get the bam file from default
BAM4NORMAL=$5   ## put a fake name if you would like it to get the bam file from default

PROCESS_DIR=$INPUT_TumorNormal_DIR  ## this will be where the output *_mod2 file will be stored

TUMOR_EXT=D  ## this may need to be changed

NEXTGEN_BIN_ROOT_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess
PERLSRC_DIR=${NEXTGEN_BIN_ROOT_DIR}/perlsrc
FIND_SUB=${PERLSRC_DIR}/FindSub.pl
COVERAGE_CODE=${NEXTGEN_BIN_ROOT_DIR}/coverage.sh
BAM_DIR=/nfs_exports/genomes/1/PCGP/BucketRaw/${DISEASE_CODE}


ls -lt ${BAM4TUMOR}
if test -s ${BAM4TUMOR}
then
  echo "find ${BAM4TUMOR}"
else
  ls -lt ${BAM4TUMOR}
  echo "no bam file ${BAM4TUMOR}"
  if test ! -d ${BAM_DIR}
  then
    echo "Fail to find bucketraw dir ${BAM_DIR}"
    exit 1
  fi
  count=`ls ${BAM_DIR}/${DONOR_NAME}_D-*.bam |wc|awk '{printf("%s", $1)}'`
  if test $count != 1
  then
    echo "fail to figure out bam file name. tumor_count=$count in ${BAM_DIR}/${DONOR_NAME}_D-*.bam"
    exit 1
  fi
  BAM4TUMOR=`ls ${BAM_DIR}/${DONOR_NAME}_D-*.bam`
fi

ls -lt ${BAM4NORMAL}
if test -s ${BAM4NORMAL}
then
  echo "find ${BAM4NORMAL}"
else
  echo "no bam file ${BAM4NORMAL}"
  if test ! -d ${BAM_DIR}
  then
    echo "Fail to find bucketraw dir ${BAM_DIR}"
    exit 1
  fi
  count=`ls ${BAM_DIR}/${DONOR_NAME}_G-*.bam |wc|awk '{printf("%s", $1)}'`
  if test $count != 1
  then
    echo "fail to figure out bam file name. tumor_count=$count in ${BAM_DIR}/${DONOR_NAME}_G-*.bam"
    exit 1
  fi
  BAM4NORMAL=`ls ${BAM_DIR}/${DONOR_NAME}_G-*.bam`
fi




TUMOR_INPUT_FILE=${INPUT_TumorNormal_DIR}/${DONOR_NAME}_tumor.txt

if test ! -s ${TUMOR_INPUT_FILE}
then
  echo "Fail to find the tumor file ${TUMOR_INPUT_FILE}"
  exit 1
fi

cd ${INPUT_TumorNormal_DIR}
awk '{printf("%s|%s\n", $2, $3)}' ${DONOR_NAME}_tumor.txt |sort -u >tumor_pos.lst
if test -s tumor_cvg.txt
then
  rm tumor_cvg.txt
fi

if test -s normal_cvg.txt
then
  rm normal_cvg.txt
fi

for i in `cat tumor_pos.lst`; do
  chr=`echo $i |cut -f1 -d"|"`
  pos=`echo $i |cut -f2 -d"|"`

  ${COVERAGE_CODE} ${BAM4TUMOR} $chr $pos $pos tumor.out
  if test -s tumor.out
  then
     sed /\,/s//\\\t/ tumor.out |sed /\,/s//\\\t/ |sed /^/s//${chr}_/ >>tumor_cvg.txt
  fi

  ${COVERAGE_CODE} ${BAM4NORMAL} $chr $pos $pos normal.out
  if test -s normal.out
  then
     sed /\,/s//\\\t/ normal.out |sed /\,/s//\\\t/ |sed /^/s//${chr}_/ >>normal_cvg.txt
  fi
done

sort +0 -1 tumor_cvg.txt >tumor_cvg.txt.sort
mv tumor_cvg.txt.sort tumor_cvg.txt

sort +0 -1 normal_cvg.txt >normal_cvg.txt.sort
mv normal_cvg.txt.sort normal_cvg.txt

awk '{printf("%s_%s\t", $2, $3); print $0}' ${DONOR_NAME}_tumor.txt |sort +0 -1 >y
join -1 1 -2 1 y tumor_cvg.txt |awk '{for(i=1; i<NF; ++i) printf("%s\t", $i); printf("%s\n", $NF)}' |sort +0 -1>${DONOR_NAME}_tumor.txt.mod
join -1 1 -2 1 ${DONOR_NAME}_tumor.txt.mod  normal_cvg.txt  |awk '{for(i=2; i<NF; ++i) printf("%s\t", $i); printf("%s\n", $NF)}' |sort +0 -1 -n -r >${DONOR_NAME}_tumor.txt.mod2
mv ${DONOR_NAME}_tumor.txt.mod2  ${DONOR_NAME}_tumor_mod.txt
rm y

disease_tag=tumor
echo "Chr	Pos	IndelSize	RepeatSize	PairCvg	TumorCvg	NormalCvg	Flanking">x
awk '{printf("%s\t%s\t%ld\t%ld\t", $2, $3, $1, $5-$4+1); printf("%s\t%s\t%s\t%s\n", $11, $13, $14, $9)}' ${DONOR_NAME}_${disease_tag}_mod.txt >>x
echo "AlleleCvg UniqCvg" >y
cut -f10 ${DONOR_NAME}_${disease_tag}_mod.txt |sed /\;/s//\\\t/g |cut -f1,2 |sed /HighQ\:/s///g |sed /HQU\:/s///g >>y
paste x y |awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\t", $1, $2, $3, $4, $9, $10); printf("%s\t%s\t%s\t%s\n", $5, $6, $7, $8)}' >${DONOR_NAME}_${disease_tag}_mod2.txt
