#!/bin/sh
INPUT_FILE=$1
MIN_GOOD_READ=$2
MIN_DISTANCE_4_RREPLICATE=$3
USE_LOW_THRESHOLD=$4

## bash-3.2$ head all_blat.lst
## chr10_112686921_blat.out
## chr10_42331711_blat.out
## chr10_42333383_blat.out
## chr10_42333393_blat.out
## chr10_42333448_blat.out
## chr10_42333449_blat.out
## chr10_42333526_blat.out
## chr10_42334142_blat.out

NEXTGEN_BIN_ROOT_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess
PERL_DIR=${NEXTGEN_BIN_ROOT_DIR}/perlsrc

## EXON_DIR=${NEXTGEN_BIN_ROOT_DIR}/GENE_EXON_REGION
EXON_DIR=/nfs_exports/apps/gnu-apps/NextGen/nextgensupport/hg19/GENE_EXON_REGION
FIND_SUB=${PERL_DIR}/FindSub.pl
CHECK_REPEAT=${PERL_DIR}/find_repeat_in_blat.pl

if test ! -d ${EXON_DIR}
then
  echo "Fail to find exon dir ${EXON_DIR}"
  exit 1
fi

for i in `cat ${INPUT_FILE}`; do
  chr=`echo $i |cut -f1 -d"_"`
  pos=`echo $i |cut -f2 -d"_"`

  echo $chr
  echo $pos
  echo "awk '{if($pos >=\$5 && $pos <=\$6) print \$0}' ${EXON_DIR}/${chr}_region.txt" >awk.x
  chmod ugo+x awk.x
  ./awk.x >exon.out
  if test -s current_gene.lst.out
  then
    rm current_gene.lst.out
  fi
  if test -s exon.out
  then
    cut -f1 exon.out |sort -u >current_gene.lst
    ${FIND_SUB} -i ${EXON_DIR}/${chr}_region.txt -c current_gene.lst -t 1 -d '\t' -n 0 -o current_gene.lst.out
  fi
    
  if test -s current_gene.lst.out
  then
    echo "${CHECK_REPEAT} -i $i -d ${MIN_DISTANCE_4_RREPLICATE} -c $chr -p $pos -b current_gene.lst.out -o ${i}.good -e  ${i}.bad -m ${USE_LOW_THRESHOLD} -s ${USE_LOW_THRESHOLD} -g ${MIN_GOOD_READ}"
    ${CHECK_REPEAT} -i $i -d ${MIN_DISTANCE_4_RREPLICATE} -c $chr -p $pos -b current_gene.lst.out -o ${i}.good -e  ${i}.bad -m ${USE_LOW_THRESHOLD} -s ${USE_LOW_THRESHOLD} -g ${MIN_GOOD_READ}
  else
    echo "${CHECK_REPEAT} -i $i -d ${MIN_DISTANCE_4_RREPLICATE} -c $chr -p $pos -o ${i}.good -e ${i}.bad -m ${USE_LOW_THRESHOLD} -s ${USE_LOW_THRESHOLD} -g ${MIN_GOOD_READ}"
    ${CHECK_REPEAT} -i $i -d ${MIN_DISTANCE_4_RREPLICATE} -c $chr -p $pos -o ${i}.good -e ${i}.bad -m ${USE_LOW_THRESHOLD} -s ${USE_LOW_THRESHOLD} -g ${MIN_GOOD_READ}
  fi

done
