#!/bin/sh

## this program is designed to create geneName2geneID.map file
## it also adds the configuration files for processing


INPUT_FILE=$1  ##validate_somatic_loci.txt
PROCESS_DIR=$2  ## ROOT directory for running the annotation

NEXTGEN_BIN_ROOT_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess
PERL_SRC=${NEXTGEN_BIN_ROOT_DIR}/perlsrc
FIND_SUB=${PERL_SRC}/FindSub.pl

SRC4CONFIG=${NEXTGEN_BIN_ROOT_DIR}/nextgensupport2
NCBI_DIR=${NEXTGEN_BIN_ROOT_DIR}/nextgensupport2/NCBI
GENE_INFO_FILE=${NCBI_DIR}/gene_info_human
RESCUE_GENEMAP_FILE=${NCBI_DIR}/all_manual_geneName2geneID.map  ## this is a list of manually curated geneName2geneID.map to avoid manually look up the records that is not in the current gene_info_human

COPY_CONFIG_SCRIPT=${NEXTGEN_BIN_ROOT_DIR}/copy_config_script

if test ! -d ${PROCESS_DIR}
then
  echo "Fail to find direcotry ${PROCESS_DIR}"
  exit 1
fi

if test ! -s ${INPUT_FILE}
then
  echo "Fail to find validate_somatic_loci.txt file ${INPUT_FILE}"
  exit 1
fi

if test ! -s ${GENE_INFO_FILE}
then
  echo "Fail to find gene info file ${GENE_INFO_FILE}"
  exit 1
fi

if test ! -s ${RESCUE_GENEMAP_FILE}
then
  echo "Fail to find rescue file ${RESCUE_GENEMAP_FILE}"
  exit 1
fi

cd ${PROCESS_DIR}
${COPY_CONFIG_SCRIPT} ${SRC4CONFIG} ${PROCESS_DIR}

if test -s all_gene.lst
then
  rm all_gene.lst
fi
cut -f1 ${INPUT_FILE}  |sort -u >all_gene.lst
if test ! -s all_gene.lst
then
  echo "Fail to find file all_gene.lst"
  exit 1
fi

## create the file geneName2geneID.map
${FIND_SUB} -i ${GENE_INFO_FILE} -c all_gene.lst -t 1 -d '\t' -n 2 -o gene.lst.out
if test -s geneName2geneID.map
then
  rm geneName2geneID.map
fi

if test -s gene.lst.out
then
  awk '{printf("%s\t%s\n", $3, $2)}' gene.lst.out >geneName2geneID.map
fi

cut -f3 gene.lst.out |sort -u >x
cat x all_gene.lst |sort |uniq -u >failed_gene.lst
if test -s failed_gene.lst
then
  ${FIND_SUB} -i ${RESCUE_GENEMAP_FILE} -c failed_gene.lst -t 1 -d '\t' -n 0 -o failed_gene.lst.out
  cat failed_gene.lst.out >>geneName2geneID.map
  cut -f1 failed_gene.lst.out |sort -u >failed_gene.lst.out.match
  cat failed_gene.lst failed_gene.lst.out.match |sort |uniq -u >u
  if test -s u
  then
    mv u failed_gene.lst
    echo "Check failed_gene.lst and look up the gene ID by going to Entrez Gene"
    exit 1
  fi
fi

## manually find out the gene's new name and store that in failed_gene2new_gene.map
cut -f8 ${INPUT_FILE}  |sort -u >sample.lst

echo "geneName2geneID.map and sample.lst files are created"
