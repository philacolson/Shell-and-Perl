#!/bin/sh

## this program checks the location where variations are found. It checks the files under directory /TCGA/nextgensupport/non_specific on lpgws506 to see whether any site has been covered by MM (either MM or MM,PM)
## it checks the size of read_length/2 as the flanking region
## the output file records the number of MM found by this approach
## this program can only be run on lpgws506
## /h1/edmonsom/bin/cov_decode -chr 19 -size 50 -start 35006400  -length 100

INPUT_VAR_FILE=$1   ##in the format of chr10:101721814. see H25_snp.lst
READ_LENGTH=$2
START_UPSTREAM=$3	## start the flanking base analysis upstream of the site
PROCESS_DIR=$4
OUTPUT_FILE=$5

NEXTGEN_BIN_ROOT_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess
PERLSRC_DIR=${NEXTGEN_BIN_ROOT_DIR}
## PROGRAM=/h1/edmonsom/bin/cov_decode
PROGRAM=${PERLSRC_DIR}/cov_decode

if test ! -d ${PROCESS_DIR}
then
  echo "Fail to find the process directory for storing the temporary file"
  exit 1
fi

if test -s ${OUTPUT_FILE}
then
  rm ${OUTPUT_FILE}
fi
touch ${OUTPUT_FILE}

TEMP_CVG_FILE=${PROCESS_DIR}/cvg.tmp.out

for i in `cat ${INPUT_VAR_FILE}`; do
  chr=`echo $i |cut -f1 -d":" |sed /^chr/s///`
  pos=`echo $i |cut -f2 -d":"`
  if test ${START_UPSTREAM} = 0
  then
    from=`echo "${READ_LENGTH} $pos" |awk '{printf("%ld", $2-$1/2)}'`
    to=`echo "${READ_LENGTH} $pos" |awk '{printf("%ld", $2+$1/2)}'`
  else
    from=`echo "${READ_LENGTH} $pos" |awk '{printf("%ld", $2-$1)}'`
    to=$pos
  fi
  ${PROGRAM} -chr $chr -size ${READ_LENGTH} -start $from -length ${READ_LENGTH} >${TEMP_CVG_FILE}
  if test -s ${TEMP_CVG_FILE}
  then
    count=`grep MM ${TEMP_CVG_FILE} |wc |awk '{printf("%ld", $1)}'`
    echo "$i $count" |awk '{printf("%s\t%s\n", $1, $2)}' >>${OUTPUT_FILE}
    rm ${TEMP_CVG_FILE}
  fi
done
