#!/bin/sh
## REPORT_FILE=/nfs_exports/genomes/1/projects/RNASEQ/PCGP/BucketIntermediate/SnpDetect/snp_check_report.txt
## BAM_FILE=/nfs_exports/genomes/1/projects/RNASEQ/PCGP/BucketRaw/SJTALL012_D-TB-09-4637P2T-reordered.bam
## CHR=1
## POS=115058052
## ALLELE=G

## BAM_FILE=/nfs_exports/genomes/1/projects/RNASEQ/PCGP/BucketRaw/SJTALL002_D-TB-05-2043-reordered.bam
## CHR=5
## POS=10311497
## ALLELE=A

## this program checks for the validity of SNVs by the following processes:
## a) are all the variants derived from clone?
## b) are the variants reside outside of 10bp from the beginning of the alignment
## c) any variants that failed to be mapped by blat to the target region
## d) any variants that turn out to be non-specific based on blat search

BAM_FILE=$1
CHR=$2
POS=$3
ALLELE=$4
REPORT_FILE=$5
IGNORE_STRAND_CHECK=$6  ## this is applicable for those that have problems with uneven forward/reverse coverage
MATCHING_GERMLINE_BAM_FILE=$7  ## there could be false negative germline calls in the original bam_20.out etc. This is an optional field and can be set as no_germline_bam_file. This flag was needed only when there is an issue related to have a low-frequency germline allele being missed from SNPFind
## germline check needs to have a couple of options: for those that the matching normal may have tumor contamination
## record matches rs number
## three posibilities for germline presence: a) tumor-contamination; b) repeats; and c) germline. b) and c) need to be rid of while
## a) needs to be retained. a) and b) may share the same profile of having low frequency
CHECK_GERMLINE_ALLELE_COUNT=$8

NEXTGEN_BIN_ROOT_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess
HG_FASTA_DIR=/user/songliu/u2/group/Qiang/Exome/reference/fasta
PERL_DIR=${NEXTGEN_BIN_ROOT_DIR}/perlsrc
SCRIPT_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess
JZ_DIR=${NEXTGEN_BIN_ROOT_DIR}
CREATE_INDEX_SCRIPT=${JZ_DIR}/CreateIndex4Fasta.scr
GET_FASTA_SCRIPT=${JZ_DIR}/FastaGetData1.prl


FILTER_CLONE=${PERL_DIR}/filter_clone_from_sam_output.pl
SAM2FASTQ=${PERL_DIR}/create_fastq_from_tweakSam.pl
RUN_BLAT_PROGRAM=${SCRIPT_DIR}/run_blat_search
## TRIM_SEQ=${SCRIPT_DIR}/trim_fastq   ## this is the program that converts fastq to fasta + doing trimming of very low quality bases
TRIM_SEQ=${SCRIPT_DIR}/trim_fastq_with_threshold
CHECK_BLAT_REPEAT=${SCRIPT_DIR}/check_all_blat_result


##  /nfs_exports/genomes/1/projects/RNASEQ/PCGP/BucketRaw/SJTALL002_D-TB-05-2043-reordered.bam 5:10311397-10311597 |/nfs_exports/apps/gnu-apps/NextGen/mrbin/java_sj_picard.sh TweakSam -G 1024 -s SILENT -p 10311497 -a A

## the fasta extracted are all mapped to the + orientation of the reads. So the double-strand check is meaningless here

echo "samtools view -b  ${BAM_FILE} ${CHR}:${POS}-${POS} |/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess/java_sj_picard.sh TweakSam -G 1024 -s SILENT -p ${POS} -a ${ALLELE} >test.out"
samtools view -b  ${BAM_FILE} ${CHR}:${POS}-${POS} |/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess/java_sj_picard.sh TweakSam -G 1024 -s SILENT -p ${POS} -a ${ALLELE} >test.out

MIN_GERMLINE_FREQ=10
MAX_GERMLINE_TOLERABLE_FREQ=20
if test ! -s test.out
then
  echo "Fail to find anything"
  exit 1
fi

grep -v ^@ test.out >test.out.clean

if test ! -s test.out.clean
then
  echo "Fail to find reads"
  echo "$CHR	${POS}	NoReads" >>${REPORT_FILE}
  exit 1
fi

${FILTER_CLONE} -i test.out.clean -o test.out.clean2 -c test.out.clone
if test ! -s test.out.clean2
then
  echo "$CHR	$POS	AllReadsClone" >>${REPORT_FILE}
  exit 1
fi
## mv test.out.clean2 test.out.clean


## need to have at least one non-clone count
count=`cut -f1 test.out.clean2  |sort -u |wc |awk '{printf("%ld", $1)}'`
if test $count -lt 1
then
  echo "$CHR	${POS}	OneClone" >>${REPORT_FILE}
  exit 1
fi

clone_count=0
if test -s test.out.clone
then
  clone_count=`cut -f1 test.out.clone  |sort -u |wc |awk '{printf("%ld", $1)}'`
  echo "clone_count=$clone_count"
fi
unique_read_count=`echo "$count $clone_count" |awk '{printf("%ld", $1+$2)}'`
if test $unique_read_count -le 1
then
  echo "$CHR	${POS}	OneClone" >>${REPORT_FILE}
  exit 1
fi


## awk '{printf("@%s_%s\n%s\n+\n%s\n", $1, $9, $10, $11)}' test.out.clean >test.fastq
## decode the flag to get the strand info
${SAM2FASTQ} -i test.out.clean -o test.fastq

## possibility of causing trouble
${TRIM_SEQ} `pwd`/test.fastq 3


## evaluate whether there is any germline sequences that match the site
## this step was used to ensure that matching alleles in germline can be extracted in case SNPFind misses them
GERMLINE_SEQ_FILE=germline.test.fastq
GERMLINE_TRIM_SEQ_FILE=${GERMLINE_SEQ_FILE}.trim.fa
if test -s ${GERMLINE_SEQ_FILE}
then
  rm ${GERMLINE_SEQ_FILE}
fi

if test -s ${GERMLINE_TRIM_SEQ_FILE}
then
  rm ${GERMLINE_TRIM_SEQ_FILE}
fi

if test -s test.out.clean
then
  rm test.out.clean
fi

if test -s ${MATCHING_GERMLINE_BAM_FILE}
then
  echo "samtools view -b  ${MATCHING_GERMLINE_BAM_FILE} ${CHR}:${POS}-${POS} |/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess/java_sj_picard.sh TweakSam -G 1024 -s SILENT -p ${POS} -a ${ALLELE}|grep -v ^@ >test.out.clean"
  samtools view -b  ${MATCHING_GERMLINE_BAM_FILE} ${CHR}:${POS}-${POS} |/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess/java_sj_picard.sh TweakSam -G 1024 -s SILENT -p ${POS} -a ${ALLELE}|grep -v ^@ >test.out.clean
  if test -s test.out.clean
  then
##    awk '{printf("@%s_%s\n%s\n+\n%s\n", $1, $9, $10, $11)}' test.out.clean >${GERMLINE_SEQ_FILE}
    if test -s ${GERMLINE_SEQ_FILE}
    then
      rm ${GERMLINE_SEQ_FILE}
    fi
    ${SAM2FASTQ} -i test.out.clean -o ${GERMLINE_SEQ_FILE}

## need to check the quality of the germline reads to be sure that if all bases are bad, do not use that germline read
    if test -s ${GERMLINE_TRIM_SEQ_FILE}
    then
      rm ${GERMLINE_TRIM_SEQ_FILE}
    fi

    if test -s ${GERMLINE_TRIM_SEQ_FILE}.qual
    then
      rm ${GERMLINE_TRIM_SEQ_FILE}.qual
    fi
    ${TRIM_SEQ} `pwd`/${GERMLINE_SEQ_FILE} 3  ## create ${GERMLINE_TRIM_SEQ_FILE}. may need to set it to 3. but try 1 first
    if test -s ${GERMLINE_TRIM_SEQ_FILE}
    then
      max_len=`grep -v ">" ${GERMLINE_TRIM_SEQ_FILE} |awk '{printf("%ld\n", length$(1))}' |sort +0 -1 -n -r |head -n1`
      echo "max_len=$max_len"
      if test $max_len -le 30
      then
        rm ${GERMLINE_TRIM_SEQ_FILE}
      fi
    fi

    if test -s ${GERMLINE_TRIM_SEQ_FILE}
    then
      rm ${GERMLINE_TRIM_SEQ_FILE}
      rm ${GERMLINE_TRIM_SEQ_FILE}.qual
      ${TRIM_SEQ} `pwd`/${GERMLINE_SEQ_FILE} 1  ## create ${GERMLINE_TRIM_SEQ_FILE}. may need to set it to 3. but try 1 first
    fi
  fi
fi
## completed collecting the possible germline sequence

SNPDetector_DIR=${NEXTGEN_BIN_ROOT_DIR}/SNPdetector
BIN_DIR=${SNPDetector_DIR}/LINUX/src
PRT_SEQLOC=${BIN_DIR}/prt_seqloc
WRITE_FETCH=${BIN_DIR}/write_fetch_file
TSIM=${BIN_DIR}/tsim
PRT_FEATURE_LOC=${BIN_DIR}/prt_feature_loc
PRT_CHMAP_ID=${BIN_DIR}/prt_chmap_id
## FIND_POLY=${BIN_DIR}/find_poly
FIND_POLY=${SNPDetector_DIR}/LINUX/new_src/find_poly
FASTA2ASN=${BIN_DIR}/fasta2asn
ATTACH_ENT=${BIN_DIR}/attatch_ent
FILTER_CLONE_ALIGNMENT=${PERL_DIR}/filter_alignment_by_clone.pl
NEW_CODE_DIR=${SNPDetector_DIR}/LINUX/new_src
CK_QUALITY_DROP=${NEW_CODE_DIR}/check_NextGen_quality_drop
FILTER_ALIGN=${NEW_CODE_DIR}/filter_align

echo "${WRITE_FETCH} -d ${HG_FASTA_DIR}/.fa"
${WRITE_FETCH} -d ${HG_FASTA_DIR}/.fa
## /nfs_exports/apps/gnu-apps/NextGen/SNPdetector/LINUX/src/write_fetch_file -d /nfs_exports/apps/gnu-apps/NextGen/nextgensupport/hg18/.fa

echo "$CHR $POS" |awk '{printf("chr%s|%ld|%ld|+", $1, $2-100, $2+100)}' >seqloc.lst
${PRT_SEQLOC} -i seqloc.lst -d ./ -b T -l 50

if test ! -s chr${CHR}.seq
then
  echo "$CHR	${POS}	NoSEQ" >>${REPORT_FILE}
  exit 1
fi

## a quick and ugly way of creating a modified seq by subsitute base 101 with the mutant allele
head -n3 chr${CHR}.seq >chr${CHR}.seq.mod
head -n4 chr${CHR}.seq |tail -n1 |awk '{printf("%s\n", substr($1, 2, 49))}' |sed /^/s//${ALLELE}/ >>chr${CHR}.seq.mod
head -n5 chr${CHR}.seq |tail -n1 >>chr${CHR}.seq.mod


germline_count=0
germline_original_count=0
## test for presence of germline sequence
if test -s ${GERMLINE_TRIM_SEQ_FILE}
then
  echo ">chr${CHR}.seq" >job.lst
  grep ">" ${GERMLINE_TRIM_SEQ_FILE} |cut -f2 -d">" >>job.lst
  ${WRITE_FETCH} -d ./.seq -l ${GERMLINE_TRIM_SEQ_FILE}
  if test -s chr${CHR}.ent
  then
    rm chr${CHR}.ent
  fi
  ${TSIM} -i job.lst -f4 -P1 -zF
  if test -s chr${CHR}.ent
  then
    ${FILTER_ALIGN} -i chr${CHR}.ent -o chr${CHR}.ent -c 0 -s 1 -u 0 -x 0 -k 0
    ${PRT_CHMAP_ID} -i chr${CHR}.ent -o out2
    germline_original_count=0
    if test -s out2
    then
      awk '{if($6 <=99 && $7 >=101 && $9 >=1) print $0}' out2 >out3
      if test -s out3
      then
        germline_original_count=`wc out3|awk '{printf("%s", $1)}'`
        echo "germline_original_count=$germline_original_count"
      fi
    fi

## make sure that the presence of the non-reference allele in germline is not bogus
## ignore the quality threshold
    echo "${CK_QUALITY_DROP}  -i chr${CHR}.ent -q ${GERMLINE_TRIM_SEQ_FILE}.qual -m 2 -a ${ALLELE} -o xxx_report.out"
    ${CK_QUALITY_DROP}  -i chr${CHR}.ent -q ${GERMLINE_TRIM_SEQ_FILE}.qual -m 2 -a ${ALLELE} -o xxx_report.out 
## modify the original germline output by getting rid of the bad germline
    ${PRT_CHMAP_ID} -i chr${CHR}.ent -o out2
    cut -f1 xxx_report.out >bad_germline.lst
    fgrep -f bad_germline.lst -v out2 >out2.good
    mv out2.good out2
    if test -s out2
    then
      if test -s out2
      then
## ignore alignment at the end. May consider doing both the high-quality analysis and this way by trimming
        awk '{if($6 <95 && $7 >105) print $0}' out2 >out3 
        if test -s out3
        then
          germline_count=`wc out3|awk '{printf("%s", $1)}'`
          if test ${germline_count} -le ${germline_original_count}
          then
            germline_count=${germline_original_count}
          fi
          echo "germline_count=$germline_count"
          if test ${CHECK_GERMLINE_ALLELE_COUNT} = 0
          then
            echo "$CHR	${POS}	GERMLINE=$germline_count" >>${REPORT_FILE}
            exit 1
          fi
        fi
     fi
    fi
  fi
fi

if test -s chr${CHR}.ent
then
  rm chr${CHR}.ent
fi



cp chr${CHR}.seq chr${CHR}.seq.original
cp chr${CHR}.seq.mod chr${CHR}.seq

if test -s chr${CHR}.ent
then
  rm chr${CHR}.ent
fi

## look for high-quality mismatch in germline sample
## in some cases, low-quality ones will get the germline trimmed off
## also if there is no germline reads that are germline reads of sufficient quality
high_quality_germline_count=0
if test -s ${GERMLINE_TRIM_SEQ_FILE}
then
  awk '{printf("@%s_%s\n%s\n+\n%s\n", $1, $9, $10, $11)}' test.out.clean >${GERMLINE_SEQ_FILE}
  if test -s ${GERMLINE_TRIM_SEQ_FILE}
  then
    rm ${GERMLINE_TRIM_SEQ_FILE}
  fi

  if test -s ${GERMLINE_TRIM_SEQ_FILE}.qual
  then
    rm ${GERMLINE_TRIM_SEQ_FILE}.qual
  fi
  ${TRIM_SEQ} `pwd`/${GERMLINE_SEQ_FILE} 3  ## create ${GERMLINE_TRIM_SEQ_FILE}. may need to set it to 3. but try 1 first
  if test -s ${GERMLINE_TRIM_SEQ_FILE}
  then
    max_len=`grep -v ">" ${GERMLINE_TRIM_SEQ_FILE} |awk '{printf("%ld\n", length$(1))}' |sort +0 -1 -n -r |head -n1`
    echo "max_len=$max_len"
    if test $max_len -le 30
    then
      rm ${GERMLINE_TRIM_SEQ_FILE}
    fi
  fi

  if test -s ${GERMLINE_TRIM_SEQ_FILE}
  then
    echo ">chr${CHR}.seq" >job.lst
    grep ">" ${GERMLINE_TRIM_SEQ_FILE} |cut -f2 -d">" >>job.lst
    ${WRITE_FETCH} -d ./.seq -l ${GERMLINE_TRIM_SEQ_FILE}
    if test -s chr${CHR}.ent
    then
      rm chr${CHR}.ent
    fi
    ${TSIM} -i job.lst -f2 -P1 -zF
    if test -s chr${CHR}.sat
    then
      cp chr${CHR}.seq.original chr${CHR}.seq
      ${FASTA2ASN} -i chr${CHR}.seq -p ./
      ${ATTACH_ENT} -t chr${CHR}.ent -a chr${CHR}.sat -o chr${CHR}.ent -d 1
      ${PRT_CHMAP_ID} -i chr${CHR}.ent -o out2
      if test -s out2
      then
        awk '{if($6 <=99 && $7 >=101 && $9 >=1) print $0}' out2 >out3
        if test -s out3
        then
          high_quality_germline_count=`wc out3|awk '{printf("%s", $1)}'`
        fi
      fi
    fi
  else
    germline_count=0
    germline_original_count=0
    echo "no germline"
  fi
fi


if test ${high_quality_germline_count} -gt ${germline_count}
then
  germline_count=${high_quality_germline_count}
fi

if test ${high_quality_germline_count} -gt ${germline_original_count}
then
  germline_original_count=${high_quality_germline_count}
fi


## reset the analysis by using the modified template as the master for tumro reads
if test -s chr${CHR}.ent
then
  rm chr${CHR}.ent
fi
cp chr${CHR}.seq.mod chr${CHR}.seq

## /nfs_exports/apps/gnu-apps/NextGen/SNPdetector/LINUX/src/prt_seqloc -i seqloc.lst -d ./ -b T

echo ">chr${CHR}.seq" >job.lst
grep ">" test.fastq.trim.fa |cut -f2 -d">" >>job.lst
## /nfs_exports/apps/gnu-apps/NextGen/SNPdetector/LINUX/src/write_fetch_file -d ./.seq -l test.fastq.trim.fa
${WRITE_FETCH} -d ./.seq -l test.fastq.trim.fa

## /nfs_exports/apps/gnu-apps/NextGen/SNPdetector/LINUX/src/tsim -i job.lst -f4 -P1 -zF
## ${TSIM} -i job.lst -f4 -P1 -zF
${TSIM} -i job.lst -f2 -P1 -zF
if test ! -s chr${CHR}.sat
then
  echo "$CHR	${POS}	NoSIMAlignment" >>${REPORT_FILE}
  exit 1
fi

mv chr${CHR}.seq.original chr${CHR}.seq
## /nfs_exports/apps/gnu-apps/NextGen/SNPdetector/LINUX/src/fasta2asn -i chr${CHR}.seq -p ./
if test -s chr${CHR}.ent
then
  rm chr${CHR}.ent
fi

${FASTA2ASN} -i chr${CHR}.seq -p ./
if test ! -s chr${CHR}.ent
then
  echo "$CHR    ${POS}  NoENT" >>${REPORT_FILE}
  exit 1
fi
  
${ATTACH_ENT} -t chr${CHR}.ent -a chr${CHR}.sat -o chr${CHR}.ent -d 1
## rm chr${CHR}.sat

## include the sites at the start. Or should we set to 5bp?? -k5?
${FIND_POLY} -d ./ -i chr${CHR}.ent -o out -a chr${CHR}.ent -sF -g3 -eF -c0 -k1 -MF -n2
## echo "${FIND_POLY} -d ./ -i chr${CHR}.ent -o out -a chr${CHR}.ent -sF -g3 -eF -c0 -k0 -MF -n2"


${PRT_FEATURE_LOC} -i chr${CHR}.ent -f 65 -c 0 -o snp_var.out -d F

if test ! -s snp_var.out
then
  echo "$CHR	${POS}	NoVar" >>${REPORT_FILE}
  exit 1
fi

SNP_VAR_COUNT=`awk '{printf("%s\n", $1)}' snp_var.out |sort -u |wc |awk '{printf("%s", $1)}'`

IS_VAR_CLUSTER=0
if test $SNP_VAR_COUNT -gt 3
then
 IS_VAR_CLUSTER=1
 echo "is_var_cluster"
 IGNORE_STRAND_CHECK=0
 MAX_GERMLINE_TOLERABLE_FREQ=10
 if test $germline_count -lt $germline_original_count
 then
  germline_count=$germline_original_count
 fi
fi

if test $germline_count -ge 1
then
  IGNORE_STRAND_CHECK=0
fi


big_align_count=0
if test ! -s snp_var.out
then
  echo "$CHR	${POS}	NotFound" >>${REPORT_FILE}
else
  awk '{if($1 == 100) print $0}' snp_var.out >xxx
  if test -s xxx
  then
    ${PRT_CHMAP_ID} -i chr${CHR}.ent -o out2
    if test ! -s out2
    then
      echo "$CHR	${POS}	NoAlignment" >>${REPORT_FILE}
    else
      align_count=`awk '{if($6 <90 && $7 >110 && $9 >=1) print $0}' out2 |cut -f1 |sort -u |wc |awk '{printf("%s", $1)}'`
      big_align_count=align_count
      if test ${align_count} -lt 1
      then
         echo "$CHR	${POS}	TooFewAlignment" >>${REPORT_FILE}
         exit 1
      fi

##      align_count=`awk '{if($6 <95 && $7 >105 && $9 >=1) print $0}' out2 |cut -f1 |sort -u |wc |awk '{printf("%s", $1)}'`
      align_count=`awk '{if($6 <=99 && $7 >=101 && $9 >=1) print $0}' out2 |cut -f1 |sort -u |wc |awk '{printf("%s", $1)}'`
      if test ${align_count} -lt 2
      then
         echo "$CHR	${POS}	TooFewAlignment" >>${REPORT_FILE}
         exit 1
      else
##        awk '{if($6 <95 && $7 >105 && $9 >=1) print $0}' out2 >out3
## use a lower stringency to allow for collecting reads located at the end
        awk '{if($6 <=99 && $7 >=101 && $9 >=1) print $0}' out2 >out3
        ${FILTER_CLONE_ALIGNMENT} -i out3 -r out3_report.out -o out3_no_clone.txt
        clone_count=`awk '{printf("%s\n", $1)}' out3_report.out |cut -f2 -d"="`
        unique_count=`awk '{printf("%s\n", $2)}' out3_report.out |cut -f2 -d"="`
        unique_strand=`awk '{printf("%s\n", $3)}' out3_report.out |cut -f2 -d"="`
        echo "unique_count=${unique_count}; unique_strand=${unique_strand}"
        if test ${unique_count} -lt 1
        then
          echo "${CHR}	${POS}	TooFewUniqueAlignment"	>>${REPORT_FILE}
        else
          total_unique_count=`echo "$clone_count $unique_count" |awk '{printf("%ld", $1+$2)}'`
          if test $total_unique_count -le 1
          then
            echo "${CHR}	${POS}	TooFewUniqueAlignment"	>>${REPORT_FILE}
            exit 1
          fi

          if test $clone_count -gt 0  ## anytime when there is a clone, there is always a multiple-strand info
          then
            unique_strand=3
          fi

          if test ${unique_strand} != 3
          then
            if test ${IGNORE_STRAND_CHECK} = 0
            then
              echo "${CHR}  ${POS}	SingleStrandUniqueAlignment" >>${REPORT_FILE}
            else
              unique_strand=3 ## fake this value so that we can proceed to the next step
            fi
          fi
          if test ${unique_strand} = 3
          then
            if test -s fail.out
            then
              rm fail.out
            fi
            if test -s xxx_report.out
            then
              rm xxx_report.out
            fi
            echo "${CK_QUALITY_DROP}  -i chr${CHR}.ent -q test.fastq.trim.fa.qual -m 20 -a ${ALLELE} -o xxx_report.out"
            ${CK_QUALITY_DROP}  -i chr${CHR}.ent -q test.fastq.trim.fa.qual -m 20 -a ${ALLELE} -o xxx_report.out |grep fail >fail.out
            if test -s fail.out
            then
              echo "$CHR	${POS}	QualityDrop"	>>${REPORT_FILE}
              exit 1
            fi

## checking the allele balance of germline
            if test $germline_count -gt 0
            then
              bad_allele_count=`wc xxx_report.out |awk '{printf("%s", $1)}'`
              germline_allele_freq=`echo "$align_count $bad_allele_count $germline_count" |awk '{printf("%ld", ($3 *100)/($1 -$2))}'`
              if test ${germline_allele_freq} -gt ${MIN_GERMLINE_FREQ}
              then
                if test ${germline_allele_freq} -gt ${MAX_GERMLINE_TOLERABLE_FREQ}
                then
                  echo "$CHR	${POS}	GERMLINE=$germline_count; GERMLINE_FREQ=$germline_allele_freq" >>${REPORT_FILE}
                  exit 1
                fi
              fi
            fi
            ${CREATE_INDEX_SCRIPT} test.fastq.trim.fa test.fastq.trim.fa.index ${JZ_DIR}
##            cut -f1 out3_no_clone.txt |sort -u >input_fastq.lst
##            ${GET_FASTA_SCRIPT} test.fastq.trim.fa test.fastq.trim.fa.index input_fastq.lst >input_fasta4_blat.seq
            cp test.fastq.trim.fa input_fasta4_blat.seq
##           ${RUN_BLAT_PROGRAM} input_fasta4_blat.seq chr${CHR}_${POS}_blat.out 1
            ${RUN_BLAT_PROGRAM} input_fasta4_blat.seq chr${CHR}_${POS}_blat.out 0
            if test ! -s chr${CHR}_${POS}_blat.out
            then
               echo "${CHR}	${POS}	NoBLAT" >>${REPORT_FILE}  ## this can happen when the mismatch is located too close to the end of the seq
            else
              if test -s chr${CHR}_${POS}_blat.out.good
              then
                rm chr${CHR}_${POS}_blat.out.good
              fi
              echo "chr${CHR}_${POS}_blat.out" >chr${CHR}_${POS}_blat.input
              echo "${CHECK_BLAT_REPEAT} chr${CHR}_${POS}_blat.input"
              MIN_DISTANCE_4_RREPLICATE=10
              MIN_GOOD_READ=2
              USE_LOW_THRESHOLD=1  ## use more relaxed threshold for those that have 1-bp mismatch or match exons
              if test $big_align_count = 1
              then
                MIN_GOOD_READ=1
                USE_LOW_THRESHOLD=0
                MIN_DISTANCE_4_RREPLICATE=20
              fi

              if test $germline_count -ge 1
              then
                USE_LOW_THRESHOLD=0
                MIN_DISTANCE_4_RREPLICATE=20
              fi

              if test $IS_VAR_CLUSTER = 1
              then
                USE_LOW_THRESHOLD=0
                MIN_DISTANCE_4_RREPLICATE=20
              fi

              echo "${CHECK_BLAT_REPEAT} chr${CHR}_${POS}_blat.input ${MIN_GOOD_READ} ${MIN_DISTANCE_4_RREPLICATE} ${USE_LOW_THRESHOLD}"
              ${CHECK_BLAT_REPEAT} chr${CHR}_${POS}_blat.input ${MIN_GOOD_READ} ${MIN_DISTANCE_4_RREPLICATE} ${USE_LOW_THRESHOLD}
              rerun_blat=0
              if test ! -s chr${CHR}_${POS}_blat.out.good
              then
                if test ! -s chr${CHR}_${POS}_blat.out.bad
                then
                  rerun_blat=1
                fi
              fi
              rerun_blat=0
              if test $rerun_blat = 1
              then
                ${RUN_BLAT_PROGRAM} input_fasta4_blat.seq chr${CHR}_${POS}_blat.out 0
                ${CHECK_BLAT_REPEAT} chr${CHR}_${POS}_blat.input ${MIN_GOOD_READ} ${MIN_DISTANCE_4_RREPLICATE} ${USE_LOW_THRESHOLD}
              fi

              if test -s chr${CHR}_${POS}_blat.out.good
              then
                cut -f1 out3_no_clone.txt |sort -u >no_clone.lst
                fgrep -f no_clone.lst chr${CHR}_${POS}_blat.out.good >no_clone_blat.out
                if test -s no_clone_blat.out
                then
                  echo "$CHR	${POS}	Found SNPcount=$SNP_VAR_COUNT Germline=$germline_count" >>${REPORT_FILE}
                else
                  echo "$CHR	${POS}	BadBLATClone" >>${REPORT_FILE}
                fi
              else
                echo "$CHR	${POS}	BadBLAT Germline=$germline_count" >>${REPORT_FILE}
              fi
## currently remove all the good and bad report file. But may be useful for debugging purpose
              rm chr${CHR}_${POS}_blat.input
              if test -s chr${CHR}_${POS}_blat.out.good
              then
                rm chr${CHR}_${POS}_blat.out.good
              fi
              if test -s chr${CHR}_${POS}_blat.out.bad
              then
                rm chr${CHR}_${POS}_blat.out.bad
              fi
            fi
          fi
        fi
      fi
    fi
  else
    echo "$CHR	${POS}	NotFound" >>${REPORT_FILE}
  fi
fi


## ${TSIM} -i job.lst -f1 -P1 -zF

