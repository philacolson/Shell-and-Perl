#!/bin/sh
## low_bam file was created without using the dbSNP database; however, high_bam was created by using dbSNP
## some low quality markers can "creep through" low_bam.out due to match to dbSNP

##This program is a waste of time, apparently.  We should just annotate ALL THE SNPS in one go.
SAMPLE_NAME=$1
DEBUG_LEVEL=$2
LOG=$3

##SCRIPT_DIR=/nfs_exports/apps/internal/scripts
##FIND_SUB=/nfs_exports/apps/gnu-apps/NextGen/perlsrc/FindSub.pl
SCRIPT_DIR=/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess
FIND_SUB=$SCRIPT_DIR/FindSub.pl
  
HIGH_BAM_FILE=${SAMPLE_NAME}_bam_high_20.out
if test ! -s ${HIGH_BAM_FILE}
then
  if [ $DEBUG_LEVEL -gt 0 ]
then
echo "Fail to find the bam file with high quality new results: ${HIGH_BAM_FILE} in add_refSNP2low_bam" >> $LOG
fi
exit 1
fi

LOW_BAM_SNP_FILE=${SAMPLE_NAME}_bam_low.out
if test ! -s ${LOW_BAM_SNP_FILE}
then
  if [ $DEBUG_LEVEL -gt 0 ]
then
echo "Fail to find low-quality compute data ${LOW_BAM_SNP_FILE} in add_refSNP2low_bam" >> $LOG
fi
exit 1
fi
grep rs ${LOW_BAM_SNP_FILE} |grep SNP >rs.lst
if test -s rs.lst
then
 if [ $DEBUG_LEVEL -gt 0 ]
then
echo "Low bam already had dbSNP" >> $LOG
fi
exit 1
fi

if test -s rs_high.sort
then
  rm rs_high.sort
fi
#First grab everything that has "rs" in the bam file (real snps, as opposed to the ones from Bambino).  Then check those for "SNP". (as opposed to deletions or insertions)  For each one, print out the name and the last column ($Number of Fields).  Sort them by the first column first and then by the next column. Then put them into rs_high.sort
grep rs ${HIGH_BAM_FILE} |grep SNP |awk '{printf("%s\t%s\n", $1, $NF)}' |sort +0 -1 >rs_high.sort
if test ! -s rs_high.sort #if the file is empty
then
 if [ $DEBUG_LEVEL -gt 0 ]
then
echo "If you see this, rs_high.sort is empty in add_refSNP2low_bam, meaning no rs found in ${HIGH_BAM_FILE}" >> $LOG #There were no real SNPs in the given HIGH_BAM_FILE - none of them were in dbSNP
fi
exit 1
fi

sort +0 -1 ${LOW_BAM_SNP_FILE} >low.sort #Grab the first and second columns from LOW_BAM_SNP_FILE and plop them in low.sort

#This combines the low sort and the rs_high.sort file and puts them into LOW_BAM_SNP_FILE.temp
#Combine the first field in the low.sort file with the first field in the rs_high.sort file
#Then, for every line, print out each column, including the last column, and put it in LOW_BAM_SNP_FILE
#Dr. Wang says "This line compares the two columns and gets rid of any doubles" (assumedly like an inner join?) 
#I think we do this because now we have a list of all the dbSNP equivalencies.
join -1 1 -2 1 low.sort rs_high.sort |awk '{for(i=1; i<NF; ++i) printf("%s\t", $i); printf("%s\n", $NF)}' >${LOW_BAM_SNP_FILE}.temp

#takes the first field from LOW_BAM_SNP_FILE and puts it in low_match.lst
cut -f1 ${LOW_BAM_SNP_FILE}.temp >low_match.lst

#FIND_SUB bounces the LOW_BAM_SNP_FILE entries off the reference list, low_match.lst to see if there is a match.  If there is NOT a match, we see the SNP again in the output file.
${FIND_SUB} -i ${LOW_BAM_SNP_FILE} -c low_match.lst -t 0 -d '\t' -n 0 -o low_left.out

## replace the current LOW_BAM_SNP_FILE with the modified version
cat low_left.out ${LOW_BAM_SNP_FILE}.temp >${LOW_BAM_SNP_FILE}
rm low_left.out ${LOW_BAM_SNP_FILE}.temp rs_high.sort low_match.lst
