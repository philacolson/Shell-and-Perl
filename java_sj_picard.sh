#!/bin/sh
# Calls java 1.6 with the classpath set to include the compbio code,
# commons-cli, and picard; also increases MaxPermSize
#
# $1 = bare (without package) runnable class name
# $2... = all others are sent as arguments to the java program

VERS=1.31
CLASS=$1
shift
#/nfs_exports/apps/compilers/jdk1.6.0_15/bin/java -Xmx2g -XX:MaxPermSize=1280M -cp /nfs_exports/apps/gnu-apps/NextGen/mrbin:/nfs_exports/apps/gnu-apps/NextGen/mrbin/sam-$VERS.jar:/nfs_exports/apps/gnu-apps/NextGen/mrbin/picard-$VERS.jar:/nfs_exports/apps/gnu-apps/NextGen/mrbin/commons-cli-1.2.jar:/nfs_exports/apps/gnu-apps/NextGen/mrbin/compiled_jars/TweakSam.jar org.stjude.compbio.sam.cmd.$CLASS $*

java -Xmx2g -XX:MaxPermSize=1280M -cp user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess:user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess/sam-$VERS.jar:user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess/picard-$VERS.jar:user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess/commons-cli-1.2.jar:user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess/compiled_jars/TweakSam.jar org.stjude.compbio.sam.cmd.$CLASS $*

