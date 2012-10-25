#!/bin/sh
##/nfs_exports/apps/compilers/jdk1.6.0_15/bin/java -cp /nfs_exports/apps/internal/java/TabToExcel/commons-cli.jar:/nfs_exports/apps/internal/java/TabToExcel/Pallas-tools-0.1.jar:/nfs_exports/apps/internal/java/TabToExcel/poi-3.6-20091214.jar org.stjude.ri.pallas.tools.TabFileToExcelViewer $*

java -cp /user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess/commons-cli.jar:/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess/Pallas-tools-0.1.jar:/user/songliu/u2/group/Qiang/Exome/scripts/snp_postprocess/snv_postprocess/poi-3.6-20091214.jar org.stjude.ri.pallas.tools.TabFileToExcelViewer $*
