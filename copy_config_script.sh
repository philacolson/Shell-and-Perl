#!/bin/sh

## this program copies the configuration files required for perform the SNP annotation analysis
SOURCE_DIR=$1           ## use /nfs_exports/apps/gnu-apps/NextGen/SeattleSNPs for default
DESTINATION_DIR=$2	## use /nfs_exports/apps/gnu-apps/NextGen/JJLungCancer

if test ! -d ${SOURCE_DIR}
then
  echo "Fail to find the source directory ${SOURCE_DIR}"
  exit 1
fi

if test ! -d ${DESTINATION_DIR}
then
  echo "Fail to find the source directory ${DESTINATION_DIR}"
  exit 1
fi

## file for running blast search
if test ! -s ${SOURCE_DIR}/human_db.lst
then
  echo "Fail to find file ${SOURCE_DIR}/human_db.lst"
  exit 1
fi
cp ${SOURCE_DIR}/human_db.lst ${DESTINATION_DIR}/.

## file for repeat mask
if test ! -s ${SOURCE_DIR}/humrep.fsa
then
  echo "Fail to find file ${SOURCE_DIR}/humrep.fsa"
  exit 1
fi
cp ${SOURCE_DIR}/humrep.fsa ${DESTINATION_DIR}/.


## configuration file for blast search
if test ! -d ${SOURCE_DIR}/config
then
  echo "Fail to find directory ${SOURCE_DIR}/config"
  exit 1
fi
cp -r ${SOURCE_DIR}/config ${DESTINATION_DIR}/.
