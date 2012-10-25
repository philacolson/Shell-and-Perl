#!/bin/sh

BIN_DIR=$3

SORT_STRING=${BIN_DIR}/sort_string
CREATE_INDEX=${BIN_DIR}/CreateIndex4Fasta.prl

if test -n "$1" && test -n "$2"
then
echo Creating unsorted index
${CREATE_INDEX} $1 temp.index
echo Sorting
## sort -t " " +0 -1 temp.index > $2
${SORT_STRING} -i temp.index -o $2
echo Deleting unsorted index
rm temp.index
else
echo usage CreateIndex4Fasta.scr InFileName OutIndexFileName 
fi
