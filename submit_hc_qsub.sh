#!/bin/bash

for manifest in $TRIOSDIR/*/manifest.txt ; do
    declare -A params
    while IFS=$'\t' read -r -a myArray
    do
      params[${myArray[0]}]=${myArray[1]}
    done < $manifest
    
    CHILDID=${params[child]}
    MOMID=${params[mom]}
    DADID=${params[dad]}
    OUTID=${params[outId]}
    BAMDIR=${params[bamDir]}
    DATADIR=${params[dataDir]}
   
    cd $DATADIR 
    qsub runHc_b38_qsub.sh $BAMDIR/${CHILDID}.bam $BAMDIR/${MOMID}.bam $BAMDIR/${DADID}.bam $DATADIR 
done