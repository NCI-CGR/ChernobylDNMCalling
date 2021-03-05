#!/bin/bash

. /etc/profile.d/modules.sh;
module load zlib/1.2.11
module load bedtools/2.27.1

for epifile in $TRIOSDOR/*/snpEff/epi.snpeff.ensembl.txt ; do
	triodir=$(dirname $epifile)
	cd $triodir
	awk 'FNR > 1 {print $1"\t"$2-2"\t"$2+1"\t"$1":"$2":"$3"/"$4}' $epifile > epi.bed
	bedtools getfasta -fi $REFFASTA -bed epi.bed -name -tab > epi.trinuc.txt
	bedtools intersect -a epi.bed -b $REFCFS -loj | cut -f1,2,3,4,5,6,7,8 > epi.humcfs.txt
done