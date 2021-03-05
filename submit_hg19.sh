#!/bin/bash

module load ucsc_userApps/v318

for bedfile in $TRIOSDIR/*/*.mdnm.relaxed.bed ; do
	triodir=$(dirname $bedfile)
	rm -rf $triodir/hg19
	mkdir -p $triodir/hg19
	cd $triodir/hg19
	awk '{print $1":"$2"-"$3}' $bedfile > hg19_input.txt
	liftOver hg19_input.txt $REFDIR/hg38ToHg19.over.chain hg19_out.txt hg19_unmapped.txt -positions
	awk '{rep="none"; "tabix $REFDIR/b38_centromeres.sorted.merged.bed.gz "$1":"$2"-"$2"|cut -f4" | getline rep; print $0"\t"rep  }' *.bedunmapped | grep -E "chr.*none" > hg19_unmapped_centromere.anno.txt
done