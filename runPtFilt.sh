#!/bin/bash
#$ -S /bin/bash
#$ -q long.q
#$ -cwd
#$ -l mem_free=20G

. /etc/profile.d/modules.sh
module load R
module load tabix/0.2.6
module load bcftools/1.4

manifest=$1
declare -A params
while IFS=$'\t' read -r -a myArray
do
 params[${myArray[0]}]=${myArray[1]}
done < $manifest

#args from manifest file
DATADIR=${params[dataDir]}
CHILDID=${params[child]}
MOMID=${params[mom]}
DADID=${params[dad]}
OUTID=${params[outId]}
BAMDIR=${params[bamDir]}

cd $DATADIR

#create ped file from manifest params
printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$OUTID" "$CHILDID" "$DADID" "$MOMID" "0" "1" > $OUTID.ped
printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$OUTID" "$DADID" "0" "0" "1" "1" >> $OUTID.ped
printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$OUTID" "$MOMID" "0" "0" "2" "1" >> $OUTID.ped

#####################################################################################################
#run PhaseByTransmission
java -Xmx16g -jar $GATK/GenomeAnalysisTK.jar -R $REFFASTA -T PhaseByTransmission -V $OUTID.raw_variants.vcf -ped $OUTID.ped -o $OUTID.pt.vcf -mvf $OUTID.mie.pt.txt

######################################################################################################
#run R script for filtering with <INFILE> and <OUTFILE> with full paths as the args
Rscript $SRCDIR/runFilt.r $DATADIR/$OUTID.mie.pt.txt $DATADIR/$OUTID.mie.pt.filt

#annotate with repeat information
awk '{rep="none"; "tabix $REFRMSK "$1":"$2"-"$2"|cut -f11" | getline rep; print $0"\t"rep  }' $OUTID.mie.pt.filt.epi.txt > $OUTID.mie.pt.filt.epi.anno.txt

#create bed file for IGV
tail -n +2 $OUTID.mie.pt.filt.mdnm.relaxed.txt | awk '{print $1"\t"$2"\t"$2}' > $OUTID.mdnm.relaxed.bed

#######################################################################################################
#ReadBackedPhasing and parent of origin
mkdir -p RBP
#filter PT VCF to add PASS to FILTER
bcftools filter $OUTID.pt.vcf -s LowQual -e 'MIN(GQ)<20 || MIN(DP)<8' > RBP/$OUTID.gqdp.filt.pt.vcf

#Create List with DNM +/-10000 bases for -L option in RBP 
tail -n +2 $OUTID.mie.pt.filt.mdnm.relaxed.txt | awk '{ if(($2-10000)<=0) print $1":1-"($2+10000); else print $1":"($2-10000)"-"($2+10000);}' > RBP/$OUTID.mdnm.relaxed.dnm.buffer.list

#Create List with DNM sites
tail -n +2 $OUTID.mie.pt.filt.mdnm.relaxed.txt | awk '{print $1":"$2"-"$2}' > RBP/$OUTID.mdnm.relaxed.dnm.list

#Call RBP and create a VCF that includes HP tag
java -Xmx16g -jar $GATK/GenomeAnalysisTK.jar -R $REFFASTA -T ReadBackedPhasing -I $BAMDIR/$CHILDID.bam --variant RBP/$OUTID.gqdp.filt.pt.vcf -o RBP/$OUTID.pt.rbp.vcf --phaseQualityThresh 20.0 -L RBP/$OUTID.mdnm.relaxed.dnm.buffer.list

#Run parent of origin script
perl $SRCDIR/extractParentOfOrigin.pl RBP/$OUTID.pt.rbp.vcf RBP/$OUTID.mdnm.relaxed.dnm.list $CHILDID RBP/$OUTID.pt.rbp.po.out