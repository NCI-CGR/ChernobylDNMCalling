#!/bin/bash
#$ -S /bin/bash
#$ -q all.q
#$ -cwd
#$ -l mem_free=4G

. /etc/profile.d/modules.sh;
module load snpEff/4.1
module load bcftools/1.9

epifile=$1

triodir=$(dirname $epifile)
rm -rf $triodir/snpEff
mkdir -p $triodir/snpEff
cd $triodir/snpEff
awk 'FNR > 1 {print $1":"$2"-"$2}' $epifile > epi.list
java -Xmx16g -jar $GATK/GenomeAnalysisTK.jar -R $REFFASTA -T SelectVariants --variant $triodir/*.pt.vcf -L epi.list -o epi.vcf

#ensembl
snpEff -dataDir $SNPEFFDATA GRCh38.79 epi.vcf | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\t%INFO/LOF\t%INFO/NMD\n' - > epi.snpeff.ensembl.txt
#refseq
snpEff -dataDir $SNPEFFDATA GRCh38.p2.RefSeq epi.vcf | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\t%INFO/LOF\t%INFO/NMD\n' - > epi.snpeff.refseq.txt