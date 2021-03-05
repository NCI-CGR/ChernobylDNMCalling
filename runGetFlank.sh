#!/bin/bash
#$ -S /bin/bash
#$ -q long.q
#$ -cwd
#$ -l mem_free=10G

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

cd $DATADIR
mkdir -p FLANKING
filename="$OUTID.flank.txt"
rm -f FLANKING/$OUTID.flank.out
while IFS=$'\t' read -r -a myArray
do
	chr="${myArray[0]}"
	pos="${myArray[1]}"

	#upstream
	interval=0
	if [ "$pos" -le 900000 ]; then
		interval="$chr:1-$pos"
	else
		interval="$chr:$((pos - 900000))-$pos"
	fi
	
	#Extract VCF with +/-200000 DNM
	java -Xmx16g -jar $GATK/GenomeAnalysisTK.jar -R $REFFASTA -T SelectVariants -V $OUTID.raw_variants.vcf -o FLANKING/$OUTID.flank.vcf -L $interval
	bcftools view -Ou -s $CHILDID,$MOMID,$DADID FLANKING/$OUTID.flank.vcf | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | grep -E "0\/1.*0\/0.*0\/1|0\/1.*0\/1.*0\/0" | awk -v position="$pos" -v outid="$OUTID" '{ print outid,"\t",position,"\t",$0,"\t",(position-$2)}' | sort -k 10,10 -V | head -2 >> FLANKING/$OUTID.flank.out

	#downstream
        interval="$chr:$pos-$((pos + 900000))"

        #Extract VCF with +/-200000 DNM
        java -Xmx16g -jar $GATK/GenomeAnalysisTK.jar -R $REFFASTA -T SelectVariants -V $OUTID.raw_variants.vcf -o FLANKING/$OUTID.flank.vcf -L $interval
        bcftools view -Ou -s $CHILDID,$MOMID,$DADID FLANKING/$OUTID.flank.vcf | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | grep -E "0\/1.*0\/0.*0\/1|0\/1.*0\/1.*0\/0" | awk -v position="$pos" -v outid="$OUTID" '{ print outid,"\t",position,"\t",$0,"\t",(position-$2)}' | sort -k 10,10 -V | head -2 >> FLANKING/$OUTID.flank.out
	
done < $filename