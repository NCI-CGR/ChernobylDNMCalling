#!/bin/bash
#$ -S /bin/bash
#$ -q xlong.q
#$ -cwd
#$ -l mem_free=25G

. /etc/profile.d/modules.sh; module load jdk/1.8.0_111

CHILDBAM=$1
MOMBAM=$2
DADBAM=$3
OUTDIR=$4

cd $OUTDIR

java -Xmx16G -Xms4G -jar $GATK/Queue.jar -jobRunner Drmaa -jobNative '-clear -q xlong.q -pe by_node 1' -jobReport $OUTDIR/variantCaller.jobReport.txt -S HaplotypeCaller.scala -R $REFFASTA -O $OUTDIR -mem 4 -nct 4 -nsc 333 -I $CHILDBAM -I $DADBAM -I $MOMBAM -D $REFDBSNP -XL $REFEXCLUDELIST -run