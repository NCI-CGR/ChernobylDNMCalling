#!/bin/bash

for epifile in $TRIOSDIR/*/*.mie.pt.filt.epi.txt ; do
	qsub runSnpEff.sh $epifile
done