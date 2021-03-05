#!/bin/bash

for manifest in $TRIOSDIR/*/manifest.txt ; do
	qsub runPtFilt.sh $manifest
done