#!/bin/bash

INDIR="/mnt/braxton/KDAR_JADE_v3/preprod"
OUTDIR="/home/littleca/kdar/cleanerKDAR/Michel/output_preProcessingForMichel"

for file in `ls $INDIR`
do
	RUN=${file:11:4}
	OUTFILE="$OUTDIR/michel_pre_$RUN.root"

	if [[ ! -f $OUTFILE ]]; then
		./preProcessingForMichel "$INDIR/$file" "$OUTFILE"
	fi
done;



