#!/bin/bash

INDIR="/home/mlf/littleca/kdar/Michel/output_preProcessingForMichel"
OUTDIR="/home/mlf/littleca/kdar/Michel/output_updateMichelPair"

for file in `ls $INDIR`
do
	RUN=${file:11:4}
	echo "$RUN"
	OUTFILE="$OUTDIR/michel_pair_$RUN.root"

	if [[ ! -f $OUTFILE ]]; then
		./update_michel_pair "$INDIR/$file" "$OUTFILE"
	fi

done;

rm $OUTDIR"/michel_pair_combined.root"
hadd $OUTDIR"/michel_pair_combined.root" $OUTDIR"/*.root"


