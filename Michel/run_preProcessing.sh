#!/bin/bash

INDIR="/gpfs/group/mlf/nu/data/littleca/KDARSelection/preprod"
OUTDIR="/home/mlf/littleca/kdar/Michel/output_preProcessingForMichel"

for file in `ls $INDIR`
do
	RUN=${file:11:4}
	OUTFILE="$OUTDIR/michel_pre_$RUN.root"

	if [[ ! -f $OUTFILE ]]; then
		echo "Pre-processing run ${RUN}"
		export INDIR
		export OUTFILE
		export file
		bsub -n 2 -o /home/mlf/littleca/kdar/Michel/preProcessingOut.txt < sub_run_preProcessing.sh
	fi
done;



