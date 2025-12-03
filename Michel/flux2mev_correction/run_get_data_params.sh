#!/bin/bash

# TODO: This currently doesn't consider if we're looking at individual subruns

RUNTYPE=	# 0: Sum (find the fit for all the runs at once. From file michel_pair_s###.root )
			# 1: Runs (find the fit for each run, using all of their subruns. From file michel_pair_####.root)
			# 2: Subruns (find the fit for each subrun. From file michel_pair_r####_s####.root)

VERSIONID=	# The name used to label the corrections to be saved to a new directory and added to json file. Will be used in downstream analysis.


INDIR="/home/littleca/kdar/Michel/output_updateMichelPair"

if [ $RUNTYPE -eq 0 ]; then
	#timestamp=$(date +%s)

	INPUT="$INDIR/michel_pair_combined.root"
	if [ -f $INPUT ]; then
		python get_data_params.py "$INPUT" "$RUNTYPE" "$VERSIONID"
	else
		echo "ERROR: $INPUT does not exist."
	fi

elif [ $RUNTYPE -eq 1 ]; then
	for file in `ls $INDIR`
	do
		RUN=${file:12:4}
		
		# Check if we're only processing individual runs
		re='^[0-9999]+$'
		if ! [[ $RUN =~ $re ]] ; then 
			continue
		fi

		echo "python get_data_params.py $INDIR/$file $RUNTYPE $VERSIONID"
		python get_data_params.py "$INDIR/$file" "$RUNTYPE" "$VERSIONID"
	done;
fi
