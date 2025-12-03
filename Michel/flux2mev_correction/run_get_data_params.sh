#!/bin/bash

# TODO: This currently doesn't consider if we're looking at individual subruns

RUNTYPE=	# 0: Sum (find the fit for all the runs at once. From file michel_pair_s###.root )
			# 1: Runs (find the fit for each run, using all of their subruns. From file michel_pair_####.root)
			# 2: Subruns (find the fit for each subrun. From file michel_pair_r####_s####.root)

VERSIONID=	# The key timestamp for this set of data, usually same as the Sum value found prior

ERICDATA=0	# Are we using Eric's data? This will be deleted later, probably

if [ $ERICDATA -eq 0 ]; then
	INDIR="/home/littleca/kdar/Michel/output_updateMichelPair"
	#OUTDIR="/home/littleca/kdar/cleanerKDAR/Michel/flux2mev_correction/output_fluxCorr"

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

else	#Using Eric's data
	INDIR="/home/littleca/kdar/Michel/old_eric_updateMichelPair"
	OUTDIR="/home/littleca/kdar/Michel/flux2mev_correction/output_fluxCorr"
	if [ $RUNTYPE -eq 0 ]; then
		#timestamp=$(date +%s)

		INPUT="$INDIR/michel_analy_combined.root"
		if [ -f $INPUT ]; then
			python get_data_params.py "$INPUT" "$RUNTYPE" "$VERSIONID"
		else
			echo "ERROR: $INPUT does not exist."
		fi

	elif [ $RUNTYPE -eq 1 ]; then

		for file in `ls $INDIR`
		do
			RUN=${file:14:4}
			
			re='^[0-9999]+$'
			if ! [[ $RUN =~ $re ]]; then
				continue
			fi

			echo "python get_data_params.py $INDIR/$file $RUNTYPE $VERSIONID"
			python get_data_params.py "$INDIR/$file" "$RUNTYPE" "$VERSIONID"
		done;
	fi
fi
