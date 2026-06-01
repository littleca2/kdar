#!/bin/bash

# TODO: This currently doesn't consider if we're looking at individual subruns

RUNTYPE=3	# 0: Sum (find the fit for all the runs at once. From file michel_pair_combined.root )
			# 1: Runs (find the fit for each run, using all of their subruns. From file michel_pair_####.root)
			# 2: Subruns (find the fit for each subrun. From file michel_pair_r####_s####.root)
			# 3: By date (find the combined fit for all runs starting during a given day)

VERSIONID=DAY	# The name used to label the corrections to be saved to a new directory and added to json file. Will be used in downstream analysis.


INDIR="/home/mlf/littleca/kdar/Michel/output_updateMichelPair"
DATE_LIST="/home/mlf/littleca/kdar/Michel/flux2mev_correction/same_day_runs.csv"
if [ $RUNTYPE -eq 0 ]; then
	#timestamp=$(date +%s)

	INPUT="$INDIR/michel_pair_combined.root"
	if [ -f $INPUT ]; then
		export RUNTYPE
		export INPUT
		export VERSIONID
		bsub -n 2 -o /home/mlf/littleca/kdar/Michel/flux2mev_correction/run_log.txt < sub_run_get_data_params.sh
	else
		echo "ERROR: $INPUT does not exist."
	fi

elif [ $RUNTYPE -eq 1 ]; then
	#DEBUGV=0
	for file in `ls $INDIR`
	do
		#if [ $DEBUGV -gt 1 ]; then
		#	break
		#fi

		RUN=${file:12:4}
		
		# Check if we're only processing individual runs
		re='^[0-9999]+$'
		if ! [[ $RUN =~ $re ]] ; then 
			continue
		fi

		INPUT="$INDIR/$file"
		export RUNTYPE
		export INPUT
		export VERSIONID
		bsub -o /home/mlf/littleca/kdar/Michel/flux2mev_correction/run_log.txt < sub_run_get_data_params.sh

		#echo "$RUN start"
		jobids=$(bjobs)
		update_jobids=$(bjobs)
		while [ "$jobids" == "$update_jobids" ]; do
			sleep 2
			update_jobids=$(bjobs)
		done
		#echo "$RUN end"

		#DEBUGV=$((DEBUGV+1))
	done;

elif [ $RUNTYPE -eq 3 ]; then
	# Find which runs start on the same day
	# Only need to run the python script if dates have been changed in run_dates.txt
	# python same_day_runs.py

	while read line; do
		#if [ $DEBUGV -gt 1 ]; then
		#	break
		#fi
		
		#runs=$(echo $line | tr " " "\n")
		#com=$(echo "$line" | sed 's/\s\+/,/g')	
		firstRun=${line%% *}
		INPUT="$INDIR/michel_pair_$firstRun.root"

		ADDRUNS=$(echo "$line" | sed "s/${firstRun}\s//")
		if [ ${ADDRUNS%% *} != $firstRun ]; then
			ADDRUNS=$(echo "$ADDRUNS" | sed 's/\s\+/,/g')
		else 
			ADDRUNS=""
		fi

		export RUNTYPE
		export INPUT
		export VERSIONID
		export ADDRUNS
		bsub -o /home/mlf/littleca/kdar/Michel/flux2mev_correction/run_log.txt < sub_run_get_data_params.sh

		#echo "$RUN start"
		jobids=$(bjobs)
		update_jobids=$(bjobs)
		while [ "$jobids" == "$update_jobids" ]; do
			sleep 2
			update_jobids=$(bjobs)
		done
		#DEBUGV=$((DEBUGV+1))
	done < $DATE_LIST

	#if [ -f $INPUT ]; then
	#	export 
	#fi
fi
