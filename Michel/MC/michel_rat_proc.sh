#!/bin/bash
RAT_FILE_LOCATION=/mnt/braxton/littleca/MichelAnalysis
JADE_GAIN_FILE=$JADEROOT/dat/highGainFitsLS.root
JADE_GAIN_RATIO_FILE=$JADEROOT/dat/gainRatiosLS.root
JADE_ARGS_STR=-g\ $JADE_GAIN_FILE\ -r\ $JADE_GAIN_RATIO_FILE

for idx in $(seq 0 500)
do 
q jade_process_data $JADE_ARGS_STR -i $RAT_FILE_LOCATION/RATFILES/TV_RAT_muminus_$idx.root -o $RAT_FILE_LOCATION/PROCFILES/TV_PROC_muminus_$idx.root

q jade_process_data $JADE_ARGS_STR -i $RAT_FILE_LOCATION/RATFILES/TV_RAT_muplus_$idx.root -o $RAT_FILE_LOCATION/PROCFILES/TV_PROC_muplus_$idx.root
done;

