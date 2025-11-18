#!/bin/bash
#usage: build_trees.py [-h] output rat proc reco index

OUTFILE_PATH=/mnt/braxton/littleca/MichelAnalysis/tree_data
RATFILE_PATH=/mnt/braxton/littleca/MichelAnalysis/RATFILES
PROC_PATH=/mnt/braxton/littleca/MichelAnalysis/PROCFILES
RECO_PATH=/mnt/braxton/littleca/MichelAnalysis/RECOFILES

ERIC_TREE_PATH=/home/marzece/KDAR_Analysis/MichelAnalysis/actual_full_volume_analysis_tree_data/

for idx in $(seq 0 500)
do
#python build_trees.py $OUTFILE_PATH/TV_TREE_muminus_$idx.root $RATFILE_PATH/TV_RAT_muminus_$idx.root $PROC_PATH/TV_PROC_muminus_$idx.root $RECO_PATH/TV_RECO_muminus_$idx.root $idx
python build_trees.py $OUTFILE_PATH/TV_TREE_muminus_$idx.root $RATFILE_PATH/TV_RAT_muminus_$idx.root $PROC_PATH/TV_PROC_muminus_$idx.root $RECO_PATH/TV_RECO_muminus_$idx.root $idx

#python build_trees.py $OUTFILE_PATH/TV_TREE_muplus_$idx.root $RATFILE_PATH/TV_RAT_muplus_$idx.root $PROC_PATH/TV_PROC_muplus_$idx.root $RECO_PATH/TV_RECO_muplus_$idx.root $idx
done;
