#!/bin/bash
PROCPATH=/mnt/braxton/littleca/MichelAnalysis/PROCFILES
RECOPATH=/mnt/braxton/littleca/MichelAnalysis/RECOFILES

for idx in $(seq 2 500)
do 
q jade_reco -m Q  -i $PROCPATH/TV_PROC_muminus_$idx.root -o $RECOPATH/TV_RECO_muminus_$idx.root
q jade_reco -m Q  -i $PROCPATH/TV_PROC_muplus_$idx.root -o $RECOPATH/TV_RECO_muplus_$idx.root
done;

