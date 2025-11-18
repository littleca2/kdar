#!/bin/bash

WORKDIR=/home/littleca/kdar/cleanerKDAR/Michel/MC
OUTDIR=/mnt/braxton/littleca/MichelAnalysis/RATFILES

for idx in $(seq 0 500)
do
	q rat ${WORKDIR}/mac/muminus_fullvolume.mac -o ${OUTDIR}/TV_RAT_muminus_$idx.root
	q rat ${WORKDIR}/mac/muplus_fullvolume.mac -o ${OUTDIR}/TV_RAT_muplus_$idx.root
done;


