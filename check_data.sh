#!/bin/bash

FILEPATH="/mnt/braxton/KDAR_JADE_v3/preprod"

for file in $(ls $FILEPATH)
do
	FILESIZE_FLOAT=$(awk '{print substr($1, 1, length($1)-1)}' <<< $(du -h $FILEPATH"/"$file))
	FILESIZE=${FILESIZE_FLOAT/.*}
	FILESUFF=$(awk '{print substr($1, length($1), length($1))}' <<< $(du -h $FILEPATH"/"$file))

#	echo "$FILESIZE $FILESUFF		$file"

	if [[ $FILESIZE -lt 5 ]] && [[ $FILESUFF == "K" ]]
	then
		echo "$FILESIZE $FILESUFF		$file"
	fi
done;
