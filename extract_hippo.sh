#!/bin/bash

dir=$1
subj_file=$2

while read -r line
do

	subj="$line"
	echo "$subj"

	mri_extract_label "$dir"/"$subj"/mri/aseg.mgz 17 "$dir"/"$subj"/mri/lh.asegHippo.mgz
	mri_extract_label "$dir"/"$subj"/mri/aseg.mgz 53 "$dir"/"$subj"/mri/rh.asegHippo.mgz

done < "$subj_file"
