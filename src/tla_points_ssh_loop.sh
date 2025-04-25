#!/bin/bash

IFS=','

# get samples_file name and determine number of samples to be processed
read name raw_path raw_samples_table raw_classes_table data_path rest < <(sed "1d" "$1")
samples_files=$data_path/$name"_samples.csv"

ncases=$(($(wc -l < $samples_files) - 1))

echo "TLA_points SSH: Processing ($ncases) samples in study <$1>" 

source /opt/anaconda3/etc/profile.d/conda.sh
conda activate tlaenv

redo=''
if (( $# < 3 )); then
    redo=$2
else
    redo=$3
fi

# run all samples in study
for (( I=0; I<$ncases; I++ ))
    do
	python src/tla_points_ssh.py $1 $I $redo
    done









