#!/bin/bash

graph=""
if [[ "$2" == "--graph" || "$3" == "--graph"]] 
then
  graph="--graph"
fi

redo=""
if [[ "$2" == "--redo" || "$3" == "--redo"]] 
then
  redo="--redo"
fi

IFS=','

# get samples_file name and determine number of samples to be processed
read name raw_path raw_samples_table raw_classes_table data_path rest < <(sed "1d" "$1")
samples_files="$raw_path/$raw_samples_table"

ncases=$(($(wc -l < $samples_files) - 1))

echo "TLA_regions: Processing ($ncases) samples in study <$1>" 

source /opt/anaconda3/etc/profile.d/conda.sh
conda activate tlaenv

# run all samples in study
for (( I=0; I<$ncases; I++ ))
do
    python src/tla_regions_run.py $1 $I $graph $redo
done   


