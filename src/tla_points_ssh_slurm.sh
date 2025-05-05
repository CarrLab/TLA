#!/bin/bash

graph=""
redo=""

if [ "$#" -eq 0 ]
then
    echo "Error: study argument filename is required"
    exit 0
else
    if [ "$#" -eq 2 ]
    then
        if [ $2 = "--graph" ]
        then
          graph="--graph"
        fi
        
        if [ $2 = "--redo" ]
        then
          redo="--redo"
        fi
    fi
    if [ "$#" -eq 3 ]
    then
        if [ $2 = "--graph" ] || [ $3 = "--graph" ] 
        then
          graph="--graph"
        fi
        
        if [ $2 = "--redo" ] || [ $3 = "--redo" ] 
        then
          redo="--redo"
        fi
    fi
fi


IFS=','
mkdir -p log

# get samples_file name and determine number of samples to be processed
read name raw_path raw_samples_table raw_classes_table data_path rest < <(sed "1d" "$1")
samples_files=$data_path/$name"_samples.csv"

ncases=$(($(wc -l < $samples_files) - 1))


echo "TLA_points SSH: Processing ($ncases) samples in study <$1>" 
# run all samples in a slum array
steps=$(sbatch --array=1-$ncases --parsable --export=STUDY=$1,GRAPH=$graph,REDO=$redo src/tla_points_ssh_sbatch.sh)


