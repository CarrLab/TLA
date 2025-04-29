#!/bin/bash
##
## TLA_points run 
##
##  This script runs TLA_points run module in a SLURM array.
##  For faster performance, samples that are already processed are skipped
##  and job requests to the Super Computer Array are adjusted accordingly. 
##  For this, several static files are created containing information about
##  partial progress of the study. IF the setup analysis fails for a group
##  of samples, it is abvisable to request more memory or time (for this, edit 
##  the `$SBATCH` parameters in the script `tla_points_setup_sample.sh`) and
##  re-run one more time to finish incomplete  samples. 
##  Repeat this process until all samples are finished. Once all 
##  samples are processed, temporary files are automatically deleted.
##
## Usage: tla_points_run_slurm ARG1 ARG2 ARG3
##
## ARGUMENTS:
##
##      ARG1: study argument file (e.g. `test_set.csv`) 
##      ARG2: graph option (i.e. "--graph" or "") 
##      ARG3: redo option (i.e. "--redo" or "") 


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

# function to count the number of lines in a file
# (accounts for files that don't end with a newline character, from Mac OS) 
numlines () {
    local eofn=$(($(tail -n 1 $1 | wc -l)))
    local nlines=$(($(wc -l < $1) - $eofn)) 
    echo $nlines
}

# get samples_file name and determine number of samples to be processed
# if some samples are already processed, they would be listed in a list 
# in the data folder with (static) name "done_samples.csv"
read name raw_path raw_samples_table raw_classes_table data_path rest < <(sed "1d" "$1")
samples_files=$data_path/$name"_samples.csv"
done_samples="$data_path/run_done_samples.csv"
sub_samples=$data_path/"run_sub_samples.csv"


if [ ! -f $done_samples ]
then
    echo "Found NO previously processed samples..."
    ncases=$(($(numlines $samples_files)))
    studtbl=$1

    echo "TLA_points_run: Processing ($ncases) samples in study <$1>" 
    mkdir -p log
    # run all samples in a slum array
    # >> use {$2, $3 = --graph, --redo} for optional graphing or redoing calculations 	
    steps=$(sbatch --array=1-$ncases --parsable --export=STUDY=$1,GRAPH=$graph,REDO=$redo src/tla_points_run_sbatch.sh)

else
    # generate a list of samples not yet processed
    echo "Found some previously processed samples..."
    
    if [[ -z "$redo" ]] 
    then
        # continue processing only cases not yet finished
        grep -w -vFf $done_samples $samples_files > $sub_samples 
        ncases=$(($(numlines $sub_samples)))
    
        if [ $ncases -eq 0 ]
        then
            ncases=$(($(numlines $samples_files)))
            echo "TLA_points_run: All ($ncases) samples are already processed in study <$1>"     
        else
            echo "TLA_points_run: Processing remaining ($ncases) samples in study <$1>" 
    
            # create a temporary parameter file pointing to sub-sample list
            studtbl=${1/'.csv'/'_sub.csv'}    
            cp $1 $studtbl	
            sed -Ei "2 s/[^,]*/run_sub_samples.csv/3" $studtbl
    
            # run all (undone) samples in a slum array
            mkdir -p log
            # >> use {$2, $3 = --graph, --redo} for optional graphing or redoing calculations 	
            steps=$(sbatch --array=1-$ncases --parsable --export=STUDY=$studtbl,GRAPH=$graph,REDO=$redo src/tla_points_run_sbatch.sh)
    
            # remove temporary file
            rm $studtbl
        fi
        
        # remove temporary file
        rm $sub_samples
        
    else
        # redo all samples
        
        # remove temporary files
        rm $done_samples
        
        echo "Redoing processing for all the samples..."
        ncases=$(($(numlines $samples_files)))
        studtbl=$1

        echo "TLA_points_run: Processing ($ncases) samples in study <$1>" 
        mkdir -p log
        # run all samples in a slum array
        # >> use {$2, $3 = --graph, --redo} for optional graphing or redoing calculations 	
        steps=$(sbatch --array=1-$ncases --parsable --export=STUDY=$1,GRAPH=$graph,REDO=$redo src/tla_points_run_sbatch.sh)
        
    fi
fi

