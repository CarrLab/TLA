#!/bin/bash

PTH=$( dirname "$0")
SRC=$PTH"/src/"

Help(){
   # Display Help
   echo "#####################################################################"
   echo "# Tumor Landscape Analysis (TLA) - Regions                          #"
   echo "#                                                                   #"
   echo "#  TLA_regions: region-level analysis of tissue compartment data    #"
   echo "#  The expected input for each sample is a raster array with        #" 
   echo "#  discrete labels representing category values (class) for each    #"
   echo "#  tissue compartment that has been segmented and classified.       #"
   echo "#                                                                   #"
   echo "#  For more info please visit: https://github.com/cisluis/TLA       #"
   echo "#                                                                   #"
   echo "# Syntax: TLA_regions [options] action study                        #"
   echo "#                                                                   #"
   echo "#  options:                                                         #"
   echo "#   -h   print this help                                            #"
   echo "#   -v   print TLA_regions modules versions                         #"
   echo "#   -l   print license statement                                    #"
   echo "#                                                                   #"
   echo "#  switches:                                                        #"
   echo "#   -s   'sbatch' mode is TRUE and TLA runs in a SLURM array,       #"
   echo "#        otherwise it runs is serial mode (for loop).               #"
   echo "#   -g   if used all sample graphs are plotted                      #"
   echo "#   -r   if used all data is re-analyzed from scratch.              #"
   echo "#        CAREFUL: any pre-existing results WILL BE OVERWRITEN       #"
   echo "#                                                                   #"
   echo "#  arguments:                                                       #"
   echo "#   action: 'run' for TLA analysis                                  #"
   echo "#           'sum' for summarizing TLA results                       #"
   echo "#           'clean' clean study folder with all analysis results    #"
   echo "#   study:  name of study set to process, corresponding with a file #"
   echo "#           <study_name.csv> with arguments for the analysis        #"
   echo "#                                                                   #"
   echo "#####################################################################"
   echo "# Copyright (C) 2025, Luis H. Cisneros <cisluis@gmail.com>          #"
   echo "#                     Arizona Cancer Evolution Center               #"
   echo "#                     Arizona State University. Tempe Arizona       #"
   echo "#                     Mayo Clinic, Rochester, Minnesota             #"
   echo "#                                                                   #"
   echo "# This software is free; you can redistribute it and/or modify it   #"
   echo "# it under the terms of the GNU General Public License as published #"
   echo "# in the URL: https://www.gnu.org/licenses/gpl-3.0.en.html          #"
   echo "#                                                                   #"
   echo "#####################################################################"
   echo
   
   exit 2
}

version_number () {
    f=$SRC$1
    v=$(grep '^__version__' $f)
    s="==> $1: version = ${v#*=}"
    echo $s
}

Versions(){
   echo "#####################################################################"
   echo " Tumor Landscape Analysis (TLA)  - Module versions:                  "
   
   echo $(version_number "tla_functions.py")
   echo $(version_number "tla_regions_run.py")
   echo $(version_number "tla_regions_sum.py")
 
   echo "#####################################################################"
   echo
   
   exit 2
}

slurm=FALSE
graph=''
redo=''

while getopts ":hvlsgr" option; do
    case "$option" in

    h)
        Help
        ;;
    v)
        Versions
        ;;
    l)
        less LICENSE
	    exit;;
    s)
        slurm=TRUE
        ;;
    g)
        graph="--graph"
        ;;
    r)
        redo="--redo"
        ;;
    *)
        Help
        ;;
    esac
done

shift $((OPTIND-1))

##############################################################################
# Program call                                                               #
##############################################################################

ARG1=${1:-''}
ARG2=${2:-'test_set.csv'}

case $ARG1 in
    run) # run TLA module
    	case $slurm in
      	   TRUE) # slurm run
             	source $SRC"tla_regions_run_slurm.sh" $ARG2 $graph $redo
             	exit;;
           FALSE) # serial run
           	source $SRC"tla_regions_run_loop.sh" $ARG2 $graph $redo
              	exit;;
    	esac
        ;;
        
    sum) # sum module
        case $slurm in
      	   TRUE) # slurm run
              	 source $SRC"tla_regions_sum_slurm.sh" $ARG2 
                 exit;;
       	   FALSE) # serial run
               	 source $SRC"tla_regions_sum_loop.sh" $ARG2
                 exit;;
    	esac
    	;;
    
   clean) # delete results
       echo "Are you sure you wish to clean cache for studies in <$ARG2>?"
       select yn in "Yes" "No"; do
           case $yn in
                Yes ) 
                    {
                       read    
                       while IFS=, read -r nam rph sam cls pth res ; do
                           echo "... cleaning data for study <$nam>"
                           datadir=${pth}
                           echo $datadir
                           rm -rf $datadir
                       done
                     } < $ARG2; 
                     exit;; 
                No ) 
                    echo "... nothing was done!"
                    exit;;
             esac
       done
       exit;;
	
   *) # unknown action
      echo -n "... ERROR: Action <<$ARG1>> is not a recognized action!"
      echo
      Help
      exit;;

esac


