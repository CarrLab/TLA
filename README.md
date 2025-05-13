# Tumor Landscape Analysis (TLA)
 
### Landscape ecology analysis methods for the investigation and characterization of digital histopathological data of tumor biopsies.

Landscape ecology is the study of relationships between populations of living organisms interacting among themselves and with the environment they inhabit. Such an environment is defined as a __landscape__ (spread in space and time) occupied by different species of organisms, and its description entails details on their spatial distributions, cohabitation, landscape morphology, configuration, composition and fragmentation in relation to the population dynamics, mobility and other ecological effects. Our goal is to implement these methodologies in the study of tissues, perceived as cellular ecologies in the context of tumor development, and observed by means of digital histopathological samples.

__TLA__ is a toolbox composed of a large set of spatial statistics, implementing functions from the landscape ecology package [pylandstats](https://github.com/martibosch/pylandstats), astronomical and GIS spatial statistics ([astropy](https://www.astropy.org/), [pysal](https://pysal.org/esda/index.html)), spatial stratified heterogeneity ([geodetector](https://cran.r-project.org/web/packages/geodetector/vignettes/geodetector.html)) and image processing methodologies ([scipy](https://scipy.org/), [scikit-image](https://scikit-image.org/)).

Please read the [TLA documentation](documentation/TLA_doc.md) for details on the different tools implemented in this pipeline.


## Getting started

__TLA__ runs in the Anaconda python distribution. A virtual environment containing all the required dependencies can be built using the environmental file `tlaenv.yml` included in this repository. 

To get started please follow these instructions:

* First install 
[Anaconda](https://docs.anaconda.com/anaconda/install/index.html)
* Clone (or download) the __TLA__ repository in a dedicated folder.
* Use the command line terminal to run the following instructions (you can also use the Anaconda Navigator app if you prefer to use a GUI)
* Create, and update, a virtual environment named __tlaenv__ from the YML file `tlaenv.yml` (__This is the simplest method__, additional info [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)).

```
> conda update conda
> conda env create -f tlaenv.yml
> conda activate tlaenv
> conda update --all -y

``` 

* Alternativelly, build the virtual environment step by step using these commands:

	1. First create the virtual environment:

		```
		> conda install -n base -c defaults 'conda>=25.3'
		> conda update conda
		> conda update --all
		> conda create -n tlaenv python=3.10
		> conda activate tlaenv
		```
	
	2. Install `pytorch`:

		```
		> conda install -y -c pytorch pytorch torchvision torchaudio 
		```
 		or
 
		``` 
		> pip3 install torch torchvision torchaudio 
		```

	3. Install all other basic dependencies:

		```
		> conda install -y -c conda-forge pandas matplotlib-base scipy tabulate swifter
		> conda install -y -c conda-forge statannot rasterio openblas geopandas
		> conda install -y -c conda-forge scikit-image scikit-learn xarray xmltodict
		> conda install -y -c conda-forge statsmodels seaborn pyqt 
		```

	4. Finally, install `pylandstats`:

		```
		> conda install -y -c conda-forge pylandstats=3.1.0  
		```

		or 

		```
		> pip3 install pylandstats==3.1.0
		```

	5. Update and create YML (in case its needed in the future):

		```
		> conda update --all -y
		> conda env export > tlaenv.yml
		```

__NOTE 1:__ depending on the platform, it might be recommended to install `mamba` right after doing `conda update conda`, and then use the command `mamba` instead of `conda` in the following steps. 

__NOTE 2:__ if you are running TLA in a SLURM server, you can follow more specific instructions in [these notes](documentation/TLA_slurm_notes.md). 


## Workspace Structure

The following elements must coexist in the same workspace folder:

* `src/` folder containing python source code and bash scripts
* `TLA` is a bash script wrapper that handles the operation of the different modules of the pipeline. 
* `test-set.csv` is an example of an analysis argument table for cell centroid data. This file can be used as a template. So any future analysis will need to have a file with this format existing in the same folder (future builds of TLA will have a more flexible usability). This table contains all the information required for the operation of the TLA pipeline, which includes __static__ path location of the raw data to be analyzed. 
* `test-set_reg.csv` is another example of an analysis argument table. This table contains all the information required for the operation of the TLA pipeline that includes regional-segmentation and cell centroids. 

Other important files found in this repository:

* `test_set.zip` is a zip file containing example data folder:  data (aka 'raw data'), samples table (with information for each sample) and classes table with general about the cell types in the study). Please decompress after downloading to your local workspace.
* `test_set_reg.zip` is a zip file containing example data folder:  data (aka 'raw data'), samples table (with information for each sample), regions data formated as numpy (`.npy`) files, and classes table with general about the cell types and region categories in the study). Please decompress this folder after downloading to your local workspace. 
* The folder `test_set_reg/rois` contains masks for slides in the `test-set-reg` example set. Some of these slides contain more than one biopsy cut and must be split into distinct samples before processing in TLA. These masks are labeled to distinguish regions coming from different biopsies. And the `test-set_reg.csv` examples are already annotated for convenience. The TLA repositopry includes a GUI python tool, called `tla_slide_splits.py`, that would help the user generate these masks for new data. In that case, fill down the `roi_file` entry in the samples table to point to the appropiate mask. If this field is enpty or not found, TLA will assume that the whole slide is conformed by a single biopsy.

## TLA Pipeline wrapper usage

### TLA CLI script

The `TLA` wrapper bash script handles all modules of analysis. For help about using this script use the following instruction in a command line terminal:

```
> ./TLA -h
```

### TLA sbatch script

For running TLA in a SLURM array you can use the `TLA` wrapper. This script will setup an array of nodes and run different samples simultaneously and then consolidate the cohort summaries. If you have access to a HPC cluster this is the most efficient way to run the TLA pipeline. 

Copy the workspace folder in your login node. The virtual environment `tlaenv` needs to be previously installed in your account (see [these notes](documentation/TLA_slurm_notes.md)). Then the `sbatch` functions in the `TLA` script will request cluster arrays and run parallel jobs for different sample analyses. The syntax is the same as for the CLI `TLA` wrapper by use of the `-s` flag option. 

__NOTE:__ Please read the [TLA user  guide](documentation/TLA_use.md) for more details on the expected format for data files and general usage of the pipeline.

## For PC users

This pipeline was developed for Linux systems, but all python scripts should be portable to Windows.  

Each TLA module has a separate Python scripts that can be run without using the `TLA` wrapper (details can be found in the documentation). But if you desire to use the convenience of this script it is advisable to enable `bash` in your system. Please follow instruction in:

* [Install Linux on Windows with WSL](https://docs.microsoft.com/en-us/windows/wsl/install)
* [How to Enable Bash in Windows 10](https://linuxhint.com/enable-bash-windows-10/)

A simple option is to run the TLA scripts directly from python. These actions should be platform independent. See python scripts section in [documentation](documentation/TLA_doc.md) for instructions.


## File structure:

├── README.md  
├── documentation/  
│   ├── TLA\_doc.md  
│   ├── TLA\_use.md  
│   └── TLA\_slurm\_notes.md  
├── src/  
│   ├── tla\_functions.py  
|   ├── tla\_points\_\*.py    
│   ├── tla\_points\_\*.sh   
│   └── tla\_regions\_\*.py  
│   └── tla\_regions\_\*.sh   
├── TLA\_points.sh  
├── TLA\_regions.sh  
├── tlaenv.yml  
├── test\_set/  
│   ├── raw_data/  
│   │   ├── cellPos/  
│   │   ├── images/  
│   │   ├── masks/   
│   │   └── rois/  
│   └── \*.csv  
├── test\_set\_reg/  
│   ├── data/  
│   │   ├── rasters/  
│   │   ├── results/  
│   │   └── \*.csv  