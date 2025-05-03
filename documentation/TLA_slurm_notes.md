
# Steps to create `tlaenv` in a HPC cluster SLURM array

This are the steps to create a virtual environment and get all the required packages for TLA to run in a SLURM array.

* First initiate an interactive node in the HPC server and load an anaconda, miniconda or mamba module before creating a virtual environment, and then activate it. __It is advisable to contact the IT administrator of the HPC system to get specific instructions for this step.__ Here is a general set of instructions, but they might be different depending of the configuration of the system. __The user should make the appropiate changes according to the administrators' advice__. 
* If `mamba` is installed, which is a package manager more efficient than `conda`, load it and create/activate a new environment:

```
> interactive
> module load mamba/latest
> mamba create --name tlaenv python=3.10
> source activate tlaenv
```

* If the server doesn't have mamba instaleld (ie. the `module load mamba/latest` command fails) you can just replace "mamba" with "conda" everywhere in the rest of the instructions, __or__ you can install mamba by doing this instead:

```
> interactive
> module load conda/latest
> conda create --name tlaenv -y -c conda-forge python=3.10 mamba
> source activate tlaenv
```

* Once `tlaenv` is activated, install pytorch and the rest of the dependencies for TLA:

```
> mamba install -y -c torch torchvision torchaudio
> mamba install -y -c conda-forge pandas matplotlib-base 
> mamba install -y -c conda-forge scipy tabulate 
> mamba install -y -c conda-forge swifter
> mamba install -y -c conda-forge statannot rasterio 
> mamba install -y -c conda-forge openblas geopandas
> mamba install -y -c conda-forge scikit-image scikit-learn 
> mamba install -y -c conda-forge xarray xmltodict
> mamba install -y -c conda-forge statsmodels seaborn pyqt
> mamba install -y -c conda-forge pylandstats=3.1.0
> mamba update --all
```

* Finally, you can save the environment YML file in case it's needed again in the future

```
> mamba env export > tlaenv.yml
```

* The conda environment is now ready to run TLA!

