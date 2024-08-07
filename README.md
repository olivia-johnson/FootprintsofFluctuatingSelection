
# Discerning the Genetic Footprints of Seasonal Fluctuating Selection: A Comparison with Established Selection Forms

This repository contains the code to replicate the analysis conducted in the aforementioned article. 

## Getting started

### The environment
The [yaml](sim_env.yml) file can be used to create a conda environment with the programs required to conduct the simulations and analysis performed in this study.

Once the repository has been cloned, the environment can be created using the following code

```ruby
conda env create -f sim_env.yml

conda activate sim_env
```

## The simulation parameters
A parameter file is used to pass values to simulations. The [example file](parameters/parameter_example.txt) can be used with the explanations in quotations removed. We recommend saving the parameter files as "group_x.txt", where x is the group number (unique identifier), to the [parameter folder](parameters). The code is structured to extract parameter values from the files in this folder.

## The simulations
All aspects of the simulations are encoded by functions found in [this file](single_locus_hpc.py). These functions are used in [single_locus_run.py](scripts/single_locus_run.py) to run simulations with intensive sampling (and intensive resource usage) and [single_locus_short.py](scripts/single_locus_short.py), simulations with minimal sampling that are much faster and require fewer resources. 
> [!NOTE]
> The following files will need the path to this directory edited:
 [python functions](scripts/single_locus_hpc.py), 
 [run full simulations](scripts/single_locus_run.py), 
  and [run short simulations](scripts/single_locus_short.py).

Once you have fixed the paths in the above Python scripts you can run the simulations using the following code
```ruby
python single_locus_short.py ${group_x} ${replicate number} ${temporary directory} ${results directory} ${selection type}
```
or
```ruby
python single_locus_run.py ${group_x} ${replicate number} ${temporary directory} ${results directory} ${selection type}
```


> [!NOTE] 
> All R scripts will need the path to this directory as well as the chosen unique parameter set identifiers to be edited. The R files are best run in RStudio where the users can edit the code as appropriate. Plots in the R code are will be saved in the [plots](plots) folder.

