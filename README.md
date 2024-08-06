
# Discerning the Genetic Footprints of Seasonal Fluctuating Selection: A Comparison with Established Selection Forms

This repository contains the code to replicate the analysis conducted in the aforementioned article. 

The [yaml](sim_env.yml) file can be used to create a conda environment with the programs required to conduct the simulations and analysis performed in this study.

Once the repository has been cloned, the environment can be created using the following code

```ruby
conda env create -f sim_env.yml

conda activate sim_env
```

A parameter file is used to pass values to simulations. The [example file](parameters/parameter_example.txt) can be used with the explanations in quotations removed. We recommend saving the parameter files as the group number (unique identifier) to the [parameter folder](parameters). The code is structured to extract parameter values from the files in this folder.

### Before you start

The following files will need the path to this directory edited:
 [python functions](scripts/single_locus_hpc.py), 
 [run full simulations](scripts/single_locus_run.py), 
  and [run short simulations](scripts/single_locus_short.py).
 
> [!NOTE]
> All R scripts will need the path to this directory as well as the chosen unique parameter set identifiers to be edited. The R files are best run in RStudio where the users can edit the code as appropriate. Plots in the R code are will be saved in the [plots](plots) folder.

