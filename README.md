
# Discerning the Genetic Footprints of Seasonal Fluctuating Selection: A Comparison with Established Selection Forms

This repository contains the code to replicate the analysis conducted in the aforementioned article. 

## Directions

### The environment
The [yaml](sim_env.yml) file can be used to create a conda environment with the programs required to conduct the simulations and analysis performed in this study.

Once the repository has been cloned, the environment can be created using the following code

```ruby
conda env create -f sim_env.yml

conda activate sim_env
```

## The simulation parameters
A parameter file is used to pass values to simulations. The [example file](parameters/fluctuating/group_x.txt) can be used with the explanations in quotations removed. We recommend saving the parameter files as ['group_x.txt'](parameters/fluctuating/group_x.txt), where x is the group number (unique identifier), to the [parameter folder](parameters) and in a folder labelled with the selection type (i.e. [fluctuating](parameters/fluctuating)). The code is structured to extract parameter values from the files in this folder and labelled in this manner. 

## The simulations
All aspects of the simulations are encoded by functions found in [this file](scripts/single_locus_hpc.py). These functions are used in [single_locus_run.py](scripts/single_locus_run.py) to run simulations with intensive sampling (and intensive resource usage) and [single_locus_short.py](scripts/single_locus_short.py), simulations with minimal sampling that are much faster and require fewer resources. These files also encode the calculation of summary statistics and in the case of the short simulations, will generate the site and haplotype frequency spectrums.
> [!IMPORTANT]
> The following files will need the path to this directory edited:
> - [single_locus_hpc.py](scripts/single_locus_hpc.py), 
> - [single_locus_run.py](scripts/single_locus_run.py),
> - [single_locus_short.py](scripts/single_locus_short.py).
>   
> The path is encoded at the top of each script for your convenience.

Once you have fixed the paths in the above Python scripts you can run the simulations using the following code, we used a temporary directory to increase speed on a cluster, if you will not be using a temporary directory you can pass the path to the results directory in it's place.

We also recommend making a folder for each selection type (i.e. hard/soft/neutral/[fluctuating](fluctuating)/balancing) and within that, a folder for each parameter set ([group](fluctuating/group_x)) as files will be deposited [in that folder](fluctuating/group_x/results_will_be_written_here.txt). Follow the links in this paragraph for an example of the file structure.

```ruby
python single_locus_short.py ${group_x} ${replicate number} ${temporary directory} ${results directory} ${selection type}
```
or
```ruby
python single_locus_run.py ${group_x} ${replicate number} ${temporary directory} ${results directory} ${selection type}
```
## Analysis
Following the successful completion of the simulation scripts, you will have the following files:
- A recombination map (rec_map_group_x.txt)
- A burn-in tree sequence (burnin_seglift_group_x_replicate.trees)
- A tree sequence from the simulation of selection (treeseq_group_x_replicate.trees)
- An allele frequency file (al_freq_group_x_replicate.txt)
- A log file that records the parameters and the number of restarts for each simulation (slimlog_group_x_replicate.txt)
- A file containing the calculated summary statistics (sim_stat_group_x_replicate.txt)
  for short simulations, you will also have:
   - A file containing the site frequency spectrum (sfs_selection_type_replicate.txt)
   - and a file containing data to calculate the haplotype frequency spectrums (hfs_group_x_replicate.txt)
       
Further analysis is conducted in R, and we recommend you conduct it in the following order:
1. [Compile summary statistic files](scripts/compile_data.R)
2. [Compare varying strengths of fluctuating selection](scripts/fluctuating_comparisons.R)
3. Generate [SFS](scripts/sfs.R) and [HFS](hfs.R)
4. [Compare differing selection types](scripts/comparing_selection_forms.R)
5. [Linear Discriminant Analysis](scripts/LDA.R)
   

> [!IMPORTANT]
> All R scripts will need the path to this directory as well as the chosen unique parameter set identifiers to be edited. The R files are best run in RStudio where the users can edit the code as appropriate. Plots in the R code will be saved in the [plots](plots) folder.

