path="~/FootprintsofFluctuatingSelection"
import os
import sys
import msprime
import pyslim
import numpy as np
import yaml
import tskit
import pandas as pd
import time
import allel
import scipy
sys.path.insert(1, '{0}scripts/'.format(path))
import single_locus_hpc
import NCD

params=sys.argv[1]
replicate = sys.argv[2]
tmpdir =str(sys.argv[3])
results_dir=str(sys.argv[4])

####  READ IN PARAMETERS
    # load in parameter file
with open('{0}parameters/{1}.txt'.format(path, params), 'r') as f:
    parameters = yaml.load(f, Loader=yaml.FullLoader)

    #set parameters from file
genomeSize = int(parameters["genomeSize"])
recRate = parameters["recRate"]
mutRate = parameters["mutRate"]
s_pop = int(parameters["s_pop"])
w_pop = int(parameters["w_pop"])
h_s = parameters["h_s"]
h_w = parameters["h_w"]
s_s = parameters["s_s"]
s_w = parameters["s_w"]
rGen=int(parameters["rGen"])
fitness_on = parameters["fitness_on"]
sum_gen = int(parameters["sum_gen"])
win_gen = int(parameters["win_gen"])
group=parameters["group"]
winSize = parameters["winSize"]
freq=int(parameters["f"])
sim_type=str(parameters["sel_type"])

start_time = time.time()

## Calculate Burnin Ne
burnin_Ne = round((sum_gen+win_gen)/(((1/s_pop)*sum_gen)+((1/w_pop)*win_gen)))

####  Simulate burnin
if os.path.exists('{0}/burnin_seglift_group_{1}_{2}.trees'.format(results_dir, group,replicate))==False:
    rate_map, sequenceSize = single_locus_hpc.recombination_map(tmpdir, group, genomeSize, recRate, winSize)
    single_locus_hpc.single_locus_burnin(tmpdir, group, replicate, sequenceSize, s_pop, burnin_Ne, recRate)

####  Simulate selection
single_locus_hpc.short_single_locus(tmpdir, results_dir, group, replicate, recRate, genomeSize, s_pop, w_pop, h_s, h_w, s_s, s_w, rGen, fitness_on, sum_gen, win_gen, freq, sim_type, winSize)

####  Run analysis
single_locus_hpc.analyse(tmpdir, results_dir, group, replicate, mutRate, genomeSize, winSize, sum_gen, win_gen, s_pop, w_pop, sim_type, s_s, s_w)

#### Calculate Site Frequency Spectrum
single_locus_hpc.sfs(sim_type, group, replicate, tmpdir, mutRate, winSize)

#### Calculate Haplotype Frequency Spectrum
single_locus_hpc.hfs(tmpdir, results_dir, group, replicate, mutRate, winSize)


print("Time = ", (time.time() - start_time))

