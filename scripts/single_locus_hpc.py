path='~/FootprintsofFluctuatingSelection/scripts/' ## Path to repository and scripts folder
import os
import sys
import msprime
import pyslim
import numpy as np
import tskit
import pandas as pd
import time
import itertools
import allel
import statistics
import scipy
sys.path.insert(1, path)
import NCD


def recombination_map(tmpdir, group, genomeSize, recRate, winSize):
    ## input group (parameter set identifier, number of chromosomes, chromosome size, recombination rate)

    ## GENERATE CHROMOSMES ##
## create recombination map for msprime and slim to simulate unlinked section
    unlinked_region = 10*winSize    
    rec_pos=[0, genomeSize-1,genomeSize, (genomeSize+unlinked_region)]
    rates=[recRate, np.log(2), recRate]
    rate_map = msprime.RateMap(position=rec_pos,rate=rates)
    rate_map
    sequenceSize=rate_map.sequence_length
    print("Sequence length is ", str(sequenceSize))
## generate slim recombination map
        # reformat to comply with recombinate map required in slim
    slim_rates=[recRate,recRate, 0.5, recRate]
    slim_pos=[0, genomeSize-1,genomeSize, (genomeSize+unlinked_region)-1]
    rec_data = {'ends': slim_pos, 'rates': slim_rates}
    slim_rec = pd.DataFrame(rec_data)
## output slim recombiation map to text file to be read into forward slim simulation
    slim_rec.to_csv("{0}/rec_map_group_{1}.txt".format(tmpdir, group), index=False, header = True, sep = "\t")

    return rate_map, sequenceSize


def single_locus_burnin(tmpdir, group, replicate, sequenceSize, s_pop, burnin_Ne, rate_map):

    ## COALESCENT BURN IN
    start_time = time.time()
    print("Burnin Ne is ", str(burnin_Ne))

    burnin = msprime.sim_ancestry(samples=s_pop, population_size=burnin_Ne, recombination_rate=rate_map, sequence_length=sequenceSize)
    ##check burnin size = genome size
    burn_length =burnin.get_sequence_length()
    if burn_length!=sequenceSize:
        print("Burnin sequence length not equal to sequence size!")
    burnin_ts = pyslim.annotate(burnin, model_type="WF", tick=1,    stage="late")
    burnin_ts.dump("{0}/burnin_seglift_group_{1}_{2}.trees".format(tmpdir,group,replicate))
    print("Time for burnin = ", (time.time()- start_time))
    
def simulate_single_locus(tmpdir, results_dir, group, replicate, recRate, genomeSize, s_pop, w_pop, h_s, h_w, s_s, s_w, rGen, fitness_on, sum_gen, win_gen, f, sim_type, winSize):

    ## FORWARD SIMULATION
    # for when uneven seasons " -d g_s=" + str(sum_gen)+ " -d g_w=" + str(win_gen)
    start_time = time.time()
    tmpdir_call = "tmpdir='" + str(tmpdir)+ "'"
    results = "results_dir='" + str(results_dir)+ "'"
    sequnceSize=genomeSize+(10*winSize)
    cmd = 'slim -d "' +str(tmpdir_call)+ '" -d "' +str(results)+ '" -d fit='+ str(fitness_on)+" -d group=" + str(group) + " -d winSize=" +
    str(winSize) +" -d f=" + str(f) + " -d h_s=" + str(h_s)+ " -d h_w=" + str(h_w)+ " -d s_s=" + str(s_s)+ " -d s_w=" + str(s_w)+
    " -d g_s=" + str(sum_gen)+" -d g_w=" + str(win_gen)+" -d replicate=" + str(replicate) + " -d GenomeSize=" + str(int(sequnceSize)) +
    " -d n_s=" + str(int(s_pop)) + " -d n_w=" + str(int(w_pop)) + " -d mut=0.0 -d rr=" + str(recRate) +   " -d rGen="+ str(rGen) +
    path+str(sim_type)+"_single_locus.slim"
    print(cmd)
    os.system(cmd)
    # #print("Time for SLiM sim = ", (time.time()- start_time))
    print("Simulations took ", (time.time()-start_time) ,  " seconds")
    
def short_single_locus(tmpdir, results_dir, group, replicate, recRate, genomeSize, s_pop, w_pop, h_s, h_w, s_s, s_w, rGen, fitness_on, sum_gen, win_gen, f, sim_type, winSize):

    ## FORWARD SIMULATION
    # for when uneven seasons " -d g_s=" + str(sum_gen)+ " -d g_w=" + str(win_gen)
    start_time = time.time()
    tmpdir_call = "tmpdir='" + str(tmpdir)+ "'"
    results = "results_dir='" + str(results_dir)+ "'"
    sequnceSize=genomeSize+(10*winSize)
    cmd = 'slim -d "' +str(tmpdir_call)+ '" -d "' +str(results)+ '" -d fit='+ str(fitness_on)+" -d group=" + str(group) + " -d winSize=" +
    str(winSize) +" -d f=" + str(f) + " -d h_s=" + str(h_s)+ " -d h_w=" + str(h_w)+ " -d s_s=" + str(s_s)+ " -d s_w=" + str(s_w)+
    " -d g_s=" + str(sum_gen)+" -d g_w=" + str(win_gen)+" -d replicate=" + str(replicate) + " -d GenomeSize=" + str(int(sequnceSize)) +
    " -d n_s=" + str(int(s_pop)) + " -d n_w=" + str(int(w_pop)) + " -d mut=0.0 -d rr=" + str(recRate) +   " -d rGen="+ str(rGen) +
    path +str(sim_type)+"_single_locus_short.slim"
    print(cmd)
    os.system(cmd)
    # #print("Time for SLiM sim = ", (time.time()- start_time))
    print("Simulations took ", (time.time()-start_time) ,  " seconds")
        
def analyse(tmpdir, results_dir,group, replicate, mutRate, genomeSize, winSize, sum_gen, win_gen, s_pop, w_pop, sim_type, s_s, s_w):
    ## input  group(parameter identifier), simulation type, simulation run), mutation rate,
    ##  number of selected loci, number of chromosomes and number of windows

     ## INPUT DATA
             # read in treesequence (ts) generated in SLiM
    slim_ts=tskit.load("{0}/treeseq_group_{1}_{2}.trees".format(tmpdir,group,replicate)).simplify()
        # extract the length of the simulate seqeunce from slim_ts
        # check number of mutations that were introduced in slim simulation
    if (slim_ts.num_sites != 1):
        print("ERROR? " + str(slim_ts.num_sites) + "introduced mutations")
    else:
        print (str(1) + " introduced mutations")

    ## SUMMARISE INDIVIDUALS - obtain metadata for individuals 'remembered' in ts
    rows_list = []
        # run through individuals (inds) in ts
    for ind in slim_ts.individuals():
        #print(ind)
        dict1 = {}
            # individual's id in ts
        dict1.update({"id" : ind.id})

        #     # genotypes individuals contained (2 nodes, each directing to a haploid genotype)
        dict1.update({"nodes": ind.nodes})

        rows_list.append(dict1)
        # convert from dictionary to data frame (fast)
    ind_met = pd.DataFrame(rows_list)

    ind_times = np.unique(slim_ts.individual_times).astype(int)


    ## ADD NEUTRAL MUTATIONS - simulations run without neutral mutations, need to put on tree to generate summary statistics removes current muts on tree
    mut_ts = msprime.sim_mutations(slim_ts, rate=mutRate, discrete_genome=False, keep=False)


    ## CALCULATE SUMMARY STATISTICS

    rows_list3 = []

        # create windows for stats to be calculated in
        # windows for tskit statistics

    win3 = np.arange(winSize/2, mut_ts.sequence_length, winSize).astype(int)
    win3 = np.append(win3, mut_ts.sequence_length)
    win3=np.insert(win3, 0, 0)
       # windows for scikit.allel statistics

    al_win3 = []  ##windows for infinite sites
    for w in range(len(win3)-1):
          if w == len(win3)-1:
              window = [win3[w]*1e10, int(mut_ts.sequence_length)*1e10]
          else:
              window = [win3[w]*1e10, win3[w+1]*1e10-1]
          al_win3.append(window)

        # cycle throught timepoints for which data has been collected
    for t in ind_times:

            # collate nodes (geotype identifiers) of indviduals at time t
        sample_ind = pyslim.individuals_alive_at(slim_ts, t)
        sample = ind_met["nodes"].iloc[sample_ind]

            ## convert to list
        samples= list(itertools.chain(*sample))

            ## create ts of just samples at time t using list
        samp_ts = mut_ts.simplify(samples = samples)

            ## create genotype array to put into scikit-allel  shape ==(var, ind)
        samp_gm=samp_ts.genotype_matrix()

            # convert genotype matrix to haplotyoe array for haplotype statistics
        h= allel.HaplotypeArray(samp_gm)
            # allele count for scikit.allel stats
        samp_ac = h.count_alleles()
            # positions of mutations in samp_ts for scikit.allel windowed_statistic function
        mut_positions = [int(var.position*1e10+1) for var in samp_ts.variants()]
                
        #### CALCULATE balancing selection NCD
        Ncd5=[]
        Ncd4=[]
        Ncd3=[]
        var=[]
        skew=[]
        kurtosis=[]
        AF=samp_ac[:,1]/sum(samp_ac[0])
        for j in range(len(win3)-1):
            win_vals=np.where((mut_positions>=win3[j]*1e10)&(mut_positions<win3[j+1]*1e10))[0]
            WAF=np.take(AF,win_vals)
            WAF=np.delete(WAF, np.where((WAF == 0.) | (WAF ==1.)))
            var.append(statistics.variance(WAF))
            skew.append(scipy.stats.skew(WAF))
            kurtosis.append(scipy.stats.kurtosis(WAF))
            MAF=np.where(WAF >0.5, 1-WAF, WAF)
            Ncd5.append(NCD.ncd(MAF, 0.5))
            Ncd4.append(NCD.ncd(MAF, 0.4))
            Ncd3.append(NCD.ncd(MAF, 0.3))
        
            # generate haplotype statistics (H1, H12, H123, H2/H1)
        hap_stats = allel.windowed_statistic(mut_positions,h,allel.garud_h, windows = al_win3)
        
            # tajimas D using tskit and branches of ts
        tajdb =  samp_ts.Tajimas_D(sample_sets=None, windows=win3, mode="branch")
        tajds =  samp_ts.Tajimas_D(sample_sets=None, windows=win3, mode="site")
        ## no. segregatiig sites
        n_seg =samp_ts.segregating_sites(sample_sets=None, windows=win3, mode="site", span_normalise=False)
        segsites =samp_ts.segregating_sites(sample_sets=None, windows=win3, mode="site")

            # wattersons theta using scikit.allel
        tw_a= allel.windowed_watterson_theta(mut_positions, samp_ac, windows=al_win3)[0]*1e10

        # tajimas D using scikit.allel
        tajda= allel.windowed_tajima_d(mut_positions, samp_ac, windows=al_win3)[0]*1e10

            # calculate diversity (tajima's pi) using tskit
        div = samp_ts.diversity(sample_sets = None, windows = win3)  ##fix windows
        ts_tests = [div, tajdb, hap_div[0]]
        for test in ts_tests:
            if len(test)!= len(win3)-1:
                print("error in test ", test, ", number of values does not match number of windows")

        ## Collate summary statics into dataframe
            # loop over windows
        for w in range(len(win3)-1):


            dict3={}
            dict3.update({"time":pyslim.slim_time(slim_ts,t)})                        ## generation
            dict3.update({"n_win":w})                       ## identifier for window
            dict3.update({"win_start" : win3[w]})            ## window start position
            dict3.update({"win_end" : win3[w+1]-1})          ## window end position
            dict3.update({"n_seg_sites" : n_seg[w]})   ## number fo segregating sites
            dict3.update({"seg_sites" : segsites[w]})   ## number fo segregating sites
            dict3.update({"diversity": div[w]})             ## diversity (tajimas pi; calculated with tsk
            dict3.update({"tajimas_d_branch":tajdb[w]})     ## tajima's D (calculated with tskit)
            dict3.update({"tajimas_d_site":tajds[w]})     ## tajima's D (calculated with tskit)
            dict3.update({"tajimas_d_allel": tajda[w]})        
            dict3.update({"theta_w_allele": tw_a[w]})        ## watterson's theta (allele with scikit allel)
            dict3.update({"ncd_5": Ncd5[w]})        ## ncd with TF=0.5
            dict3.update({"ncd_4": Ncd4[w]})        ## ncd with TF=0.4
            dict3.update({"ncd_3": Ncd3[w]})        ## ncd with TF=0.3
            dict3.update({"variance": var[w]})        ## sample variance
            dict3.update({"kurtosis": kurtosis[w]})        ## kurtosis of distirbution af
            dict3.update({"skew": skew[w]})        ## skew of distribution of af
            dict3.update({"H1": hap_stats[0][w][0]})         ## H1
            dict3.update({"H12":hap_stats[0][w][1]})         ## H12
            dict3.update({"H123": hap_stats[0][w][2]})       ## H123
            dict3.update({"H2H1": hap_stats[0][w][3]})       ## H2/H1

            rows_list3.append(dict3)

            # convert dictionary to datafram
    ts_stats = pd.DataFrame(rows_list3)

            # write statistic df to text file
    ts_stats.to_string(buf = "{2}/sim_stat_{0}_{1}.txt".format(group,replicate,results_dir), index=False) 

def sfs(sim_type, group, run, tmpdir, mutRate, winSize):  ## generate SFS
    slim_ts=tskit.load("{0}/{1}/group_{2}/treeseq_group_{2}_{3}.trees".format(tmpdir,sim_type,group,run)).simplify()
    mut_ts = msprime.sim_mutations(slim_ts, rate=mutRate, discrete_genome=False, keep=False)
    midpoint=(slim_ts.sequence_length-1e5)/2

    rows_list = []
        # run through individuals (inds) in ts
    for ind in slim_ts.individuals():
        #print(ind)
        dict1 = {}
            # individual's id in ts
        dict1.update({"id" : ind.id})

             # genotypes individuals contained (2 nodes, each directing to a haploid genotype)
        dict1.update({"nodes": ind.nodes})

        rows_list.append(dict1)
        # convert from dictionary to data frame (fast)
    ind_met = pd.DataFrame(rows_list)

    ind_times = np.unique(slim_ts.individual_times).astype(int)
    # create windows
    start_pos=(2500000-5000)-(((2500000-(winSize/2))//winSize)*winSize)

    win3 = np.arange(start_pos, mut_ts.sequence_length, winSize).astype(int)
    win=win3[win3>=(midpoint-5000)]
    win3=win[win<=(midpoint+5000)]
    win3=np.insert(win3, 0, 0)
    win3 = np.append(win3, mut_ts.sequence_length)
    # define sampling time based on selection type (Only works on short format simulations with minimal sampling)
    if (sim_type=='neutral')==True:
        samp_times=[10] # Final summer generation at the end of the simulation
    if (sim_type=='balancing')==True:
        samp_times=[(max(ind_times)-10), 10] # When first reach equilibrium and final summer generation at the end of the simulation
    if (sim_type=='fluctuating')==True:
        samp_times=[(max(ind_times)-10), 10]# When first reach equilibrium and final summer generation at the end of the simulation
    if (sim_type=='hard')==True:
        samp_times=[max(ind_times)] # When fixation is reached
    if (sim_type=='soft')==True:
        samp_times=[max(ind_times)] # When fixation is reached

    for t in samp_times:
    
        # collate nodes (genotype identifiers) of indviduals at time t
        sample_ind = pyslim.individuals_alive_at(slim_ts, t)
        sample = ind_met["nodes"].iloc[sample_ind[:100]]

            ## convert to list
        samples= list(itertools.chain(*sample))
        slim_time=pyslim.slim_time(slim_ts, t)

            ## create ts of just samples at time t using list
        samp_ts = mut_ts.simplify(samples = samples)
            ## calculate sfs
        sfs_folded =  samp_ts.allele_frequency_spectrum(sample_sets=None, windows=win3, span_normalise=False, polarised=False)
        sfs_unfolded =  samp_ts.allele_frequency_spectrum(sample_sets=None, windows=win3, span_normalise=False, polarised=True)
            ## form table and export file
        dataset = pd.DataFrame({'Folded': sfs_folded[1, 1:-1], 'Unfolded': sfs_unfolded[1, 1:-1]})
        dataset.to_string(buf = "{2}/sfs_{0}_{1}_{3}.txt".format(sim_type,run,tmpdir,slim_time), index=False)

    if sim_type=="fluctuating":  ## SFS also calculated for window 250kb away from selected site for fluctuating selection
         start_pos=(2750000-5000)-(((2750000-(winSize/2))//winSize)*winSize)

        win3 = np.arange(start_pos, mut_ts.sequence_length, winSize).astype(int)
        win=win3[win3>=(midpoint-5000)]
        win3=win[win<=(midpoint+5000)]
        win3=np.insert(win3, 0, 0)
        win3 = np.append(win3, mut_ts.sequence_length)
        for t in samp_times:
        
            # collate nodes (genotype identifiers) of indviduals at time t
            sample_ind = pyslim.individuals_alive_at(slim_ts, t)
            sample = ind_met["nodes"].iloc[sample_ind[:100]]

                ## convert to list
            samples= list(itertools.chain(*sample))
            slim_time=pyslim.slim_time(slim_ts, t)

                ## create ts of just samples at time t using list
            samp_ts = mut_ts.simplify(samples = samples)
            sfs_folded =  samp_ts.allele_frequency_spectrum(sample_sets=None, windows=win3, span_normalise=False, polarised=False)
            sfs_unfolded =  samp_ts.allele_frequency_spectrum(sample_sets=None, windows=win3, span_normalise=False, polarised=True)
            
            dataset = pd.DataFrame({'Folded': sfs_folded[1, 1:-1], 'Unfolded': sfs_unfolded[1, 1:-1]})
            dataset.to_string(buf = "{2}/sfs_250kb_{0}_{1}_{3}.txt".format(sim_type,run,tmpdir,slim_time), index=False)


def hfs(tmpdir, results_dir,group, sim_run, mutRate, winSize):  ## Generate haplotype frequency spectrum
    ## input  group(parameter identifier), simulation type, simulation run), mutation rate,
    ##  number of selected loci, number of chromosomes and number of windows

     ## INPUT DATA
             # read in treesequence (ts) generated in SLiM
    slim_ts=tskit.load("{0}/treeseq_group_{1}_{2}.trees".format(tmpdir,group,sim_run)).simplify()
        # extract the length of the simulate seqeunce from slim_ts
        # check number of mutations that were introduced in slim simulation
    if (slim_ts.num_sites != 1):
        print("ERROR? " + str(slim_ts.num_sites) + "introduced mutations")
    else:
        print (str(1) + " introduced mutations")

    ## SUMMARISE INDIVIDUALS - obtain metadata for individuals 'remembered' in ts
    rows_list = []
        # run through individuals (inds) in ts
    for ind in slim_ts.individuals():
        #print(ind)
        dict1 = {}
            # individual's id in ts
        dict1.update({"id" : ind.id})

        #     # genotypes individuals contained (2 nodes, each directing to a haploid genotype)
        dict1.update({"nodes": ind.nodes})

        rows_list.append(dict1)
        # convert from dictionary to data frame (fast)
    ind_met = pd.DataFrame(rows_list)

    ind_times = np.unique(slim_ts.individual_times).astype(int)


    ## ADD NEUTRAL MUTATIONS - simulations run without neutral mutations, need to put on tree to generate summary statistics removes current muts on tree
    mut_ts = msprime.sim_mutations(slim_ts, rate=mutRate, discrete_genome=False, keep=False)

    ## CALCULATE HFS

    rows_list3 = []

        # create windows for stats to be calculated in
        # windows for tskit statistics

    win3 = np.arange(winSize/2, mut_ts.sequence_length, winSize).astype(int)
    win3 = np.append(win3, mut_ts.sequence_length)
    win3=np.insert(win3, 0, 0)

    al_win3 = []  ##windows for infinite sites
    for w in range(len(win3)-1):
          if w == len(win3)-1:
              window = [win3[w]*1e10, int(mut_ts.sequence_length)*1e10]
          else:
              window = [win3[w]*1e10, win3[w+1]*1e10-1]
          al_win3.append(window)

        # cycle throught timepoints for which data has been collected
    for t in ind_times:

            # collate nodes (geotype identifiers) of indviduals at time t
        sample_ind = pyslim.individuals_alive_at(slim_ts, t)
        sample = ind_met["nodes"].iloc[sample_ind]

            ## convert to list
        samples= list(itertools.chain(*sample))

            ## create ts of just samples at time t using list
        samp_ts = mut_ts.simplify(samples = samples)

        for w in range(len(win3)-1):

            inter=np.array([win3[w],win3[w+1]-1])
            inter=np.reshape(inter, (1,2))

            win_ts=samp_ts.keep_intervals(inter, simplify=False)

            ## create genotype array to put into scikit-allel  shape ==(var, ind)
            samp_gm=win_ts.genotype_matrix()

            # convert genotype matrix to haplotyoe array for haplotype statistics
            h= allel.HaplotypeArray(samp_gm)

            haps=np.swapaxes(h,0,1) ##swap axis so each row is a hap

            strings=[np.array_str(a) for a in haps] ## make each hap into a string

            unique_haps=np.unique(strings)  ## isolate unique haplotypes

            for i in range(len(unique_haps)):
                hap=unique_haps[i]
                count = sum(np.char.count(strings, hap))
                dict3={}
                dict3.update({"time":pyslim.slim_time(slim_ts,t)})                        ## generation
                dict3.update({"n_win":w})                       ## identifier for window
                dict3.update({"win_start" : win3[w]})            ## window start position
                dict3.update({"win_end" : win3[w+1]-1})          ## window end position
                dict3.update({"hap_id" : i})          ## identifier for unique haplotype
                dict3.update({"count" : count})          ## count of occurances of that haplotype
                rows_list3.append(dict3)

            # convert dictionary to datafram
    ts_stats = pd.DataFrame(rows_list3)

            # write  df to text file
    ts_stats.to_string(buf = "{2}/hfs_{0}_{1}.txt".format(group,sim_run,results_dir), index=False)

            
        





            
        






