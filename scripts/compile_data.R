
## Code to compile output from python scripts into an R data file that can be used for further anaylsis ##


## Load packages and install any missing from the user's R library
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table)

## Set path to directory
path='/Path/to/FootprintsofFluctuatingSelection/'
## Set selection type 
sim_type="hard" 
## Set working directory
setwd(paste0(path, sim_type, "/"))

## Vector of unique parameter IDs
groups=c(1:6)  

## Create Allele frequency and summary staistic data tables
freq_data = NULL
sum_stats = NULL

## A loop that will compile data by replicateing through unique paramater identifiers, 
## extract the parameters files, allele frequency files and summary statsics files 
## and compile them into two data.tables (freq_data and sum_stats)

for (g in groups){
  print(g)
  
  f_list  <- list.files(path =paste0(path, sim_type,"/group_", g, "/"),pattern ="al_freq_")

  s_list  <- list.files(path =paste0(path, sim_type,"/group_", g, "/"),pattern ="sim_stat_")

  parameters <- fread(file=paste0(path, "parameters/", sim_type,"/group_",g, ".txt"), sep = ":")
  setkey(parameters, V1)
  
  fiton=parameters["fitness_on", V2]
  if (fiton==0){
    h_s = 0
    h_w = 0
    s_s = 0
    s_w = 0
  } else{
    h_s = parameters["h_s", V2]
    h_w = parameters["h_w", V2]
    s_s = parameters["s_s", V2]
    s_w = parameters["s_w", V2]
  }
  sum_gen =parameters["sum_gen", V2]
  win_gen =parameters["win_gen", V2] 
  if (sum_gen==win_gen){
    gen_s = "EG"
  } else{
    gen_s = "UG"
  }
  sum_pop =parameters["s_pop", V2]
  win_pop =parameters["w_pop", V2]
  if (sum_pop==win_pop){
    if (sum_pop==10000){
      pop_s = "CP10k"
    }else{
    pop_s = "CP1k"}
  } else{
    pop_s = "FP"
  }
  mr = parameters["mutRate", V2]
  rr = parameters["recRate", V2]
  nWin = parameters["winpChrom", V2]
  genomeSize = parameters["genomeSize", V2]
  
  if (s_s==s_w){ g_label = ifelse(fiton==0, paste( pop_s, gen_s, "No Fitness", sep="_"), paste(pop_s, gen_s, "h", h_s,"s", s_s, sep="_"))}else{
    eq_freq=s_s/(s_s+s_w)
    g_label=paste(pop_s, gen_s, "h", h_s,"s", round(s_s, digits=3),"eq", eq_freq, sep="_")}
  # collate al_freq files
  for (i in 1:length(f_list)){
     filename = f_list [i]
    al_freq  = fread(file = paste0(path, sim_type,"/group_",g, "/",filename))
    al_freq[1, Gen:=1]
    al_freq[,replicate:=i]
    al_freq[,label:=g_label]
    al_freq[,group:=g]
    al_freq [,h_s:=h_s]
    al_freq [,h_w:=h_w]
    al_freq [,s_s:=s_s]
    al_freq [,s_w:=s_w]
    al_freq [,fit:=fiton]
    al_freq [,s_gen:=sum_gen]
    al_freq [,w_gen:=win_gen]
    al_freq [,s_pop:=sum_pop]
    al_freq [,w_pop:=win_pop]
    al_freq[, pop_season := pop_s]
    al_freq[, gen_season := gen_s]
    freq_data  = rbind(freq_data , al_freq )
  }

  # collate sum_stat files
  for (i in 1:length(s_list)){
     filename = s_list [i]
    s_stat  = fread(file = paste0(path, sim_type,"/group_",g, "/",filename))
    s_stat[,replicate:=i]
    s_stat[,label:=g_label]
    s_stat[,group:=g]
    s_stat [,h_s:=h_s]
    s_stat [,h_w:=h_w]
    s_stat [,s_s:=s_s]
    s_stat [,s_w:=s_w]
    s_stat [,fit:=fiton]
    s_stat [,s_gen:=sum_gen]
    s_stat [,w_gen:=win_gen]
    s_stat [,s_pop:=sum_pop]
    s_stat [,w_pop:=win_pop]
    s_stat[, pop_season := pop_s]
    s_stat[, gen_season := gen_s]
    s_stat[, mutRate:=mr]
    s_stat[, recRate:=rr]
    s_stat[, genomeSize := genomeSize]
    sum_stats = rbind(sum_stats, s_stat,fill=TRUE)
  }
}

freq_data [, mean_freq:=mean(mut_freq), by=c("Gen", "mut_pos", "group")] ## calculate mean allele frequency
freq_data [, block :=paste0(group, "_", replicate)] ## used to isolate lines in plots
freq_data[gen_season == "EG", gen_year:=Gen%%20, by=c("label", "Gen")] ## the generation in a seasonal cycle/year


sum_stats[, midpoint:= (win_end-win_start)/2 + win_start] ## midpoint of a window on which statistics are calculated
sum_stats[, block :=paste0(group, "_", replicate)] ## used to isolate lines in plots
sum_stats[,burninNe:=((s_gen+w_gen)/(((1/s_pop)*s_gen)+((1/w_pop)*w_gen))), by="group"] ## harmonic mean Ne used for burn-in
sum_stats[, exp:=4*as.numeric(burninNe)*as.numeric(mutRate), by="group"] ## Expected theta (4Neu)
setnames(sum_stats, "time", "Gen")
sum_stats[gen_season == "EG", gen_year:=Gen%%20, by=c("label", "Gen")]
sum_stats[gen_season == "UG", gen_year:=Gen%%12, by=c("label", "Gen")] 
sum_stats[gen_season == "EG", season:=ifelse((gen_year<11 & gen_year>0), "summer", "winter"), by=c("Gen")] ## season (summer/winter)
sum_stats[gen_season == "UG", season:=ifelse((gen_year<11 & gen_year>0), "summer", "winter"), by=c("Gen")]

sum_stats[, edges:=ifelse(n_win==0 |n_win==500|n_win==510, T, F), by="n_win"] ## identified windows where linkage breaks are/ half windows
sum_stats[,linkage:=ifelse(midpoint>genomeSize, "unlinked", "linked"), by=c("group", "n_win")] ## Identify windows linked to selected site
sum_stats[, dist:=n_win-250]  ## distance from the selected site

## Save as R data objects
save(sum_stats, file=paste0(path,"single_locus_", sim_type,".RData"))
save(freq_data, file=paste0(path,"af_single_locus_", sim_type,".RData"))
