### Haplotype Frequency Spectrum ###

if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, ggplot2)

## set path and working directory
path='/Path/to/FootprintsofFluctuatingSelection/'
setwd(path)

types=c("hard", "soft", "neutral", "fluctuating", "heterozygote") # selection types
groups=c(2,2,1,2,2) # unique parameter set identifiers (1 for each selection type)
hdata=NULL
for (s in 1:length(types)){ # loop through selection types
  sim_type=types[s]
  g=groups[s]
  # extract haplotype frequency files
  hfs_list  <- list.files(path =paste0(path, sim_type, "/group_", g, "/"),pattern =paste0("hfs_",g)) 
  
  for (i in 1:50){ # loop through 50 replicates
    rep=fread(paste0(path, sim_type, "/group_", g, "/",hfs_list[i])) # read in haplotype file
    rep[, rep:=i] # label replicate
    eq=rep[, min(time)] # extract earliest timepoint to subset data (should only be used on short sims so earliest is time of equilibrium or fixation)
    eq_hap=rep[time==eq&n_win==250] # subset early equilirbium/fixation data
    eq_hap[, sim_type:=sim_type] # add selection type label
    hdata=rbind(hdata,eq_hap) # bind to hdata object
    if (sim_type=="fluctuating"|sim_type=="balancing"|sim_type=="neutral"){ ## for selection types sampled twice extract long-term data
      lt_hap=rep[time==96050&n_win==250]
      lt_hap[, sim_type:=paste0("lt_",sim_type)]
      hdata=rbind(hdata, lt_hap)}
      
    }
  }

# label selection types
hdata[sim_type=="balancing", label:="Early Equilibrium HA"]
hdata[sim_type=="wittmann_unlinked", label:="Early Equilibrium FS"]
hdata[sim_type=="lt_balancing", label:="Long-Term HA"]
hdata[sim_type=="lt_wittmann_unlinked", label:="Long-Term FS"]
hdata[sim_type=="lt_neutral", label:="Neutral"]
hdata[sim_type=="hard", label:="Hard Sweep"]
hdata[sim_type=="soft", label:="Soft Sweep"]

# identify the count of the 5th most common haplotype
top5=hdata[,tail(sort(count),5), by=c("sim_type", "rep")][, min(V1), by=c("sim_type", "rep")]

# extract the top 5 most common haplotypes for all selection types, timepoitns and replicates
top.5=NULL
  for (s in unique(data$sim_type)){
    for (i in 1:50){
 top=hdata[sim_type==s&count>=top5[sim_type==s& rep==i]$V1 & rep==i] 
 setorder(top, -count)
 top[, id:=1:length(count)]
 top.5=rbind(top.5,top)}
  }
# average count to get frequency of the haplotypes
top.5[,mean:=(mean(count)/50)/200, by=c("sim_type", "id")]

# plot most common haplotype frequency spectrum
hapfreq=ggplot(top.5[id<6], aes(x=id, y=mean, fill=label))+geom_col()+theme(legend.position = "none")+
  labs(x="Most Common Haplotypes", y="Average Haplotype Frequency")+ 
  facet_wrap(~factor(label, c("Neutral", "Early Equilibrium FS","Long-Term FS","Hard Sweep", "Soft Sweep", "Early Equilibrium HA","Long-Term HA"))) 
  
ggsave("plots/FigS4.pdf", plot=hapfreq)
ggsave("plots/FigS4.jpg", plot=hapfreq)
