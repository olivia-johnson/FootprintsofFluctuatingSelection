## Comparing all selection forms

if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, ggplot2, stringr, egg, dplyr, rstatix, ggpubr, lsr, ggnewscale)

path="/Path/to/FootprintsofFluctuatingSelection/"
setwd(path)
## Labels used to subset selection data
p_label="CP10k_EG_h_0_s_0.5"  ## parameter label for positive selection sims
w_label="CP10k_EG_h_0.6_s_0.5" ## parameter label for fluctutaing selection sims

load(file=paste0(path,"single_locus_neutral.RData")) ## load R data object
sum_stats=sum_stats[label=="CP10k_EG_No Fitness"] ## subset for conditions of interest
bins=seq(from=0,to= 25*((max(sum_stats$n_seg_sites)%/%25)+1), by=25)# create bins of number of segregating sites
sum_stats[, bin:=cut(n_seg_sites, bins)]# bin windows  by number of segregating sites
#calculate mean and sd to use to standardise NCD based on number of segregating sites per window, as suggested in Bitarello et al. 2018.
ncd_stand=sum_stats[, .(m.ncd3=mean(ncd_3), m.ncd4=mean(ncd_4), m.ncd5=mean(ncd_5), sd.ncd3=sd(ncd_3), 
sd.ncd4=sd(ncd_4), sd.ncd5=sd(ncd_5)), by="bin"]
sum_stats=merge(sum_stats, ncd_stand, by="bin") ## merge to full dataset

sum_stats[,ncd5_z:=(ncd_5-m.ncd5)/sd.ncd5, by=c("bin", "ncd_5")] # standardise NCD (TF=0.5)
sum_stats[,ncd4_z:=(ncd_4-m.ncd4)/sd.ncd4, by=c("bin", "ncd_4")] # standardise NCD (TF=0.4)
sum_stats[,ncd3_z:=(ncd_3-m.ncd3)/sd.ncd3, by=c("bin", "ncd_3")] # standardise NCD (TF=0.3)

sum_stats[, year:=ifelse(gen_season=="EG",Gen%/%20, Gen%/%12), by="Gen"] # Determine the year/seasonal cycle of each generation
gens=c(96040:96060) ## Finals sampling generations
neutral=sum_stats[Gen %in% gens & edges==F &(gen_year==5|gen_year==10|gen_year==15|gen_year==0)] ## subset to just generations of interest
neutral[, sim_type:=ifelse(Gen<2000, "Early Neutral", "Neutral"), by="Gen"] ## add sim_type label



#load in soft sweep data, recalculate ncd, sample generation following fixation of allele
load(file=paste0(path,"single_locus_soft.RData"))
bins=seq(from=0,to= 25*((max(sum_stats$n_seg_sites)%/%25)+1), by=25)
sum_stats[, bin:=cut(n_seg_sites, bins)]
sum_stats=merge(sum_stats, ncd_stand, by="bin")

sum_stats[,ncd5_z:=(ncd_5-m.ncd5)/sd.ncd5, by=c("bin", "ncd_5")]
sum_stats[,ncd4_z:=(ncd_4-m.ncd4)/sd.ncd4, by=c("bin", "ncd_4")]
sum_stats[,ncd3_z:=(ncd_3-m.ncd3)/sd.ncd3, by=c("bin", "ncd_3")]
sum_stats[,ncd_z:=(ncd-m.ncd)/sd.ncd, by=c("bin", "ncd")]

sum_stats[, year:=ifelse(gen_season=="EG",Gen%/%20, Gen%/%12), by="Gen"]

## conditional sampling based on time of allele fixation for each replicate
load(file=paste0(path,"af_single_locus_soft.RData"))
fix_times = freq_data[, sum(mut_freq)==1, by=c("Gen","label", "run")][V1==TRUE, min(Gen), by=c("label", "run")]
soft_data=NULL
start_values=c(cumprod(c(20, rep.int(2, 9))), 12000, 24000, 36000,48000,60000,72000,84000,96000)
  for (j in 1:50){ ## loop through 50 replicates
    fixed=fix_times[label==p_label & run==j, V1]
    values=sort(unique(sum_stats[label==p_label & run==j]$Gen))
    if (fixed %in% values){
      if ((fixed) %in% values ){
        
        sum_soft=sum_stats[label==p_label & run==j & Gen==fixed & edges==F ]
      }else{
        val=start_values[start_values>fixed]
        nearest=min(val)
        sum_soft=sum_stats[label==p_label & run==j & Gen==nearest  &edges==F ]
      }}else {
        val1=min(values[values>fixed])
        val2=min(start_values[start_values>fixed])
        val=min(val1, val2)
        nearest=min(val)
        sum_soft=sum_stats[label==p_label & run==j & Gen==nearest&edges==F]
        
      }
    soft_data=rbind(sum_soft,soft_data)
}
soft=soft_data[, sim_type:="Soft Sweep"] ## add label
soft[,gen_year:=10] ## add gen year that will be used to subset later on
soft_af=freq_data[label==p_label] ## subset allele frequency data
soft_af[, sim_type:="Soft Sweep"]  ## add label
soft_af[, year:=ifelse(gen_season=="EG",Gen%/%20, Gen%/%12), by="Gen"] 
soft_af[,gen_year:=10]


#load in hard sweep data, recalculate ncd, sample generation following fixation of allele (the same as for soft sweeps)
load(file=paste0(path,"single_locus_hard.RData"))

bins=seq(from=0,to= 25*((max(sum_stats$n_seg_sites)%/%25)+1), by=25)
sum_stats[, bin:=cut(n_seg_sites, bins)]
sum_stats=merge(sum_stats, ncd_stand, by="bin")

sum_stats[,ncd5_z:=(ncd_5-m.ncd5)/sd.ncd5, by=c("bin", "ncd_5")]
sum_stats[,ncd4_z:=(ncd_4-m.ncd4)/sd.ncd4, by=c("bin", "ncd_4")]
sum_stats[,ncd3_z:=(ncd_3-m.ncd3)/sd.ncd3, by=c("bin", "ncd_3")]
sum_stats[,ncd_z:=(ncd-m.ncd)/sd.ncd, by=c("bin", "ncd")]

sum_stats[, year:=ifelse(gen_season=="EG",Gen%/%20, Gen%/%12), by="Gen"]


load(file=paste0(path,"af_single_locus_hard.RData"))
fix_times = freq_data[, sum(mut_freq)==1, by=c("Gen","label", "run")][V1==TRUE, min(Gen), by=c("label", "run")]
hard_data=NULL
start_values=c(cumprod(c(20, rep.int(2, 9))), 12000, 24000, 36000,48000,60000,72000,84000,96000)

for (j in 1:50){
  fixed=fix_times[label==p_label & run==j, V1]
  values=sort(unique(sum_stats[label==p_label & run==j]$Gen))
  if (fixed %in% values){
    if ((fixed) %in% values ){
      
      sum_hard=sum_stats[label==p_label & run==j & Gen==fixed & edges==F ]
    }else{
      val=start_values[start_values>fixed]
      nearest=min(val)
      sum_hard=sum_stats[label==p_label & run==j & Gen==nearest  &edges==F ]
    }}else {
      val1=min(values[values>fixed])
      val2=min(start_values[start_values>fixed])
      val=min(val1, val2)
      nearest=min(val)
      sum_hard=sum_stats[label==p_label & run==j & Gen==nearest&edges==F]
      
    }
  
  hard_data=rbind(sum_hard,hard_data)
}
hard=hard_data[, sim_type:="Hard Sweep"]
hard[,gen_year:=10]
hard_af=freq_data[label==p_label]
hard_af[, sim_type:="Hard Sweep"]
hard_af[, year:=ifelse(gen_season=="EG",Gen%/%20, Gen%/%12), by="Gen"]
hard_af[,gen_year:=10]

#load in fluctuating data, recalculate ncd, sample one complete year of data (1 summer/winter cycle) following at end of sim (96040 gens; long-term) or early equilibrium
load(file=paste0(path,"single_locus_wittmann_unlinked.RData"))

bins=seq(from=0,to= 25*((max(sum_stats$n_seg_sites)%/%25)+1), by=25)
sum_stats[, bin:=cut(n_seg_sites, bins)]
sum_stats=merge(sum_stats, ncd_stand, by="bin")

sum_stats[,ncd5_z:=(ncd_5-m.ncd5)/sd.ncd5, by=c("bin", "ncd_5")]
sum_stats[,ncd4_z:=(ncd_4-m.ncd4)/sd.ncd4, by=c("bin", "ncd_4")]
sum_stats[,ncd3_z:=(ncd_3-m.ncd3)/sd.ncd3, by=c("bin", "ncd_3")]
sum_stats[,ncd_z:=(ncd-m.ncd)/sd.ncd, by=c("bin", "ncd")]

sum_stats[, year:=ifelse(gen_season=="EG",Gen%/%20, Gen%/%12), by="Gen"]

load(file=paste0(path,"af_single_locus_wittmann_unlinked.RData"))

wittmann=sum_stats[label==w_label& Gen >=96040 & Gen<96060] ## sample long-term timepoint
wittmann[, sim_type:="Long-Term Fluctuating"] ## add label
wittmann_af=freq_data[label==w_label] ## subset allele frequency data
wittmann_af[, sim_type:="Fluctuating"]
wittmann_af[, year:=ifelse(gen_season=="EG",Gen%/%20, Gen%/%12), by="Gen"]
e.year=wittmann_af[gen_year==5|gen_year==15, mean(mut_freq), by=c("year", "run")][V1>0.45, min(year), by="run"] ## determine year when simulation reaches equilibrium (when mean frequency acorss the season is > 0.45)
## conditional sampling to extract data for the year following the alelle frequency reaching equilibrium
wittmanne=NULL
for (j in 1:50){
  eq=e.year[ run==j, V1]
  values=sort(unique(sum_stats[label==w_label & run==j]$year))
  if (eq %in% values){
      sum_wit=sum_stats[label==w_label & run==j & year==eq &edges==F ]
  if (length(sum_wit[,unique(Gen)])!=20){
    val=values[values>eq]
    nearest=min(val)
    sum_wit=sum_stats[label==w_label & run==j & year==nearest &edges==F ]}}else{
      val=values[values>eq]
      nearest=min(val)
      sum_wit=sum_stats[label==w_label & run==j & year==nearest &edges==F ]}
  wittmanne=rbind(sum_wit, wittmanne)}
if (sum(wittmanne[,length(unique(Gen)), by="run"]$V1==20)==50){print("Full season")}else{print("Not full season! - Wittmann")}
wittmanne[, sim_type:="Early Equilibrium Fluctuating"]

#load in balancing data, recalculate ncd, sample one complete year of data (1 summer/winter cycle) following at end of sim (96040 gens; long-term) or early equilibrium

load(file=paste0(path,"single_locus_balancing.RData"))

bins=seq(from=0,to= 25*((max(sum_stats$n_seg_sites)%/%25)+1), by=25)
sum_stats[, bin:=cut(n_seg_sites, bins)]
sum_stats=merge(sum_stats, ncd_stand, by="bin")

sum_stats[,ncd5_z:=(ncd_5-m.ncd5)/sd.ncd5, by=c("bin", "ncd_5")]
sum_stats[,ncd4_z:=(ncd_4-m.ncd4)/sd.ncd4, by=c("bin", "ncd_4")]
sum_stats[,ncd3_z:=(ncd_3-m.ncd3)/sd.ncd3, by=c("bin", "ncd_3")]

sum_stats[, year:=ifelse(gen_season=="EG",Gen%/%20, Gen%/%12), by="Gen"]

load(file=paste0(path, "af_single_locus_balancing.RData"))

balancing=sum_stats[label=="CP10k_EG_h_0_s_0.1"&group==2& Gen >=96040 & Gen<96060] ## subset s = 0.1 sims
balancing[, sim_type:="Long-Term Balancing"]
balancing_af=freq_data[label=="CP10k_EG_h_0_s_0.1"]
balancing_af[, sim_type:="Balancing"]
balancing_af[, year:=ifelse(gen_season=="EG",Gen%/%20, Gen%/%12), by="Gen"]
## conditional sampling to extract data for year after equilibrium frequency (>0.49) is reached
e.year=balancing_af[gen_year==5|gen_year==15, mean(mut_freq), by=c("year", "run")][V1>0.49, min(year), by="run"]
balancinge=NULL
for (j in 1:50){
  eq=e.year[ run==j, V1]
  values=sort(unique(sum_stats[label=="CP10k_EG_h_0_s_0.1" & run==j]$year))
  if (eq %in% values){
    sum_bal=sum_stats[label=="CP10k_EG_h_0_s_0.1" & run==j & year==eq &edges==F ]
    if (length(sum_bal[,unique(Gen)])!=20){
      val=values[values>eq]
      nearest=min(val)
      sum_bal=sum_stats[label=="CP10k_EG_h_0_s_0.1" & run==j & year==nearest &edges==F ]}}else{
        val=values[values>eq]
        nearest=min(val)
        sum_bal=sum_stats[label=="CP10k_EG_h_0_s_0.1" & run==j & year==nearest &edges==F ]}
  balancinge=rbind(sum_bal, balancinge)}

 if (sum(balancinge[,length(unique(Gen)), by="run"]$V1==20)==50){print("Full season")}else{print("Not full season! - Balancing")}

balancinge[, sim_type:="Early Equilibrium Balancing"]

## Merge datasets into a single complete dataset
sum_stat=rbind(balancing, wittmann, soft, hard, neutral,balancinge, wittmanne)
freq_data=rbind(balancing_af, wittmann_af, soft_af, hard_af)

sum_stat[, dist:=n_win-250] ## calculate distance from selected site
sum_stat[, block:=paste0(run,"_", Gen), by=c("sim_type", "run", "Gen")] ## a groupng variable to separate replicate lines in pltos

data=sum_stat[edges==F&(gen_year==10)]#subset end of summer generation
## determine timepoint in season
data[gen_year==0 | gen_year==10, fluctuation:=ifelse(gen_year==10, "End Summer", "End Winter"), by="gen_year"] 
data[gen_year==5 | gen_year==15, fluctuation:=ifelse(gen_year==5, "Mid Summer", "Mid Winter"), by="gen_year"]
## calculate mean of statistics for each sim_type, window and timepoint
data[, rel_div:=diversity/exp, by=c("sim_type", "run", "Gen")]
data[, mean_tajd:=mean(tajimas_d_branch), by=c( "sim_type","midpoint", "fluctuation")]
data[, mean_div:=mean(diversity), by=c( "sim_type","midpoint", "fluctuation")]
data[, mean_thetaw:=mean(theta_w_allele), by=c("sim_type", "midpoint", "fluctuation")]
data[, mean_H1:=mean(H1), by=c( "sim_type","midpoint", "fluctuation")]
data[, mean_H12:=mean(H12), by=c("sim_type", "midpoint", "fluctuation")]
data[, mean_H123:=mean(H123), by=c("sim_type", "midpoint", "fluctuation")]
data[, mean_H2H1:=mean(H2H1), by=c( "sim_type","midpoint", "fluctuation")]
data[, mean_ncd_5:=mean(ncd_5), by=c( "sim_type","midpoint", "fluctuation")]
data[, mean_ncd5_z:=mean(ncd5_z), by=c( "sim_type","midpoint", "fluctuation")]
data[, mean_ncd_4:=mean(ncd_4), by=c( "sim_type","midpoint", "fluctuation")]
data[, mean_ncd4_z:=mean(ncd4_z), by=c( "sim_type","midpoint", "fluctuation")]
data[, mean_ncd_3:=mean(ncd_3), by=c( "sim_type","midpoint", "fluctuation")]
data[, mean_ncd3_z:=mean(ncd3_z), by=c( "sim_type","midpoint", "fluctuation")]
data[, mean_var:=mean(variance), by=c( "sim_type","midpoint", "fluctuation")]
data[, mean_skew:=mean(skew), by=c( "sim_type","midpoint", "fluctuation")]
data[, mean_kurtosis:=mean(kurtosis), by=c( "sim_type","midpoint", "fluctuation")]

### Compare statistics at central window for selection types using t-test and Hockberg correction
significance=NULL
for (i in c(10)){ ## can cycle through seasonal timepoint but set to end of summer
  sig_data=data[n_win==250]
  setorder(sig_data, -sim_type)
  sig = as.data.table(compare_means(c(diversity, theta_w_allele,tajimas_d_site, H1, H12, H2H1,H123, ncd_3, ncd_4, ncd_5, ncd5_z,
              ncd4_z,ncd3_z, variance, skew, kurtosis)~sim_type, data = sig_data, method = "t.test", p.adjust.method="hochberg"))
  sig[, gen_year:=i]
  significance=rbind(significance, sig)
  }
print(sig_data[, .N, by="sim_type"])
## adjust stars to represent Hochberg adjusted p-value
significance=significance%>%mutate(p.adj.sig = ifelse(p.adj < 0.001, "***",
                                                       ifelse(p.adj < 0.01, "**",
                                                              ifelse(p.adj < 0.05,
                                                                     "*", "ns"))))

## Calculate Cohen's D from significance dataset
comps=unique(significance[, .(group1, group2, .y.=as.character(.y.))]) ##.y. is the statistic name
comps[, cohens_d:=0]
coh_data=data[dist==0]
for (i in 1:length(comps$group1)){
  group.1=comps$group1[i]
  group.2=comps$group2[i]
  statistic=comps$.y.[i]
  g1 = coh_data[sim_type==group.1, .SD, .SDcol=statistic][[1]]
  g2 =coh_data[sim_type==group.2, .SD, .SDcol=statistic][[1]]
  cd=cohensD(g1,g2)
  comps$cohens_d[i]=cd
}
## Merge Cohen's D to significance dataset
significance=merge(significance, comps, by = c("group1", "group2", ".y."))

## Investigate significance dataset
##Example comparing fluctuating with positive selection
significance[ p.adj<0.05 & gen_year==10 &(group1=="Early Equilibrium Fluctuating" |group1=="Long-Term Fluctuating" |
  group2=="Early Equilibrium Fluctuating" |group2=="Long-Term Fluctuating")& (group1=="Soft Sweep" |group1=="Hard Sweep" |
  group2=="Soft Sweep" |group2=="Hard Sweep"), .N, by=c(".y.")] ## number of significant comparisons

significance[ p.adj<0.05 &(group1=="Long-Term Fluctuating" |group2=="Long-Term Fluctuating")& (group1=="Hard Sweep" |
  group2=="Hard Sweep"), .(max(p.adj), min(cohens_d)), by=".y."]  ## Maximum P value and Minimum Cohen's D for comparisons

## Plots ##

## labels
supp.labs <- c("Nucleotide Diversity", "Watterson's Theta","Tajima's D",  "Garud's H1", "Garud's H12", "Garud's H2/H1","Garud's H123", "NCD (TF=0.5)","NCD (TF=0.4)","NCD (TF=0.3)", "Standardised NCD (TF=0.5)","Standardised NCD (TF=0.4)","Standardised NCD (TF=0.3)",  "Variance", "Skew", "Kurtosis")
names(supp.labs) <- c("diversity", "theta_w_allele","tajimas_d_site",  "H1", "H12", "H2H1","H123","ncd_5","ncd_4","ncd_3", "ncd5_z", "ncd4_z","ncd3_z", "variance", "skew", "kurtosis")
col.names <- c("diversity", "theta_w_allele","tajimas_d_site",  "H1", "H12", "H2H1","H123","ncd_5","ncd_4","ncd_3", "ncd5_z", "ncd4_z","ncd3_z", "variance", "skew", "kurtosis")
data_cols=c("dist", "diversity", "theta_w_allele","tajimas_d_site",  "H1", "H12", "H2H1","H123","ncd_5","ncd_4","ncd_3", "ncd5_z","ncd4_z","ncd3_z",  "variance", "skew", "kurtosis", "sim_type", "run", "gen_year", "fluctuation")

## Figure 5 (S10, S11)
  ## columns to subset data
data_cols=c("dist", "diversity", "theta_w_allele","tajimas_d_site", "H1", "H12", "H123","H2H1","ncd_5", "ncd5_z", "variance", "skew", "kurtosis", "sim_type", "run", "gen_year", "fluctuation")

focal=data[dist==0 , .SD, .SDcols=data_cols] ## subset data
all_comp=melt(focal, id.vars=c("sim_type", "run", "gen_year", "dist", "fluctuation"), variable.name="statistic") ## change shape of tables
plot=ggplot(all_comp[gen_year==10], aes(x=sim_type, y=value))+ theme_bw()+
  geom_boxplot(aes(fill=sim_type), col="gray30") +
  scale_fill_manual(values = c("Neutral" = "#C49A00","Long-Term Balancing"="#F8766D","Early Equilibrium Balancing"="#FB61D7","Long-Term Fluctuating"="#00C094","Early Equilibrium Fluctuating"="#53B400", "Hard Sweep"="#A58AFF", "Soft Sweep"="#00B6EB")) +
  facet_wrap("statistic", scale="free_y", labeller = labeller(statistic = supp.labs), ncol = 3)+  
  scale_x_discrete(limits = c("Neutral", "Early Equilibrium Balancing", "Long-Term Balancing", "Early Equilibrium Fluctuating","Long-Term Fluctuating","Hard Sweep", "Soft Sweep"))+
  labs( y="Value of Statistic", x="Type of Selection")+theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plot
ggexport(plot, filename="plots/boxplot_s01.pdf", width = 8.5, height=9)
ggsave(plot, filename="plots/boxplot_s01.jpg", width = 8.5, height=9)


## Figure 6 - Comparisons between balancing and fluctuating selection

plot_dist=0.25 ## distance on either side of selected site to plot
bdata=data[sim_type=="Early Equilibrium Balancing"& gen_year==10] ## subset data to just early equilibrium balancing
fdata=data[sim_type=="Early Equilibrium Fluctuating"& gen_year==10] ## subset data to just early equilibrium fluctuating
f.col="#53B400" ## colour for early equilibrium fluctuating
b.col="#FB61D7" ## colour for early equilibrium balancing

  ## Early equilibrium panels
f_taj = ggplot(fdata, aes(x=dist/100))+
  geom_line(aes(y=tajimas_d_site,group=block), alpha=0.1, col=f.col)+
  geom_line(aes(y=mean_tajd), col="black", size=0.8)+
  theme_bw()+
  geom_vline(data = subset(fdata),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  facet_wrap("sim_type")+
  theme(axis.title.x.bottom = element_blank(), axis.title.y.left = element_blank())+
  coord_cartesian(y=c(-1, 1)) +
  xlim(-plot_dist,plot_dist) 
b_taj = ggplot(bdata, aes(x=dist/100))+
  geom_line(aes(y=tajimas_d_site,group=block), alpha=0.1, col=b.col)+
  geom_line(aes(y=mean_tajd), col="black", size=0.8)+
  theme_bw()+
  geom_vline(data = subset(bdata),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  facet_wrap("sim_type")+
  theme(axis.title.x.bottom = element_blank())+
  coord_cartesian(y=c(-1, 1)) +
  ylab("Tajima's D")+
  xlim(-plot_dist,plot_dist) 
com_taj = ggarrange(b_taj, f_taj, nrow=1)
f_H1 = ggplot(fdata, aes(x=dist/100))+
  geom_vline(data = subset(fdata),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=H1,group=block), alpha=0.1, col=f.col)+
  geom_line(aes(y=mean_H1), col="black", size=0.8)+
  theme_bw()+
  facet_wrap("sim_type")+  theme(axis.title.x.bottom = element_blank(), axis.title.y.left = element_blank())+
  coord_cartesian(x=c(-plot_dist,plot_dist), y=c(0.002, 0.075))
b_H1 = ggplot(bdata, aes(x=dist/100))+
  geom_vline(data = subset(bdata),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=H1,group=block), alpha=0.1, col=b.col)+
  geom_line(aes(y=mean_H1), col="black", size=0.8)+
  theme_bw()+
  facet_wrap(~sim_type)+
  theme(axis.title.x.bottom = element_blank(),)+
  ylab("Garud's H1")+
  coord_cartesian(x=c(-plot_dist,plot_dist), y=c(0.002, 0.075))
com_H1 = ggarrange(b_H1, f_H1, nrow=1)

  ## long-term panels

bdata=data[sim_type=="Long-Term Balancing"& gen_year==10] ## subset data to just long-term balancing
fdata=data[sim_type=="Long-Term Fluctuating"& gen_year==10] ## subset data to just long-term fluctuating
f.col="#00C094" ## colour for long-term fluctuating
b.col="#F8766D" ## colour for long-term fluctuating

f_div = ggplot(fdata, aes(x=dist/100))+
  geom_vline(data = subset(fdata),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=diversity,group=block), alpha=0.1, col=f.col)+
  geom_line(aes(y=mean_div), col="black", size=0.8)+
  theme_bw()+
  facet_wrap("sim_type")+
  theme(axis.title.x.bottom = element_blank(), axis.title.y.left = element_blank())+
  coord_cartesian(x=c(-.5,.5), y=c(0.0025, 0.0055))
b_div = ggplot(bdata, aes(x=dist/100))+
  geom_vline(data = subset(bdata),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=diversity,group=block), alpha=0.1, col=b.col)+
  geom_line(aes(y=mean_div), col="black", size=0.8)+
  theme_bw()+
  facet_wrap("sim_type")+
  theme(axis.title.x.bottom = element_blank(),)+
  ylab("Nucleotide Diversity")+
  coord_cartesian(x=c(-.5,.5), y=c(0.0025, 0.0055))
com_div = ggarrange(b_div, f_div)

f_ncd5 = ggplot(fdata, aes(x=dist/100))+
  geom_line(aes(y=ncd5_z,group=block), alpha=0.1, col=f.col)+
  geom_line(aes(y=mean_ncd5_z), col="black", size=0.8)+
  theme_bw()+
  geom_vline(data = subset(fdata),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  facet_wrap("sim_type")+
  theme(axis.title.x.bottom = element_blank(), axis.title.y.left = element_blank())+
  xlim(-plot_dist,plot_dist) 
b_ncd5 = ggplot(bdata, aes(x=dist/100))+
  geom_line(aes(y=ncd5_z,group=block), alpha=0.1, col=b.col)+
  geom_line(aes(y=mean_ncd5_z), col="black", size=0.8)+
  theme_bw()+
  geom_vline(data = subset(bdata),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  facet_wrap("sim_type")+
  theme(axis.title.x.bottom = element_blank(),)+
  ylab("Standardised NCD (TF=0.5)")+
  xlim(-plot_dist,plot_dist) 
com_ncd5 = ggarrange(b_ncd5, f_ncd5, nrow=1)
  ## combine significant panels
sig=ggarrange(com_taj, com_div, ncol=1)
sig=annotate_figure(sig, left="Significantly Different")
  ## combine non-significant panels
nonsig=ggarrange(com_H1, com_ncd5, ncol=1)
nonsig=annotate_figure(nonsig, left="Not Significantly Different")
  ## combine all panels
div_fb=ggarrange(sig, nonsig, ncol=1)
div_fb=annotate_figure(div_fb, bottom="Distance from selected site (Mb)")
ggexport(div_fb, filename="plots/Fig6.pdf", width=7, height = 8)
ggsave(filename="plots/Fig6.jpg", div_fb, width=7, height = 8)

## SI heatmaps
## breaks for p-value scale 
col.labs <- data.frame(col.names, supp.labs, row.names = NULL)
col.labs$supp.labs<-factor(col.labs$supp.labs, levels=c(col.labs$supp.labs))

heatdata=significance[.y.!="ncd3_z" &.y.!="ncd4_z"&.y.!="ncd_3" &.y.!="ncd_4"] ## exclude NCD (TF = 0.4 & 0.3)
heatdata=merge(heatdata, col.labs, by.x = ".y.", by.y="col.names")
heat=ggplot()+
  geom_tile(data=heatdata[p.adj<0.05],aes(x=group1, y=group2, fill=p.adj))+
  scale_fill_gradientn(trans="log10", na.value = "white",colours=c("orange", "darkmagenta", "navyblue"))+
  geom_text(data=heatdata[p.adj>0.05],aes(x=group1, y=group2,label=p.adj.sig), col="black")+
  labs(fill="Adjusted p-value \n(p < 0.05)")+
  new_scale_fill()+
  geom_tile(data=heatdata, aes( x=group2, y=group1, fill=cohens_d))+
  labs(fill="Cohen's D")+
  facet_wrap("supp.labs", ncol=3)+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1), axis.title = element_blank())
heat
ggsave(plot=heat, filename="plots/heatmap_s05.jpg", height = 9, width=8.5)
ggsave(plot=heat, filename="plots/heatmap_s05.pdf", height = 9, width=8.5)
