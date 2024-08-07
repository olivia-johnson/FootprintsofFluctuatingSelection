## Linear Discriminant Analysis

if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, ggplot2, stringr, egg, dplyr, ggpubr, MASS, klaR)

path='/Path/to/FootprintsofFluctuatingSelection/'

### Collate additional simulation replicates for LDA ##

sim_type="neutral" # selection type to collate
setwd(paste0(path, sim_type, "/"))

groups=c(1:6) ## unqiue parameter set IDs to cycle through

freq_data = NULL
sum_stats = NULL

for (g in groups){
  print(g)
  
  f_list  <- list.files(path =paste0(path, sim_type,"/group_", g, "/"),pattern ="al_freq_")
  
  s_list  <- list.files(path =paste0(path, sim_type,"/group_", g, "/"),pattern ="sim_stat_")
  
  parameters <- fread(file=paste0(path, sim_type,"/group_",g, ".txt"), sep = ":")
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
    al_freq[,run:=i]
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
    s_stat[,run:=i]
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

freq_data [, mean_freq:=mean(mut_freq), by=c("Gen", "mut_pos", "group")] ## calculate mean allele feqyency for each mutation
freq_data [, block :=paste0(group, "_", run)] ## grouping to separate replicates
freq_data[gen_season == "EG", gen_year:=Gen%%20, by=c("label", "Gen")] ## assign point in the seasonal cycle/year


sum_stats[, midpoint:= (win_end-win_start)/2 + win_start] ## midpoint of the window
sum_stats[, block :=paste0(group, "_", run)] ## grouping to separate replicates
sum_stats[,burninNe:=((s_gen+w_gen)/(((1/s_pop)*s_gen)+((1/w_pop)*w_gen))), by="group"] ## calculate the harmonic mean Ne used for the burn-in
sum_stats[, exp:=4*as.numeric(burninNe)*as.numeric(mutRate), by="group"] ## expected theta (4Neu)
setnames(sum_stats, "time", "Gen")
sum_stats[gen_season == "EG", gen_year:=Gen%%20, by=c("label", "Gen")] ## assign point in the seasonal cycle/year
sum_stats[gen_season == "UG", gen_year:=Gen%%12, by=c("label", "Gen")] ## assign point in the seasonal cycle/year
sum_stats[gen_season == "EG", season:=ifelse((gen_year<11 & gen_year>0), "summer", "winter"), by=c("Gen")] # assign season
sum_stats[gen_season == "UG", season:=ifelse((gen_year<11 & gen_year>0), "summer", "winter"), by=c("Gen")] # asign season

sum_stats[, edges:=ifelse(n_win==0 |n_win==500|n_win==510, T, F), by="n_win"] ## label half windows and wndows over recombination break
sum_stats[, dist:=n_win-250] ## distance from selected site

## save R data objects
save(sum_stats, file=paste0(path,"200_single_locus_", sim_type,".RData"))
save(freq_data, file=paste0(path, "200_af_single_locus_", sim_type,".RData"))


# Extract data from replicates for LDA#

w_label="CP10k_EG_h_0.6_s_0.5" ## label of fluctuating selection parameters to be used in LDA

## load neutral data
load(file=paste0(path,"200_single_locus_neutral.RData"))
## values for NCD standardisation (as in previous R scripts)
ggplot(sum_stats)+ geom_density(aes(x=n_seg_sites, y=..density..))
bins=seq(from=0,to= 25*((max(sum_stats$n_seg_sites)%/%25)+1), by=25)
sum_stats[, bin:=cut(n_seg_sites, bins)]
ncd_stand=sum_stats[, .(m.ncd3=mean(ncd_3), m.ncd4=mean(ncd_4), m.ncd5=mean(ncd_5), sd.ncd3=sd(ncd_3), sd.ncd4=sd(ncd_4), sd.ncd5=sd(ncd_5), m.ncd=mean(ncd), sd.ncd=sd(ncd)), by="bin"]
sum_stats=merge(sum_stats, ncd_stand, by="bin")

sum_stats[,ncd5_z:=(ncd_5-m.ncd5)/sd.ncd5, by=c("bin", "ncd_5")]
sum_stats[,ncd4_z:=(ncd_4-m.ncd4)/sd.ncd4, by=c("bin", "ncd_4")]
sum_stats[,ncd3_z:=(ncd_3-m.ncd3)/sd.ncd3, by=c("bin", "ncd_3")]
sum_stats[,ncd_z:=(ncd-m.ncd)/sd.ncd, by=c("bin", "ncd")]

sum_stats[, year:=ifelse(gen_season=="EG",Gen%/%20, Gen%/%12), by="Gen"]

neutral=sum_stats[Gen >=96040 & Gen<96060 & edges==F &(gen_year==5|gen_year==10|gen_year==15|gen_year==0)& group==1] ## subset neutral data across the season
neutral[, sim_type:="Neutral"]

## load balancing data
load(file=paste0(path,"200_single_locus_balancing.RData")) 

## standardise NCD
bins=seq(from=0,to= 25*((max(sum_stats$n_seg_sites)%/%25)+1), by=25)
sum_stats[, bin:=cut(n_seg_sites, bins)]
sum_stats=merge(sum_stats, ncd_stand, by="bin")
sum_stats[, `:=` (ncd=as.numeric(ncd),ncd_5=as.numeric(ncd_5),ncd_4=as.numeric(ncd_4),
                  ncd_3=as.numeric(ncd_3),skew=as.numeric(skew),kurtosis=as.numeric(kurtosis),
                  variance=as.numeric(variance))]

sum_stats[,ncd5_z:=(ncd_5-m.ncd5)/sd.ncd5, by=c("bin", "ncd_5")]
sum_stats[,ncd4_z:=(ncd_4-m.ncd4)/sd.ncd4, by=c("bin", "ncd_4")]
sum_stats[,ncd3_z:=(ncd_3-m.ncd3)/sd.ncd3, by=c("bin", "ncd_3")]
sum_stats[,ncd_z:=(ncd-m.ncd)/sd.ncd, by=c("bin", "ncd")]

sum_stats[, year:=ifelse(gen_season=="EG",Gen%/%20, Gen%/%12), by="Gen"]

load(file=paste0(path,"200_af_single_locus_balancing.RData"))

## subset long-term and early equilibrium data as in previous R scripts
balancing=sum_stats[label=="CP10k_EG_h_0_s_0.1" & Gen >=96040 & Gen<96060]
balancing[, sim_type:="Long-Term\nBalancing"]
balancing_af=freq_data[label=="CP10k_EG_h_0_s_0.1"]
balancing_af[, sim_type:="Balancing"]
balancing_af[, year:=ifelse(gen_season=="EG",Gen%/%20, Gen%/%12), by="Gen"]
e.year=balancing_af[gen_year==5|gen_year==15, mean(mut_freq), by=c("year", "run")][V1>0.49, min(year), by="run"]
balancinge=NULL
for (j in 1:200){ ## loop through 200 replicates
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
if (sum(balancinge[,length(unique(Gen)), by="run"]$V1==20)==200){print("Full season")}else{print("Not full season! - Balancing")}
balancinge[, sim_type:="Early Equilibrium\nBalancing"]


#load in fluctuating data, recalculate ncd, sample one complete year of data (1 summer/winter cycle) following at end of sim (96040 gens) or early
load(file=paste0(path, "200_single_locus_wittmann_unlinked.RData"))

bins=seq(from=0,to= 25*((max(sum_stats$n_seg_sites)%/%25)+1), by=25)
sum_stats[, bin:=cut(n_seg_sites, bins)]
sum_stats=merge(sum_stats, ncd_stand, by="bin")

sum_stats[,ncd5_z:=(ncd_5-m.ncd5)/sd.ncd5, by=c("bin", "ncd_5")]
sum_stats[,ncd4_z:=(ncd_4-m.ncd4)/sd.ncd4, by=c("bin", "ncd_4")]
sum_stats[,ncd3_z:=(ncd_3-m.ncd3)/sd.ncd3, by=c("bin", "ncd_3")]
sum_stats[,ncd_z:=(ncd-m.ncd)/sd.ncd, by=c("bin", "ncd")]

sum_stats[, year:=ifelse(gen_season=="EG",Gen%/%20, Gen%/%12), by="Gen"]

load(file=paste0(path,"200_af_single_locus_wittmann_unlinked.RData"))

gens=c(10,110,650, 2570, 5130, 12010, 24010, 36010,48010,60010,72010,84010, 96045, 96050, 69055,96060)

wittmann=sum_stats[label==w_label & Gen>=96040 & Gen<96060]#Gen >=96040 & Gen<96060]

wittmann[, sim_type:=ifelse(Gen>96000, "Long-Term\nFluctuating", paste0("Fluctuating_", Gen)), by="Gen"]
wittmann_af=freq_data[label==w_label]
wittmann_af[, sim_type:="Fluctuating"]
wittmann_af[, year:=ifelse(gen_season=="EG",Gen%/%20, Gen%/%12), by="Gen"]
e.year=wittmann_af[gen_year==5|gen_year==15, mean(mut_freq), by=c("year", "run")][V1>0.45, min(year), by="run"]
wittmanne=NULL
for (j in 1:200){
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
if (sum(wittmanne[,length(unique(Gen)), by="run"]$V1==20)==200){print("Full season")}else{print("Not full season! - Wittmann")}
wittmanne[, sim_type:="Early Equilibrium\nFluctuating"]

remove(wittmann_af, balancing_af)  ## remove allele frequency files to make space

## add sel_type label
neutral[, sel_type:="Neutral"]
balancing[, sel_type:="Balancing"]
balancinge[, sel_type:="Balancing"]
wittmann[, sel_type:="Fluctuating"]
wittmanne[, sel_type:="Fluctuating"]

## Merge dataset
data=rbind(balancing,  wittmann[sim_type=="Long-Term\nFluctuating"], neutral,balancinge, wittmanne)

### RUN STEP-WISE LDAS ###

## s = 0.1

# Generate training and test datasets of 150 and 50 replicates, respectively.
train.data <- data[ n_win==250 & gen_year==10& run>50]
test.data <- data[ n_win==250 & gen_year==10&run<=50]

# Check correct number of replicates
train.data[, .N, by="sim_type"]
test.data[, .N, by="sim_type"]

# Assign number of categories (simulation types)
numberofCat=length(unique(train.data$sim_type))

# Formula for step-wise analysis (all statistis to be considered)
formulaAll=sim_type~diversity+tajimas_d_site+theta_w_allele+ncd_5+ncd_4+ncd_3+ncd5_z+ncd4_z+ncd3_z+variance+skew+kurtosis+H1+H12+H123+H2H1

# Generate formula for stepwise LDA
greedy.wilks(formulaAll,data=train.data, niveau = 0.05) 

# Copy formula from output of above function
formulaStepwise=sim_type ~ H2H1 + diversity + kurtosis + H1 + H123 + tajimas_d_site +  ncd_3 + variance

# Generate LDA model from training dataset using stepwise formula
model <- lda(formulaStepwise,data=train.data, method="moment")

# Predict selection types for test dataset using LDA model
predictions <- model %>% predict(test.data, CV=TRUE)

# Table of predictions
pt=as.data.table(table(pred=predictions$class, true=test.data$sim_type))

# Calculate accuracy proportions (N correct calls /50)
pt[pred=="Equilibrium\nFluctuating" & true=="Equilibrium\nFluctuating", N]/50
pt[pred=="Long-Term Fluctuating" & true=="Long-Term Fluctuating", N]/50
pt[pred=="Equilibrium Balancing" & true=="Equilibrium Balancing", N]/50
pt[pred=="Long-Term Balancing" & true=="Long-Term Balancing", N]/50
pt[pred=="Neutral" & true=="Neutral", N]/50
pt[, Freq:=N/50, by=c("pred", "true")]

overall_cw=mean(predictions$class==test.data$sim_type)

# Generate panel for confusion matrix figure
cw01=ggplot(pt, aes(pred,true, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#F8766D") +
  labs( title="Central window") + 
  theme(legend.position = "None", axis.title.x = element_blank(), axis.title.y = element_blank() )

# vector or statistics
stats=c("diversity","tajimas_d_site","theta_w_allele","ncd_5","ncd_4","ncd_3","ncd5_z","ncd4_z","ncd3_z",
        "variance","skew","kurtosis","H1","H12","H123","H2H1")

# LDA using data from the end of seasons
train.data <- data[ n_win==250 & (gen_year==10 |gen_year==0) & run>50] ## subset training dataset using sampled from end of seasons
train.data=dcast(train.data, sel_type+sim_type+year+run+label ~ gen_year, value.var = stats) ## transform dataset for LDA
test.data <- data[ n_win==250 & (gen_year==10 |gen_year==0)&run<=50] # subset test dataset using sampled from end of seasons
test.data=dcast(test.data, sel_type+sim_type+year+run+label ~ gen_year, value.var = stats) ## transform dataset for LDA

# Number of categories
numberofCat=length(unique(train.data$sim_type))

# Formual of all stistics to be used for step-wise grouping
formulaAll=sim_type~diversity_0+tajimas_d_site_0+theta_w_allele_0+ncd_5_0+ncd_4_0+ncd_3_0+ncd5_z_0+ncd4_z_0+ncd3_z_0+
  variance_0+skew_0+kurtosis_0+H1_0+H12_0+H123_0+H2H1_0+diversity_10+tajimas_d_site_10+theta_w_allele_10+ncd_5_10+ncd_4_10+
  ncd_3_10+ncd5_z_10+ncd4_z_10+ncd3_z_10+variance_10+skew_10+kurtosis_10+H1_10+H12_10+H123_10+H2H1_10

# Run stepwise analysis
greedy.wilks(formulaAll,data=train.data, niveau = 0.05) 

# Copy formula from output of above function
formulaStepwise=sim_type ~ H2H1_10 + diversity_10 + kurtosis_10 + H1_10 + H123_10 + 
  H2H1_0 + diversity_0 + H1_0 + H123_0 + variance_10 + variance_0

# Generate LDA model on training data using statistics from stepwise model
model <- lda(formulaStepwise,data=train.data, method="moment")

# Predict selection types for test dataset using LDA model
predictions <- model %>% predict(test.data, CV=TRUE)

# Prediciton table
pt=as.data.table(table(pred=predictions$class, true=test.data$sim_type))

# Calculate accuracy proportions (N correct calls /50)
pt[pred=="Equilibrium Fluctuating" & true=="Equilibrium Fluctuating", N]/50
pt[pred=="Long-Term Fluctuating" & true=="Long-Term Fluctuating", N]/50
pt[pred=="Equilibrium Balancing" & true=="Equilibrium Balancing", N]/50
pt[pred=="Long-Term Balancing" & true=="Long-Term Balancing", N]/50
pt[pred=="Neutral" & true=="Neutral", N]/50
pt[, Freq:=N/50, by=c("pred", "true")]

overall_esw=mean(predictions$class==test.data$sim_type)

# Generate panel for confusion matrix figure
esw01=ggplot(pt, aes(pred,true, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#F8766D") +
  labs(title="Ends of seasons") + 
  theme(legend.position = "None", axis.title.x = element_blank(), axis.title.y = element_blank() )

## LDA using data from mid and end of summer season
train.data <- data[ n_win==250 & (gen_year==10 |gen_year==5)& run>50 ]
train.data=dcast(train.data, sel_type+sim_type+year+run+label ~ gen_year, value.var = stats)
test.data <- data[ n_win==250 & (gen_year==10 |gen_year==5)&run<=50]
test.data=test.data=dcast(test.data, sel_type+sim_type+year+run+label ~ gen_year, value.var = stats)

numberofCat=length(unique(train.data$sim_type))

formulaAll=sim_type~diversity_5+tajimas_d_site_5+theta_w_allele_5+ncd_5_5+ncd_4_5+ncd_3_5+ncd5_z_5+ncd4_z_5+ncd3_z_5+
  variance_5+skew_5+kurtosis_5+H1_5+H12_5+H123_5+H2H1_5+diversity_10+tajimas_d_site_10+theta_w_allele_10+ncd_5_10+ncd_4_10+
  ncd_3_10+ncd5_z_10+ncd4_z_10+ncd3_z_10+variance_10+skew_10+kurtosis_10+H1_10+H12_10+H123_10+H2H1_10

greedy.wilks(formulaAll,data=train.data, niveau = 0.05) 

formulaStepwise=sim_type ~ H2H1_10 + diversity_10 + kurtosis_10 + H1_10 + H123_10 + 
  H2H1_5 + tajimas_d_site_10 + ncd_3_10 + variance_10

model <- lda(formulaStepwise,data=train.data, method="moment")

predictions <- model %>% predict(test.data, CV=TRUE)
pt=as.data.table(table(pred=predictions$class, true=test.data$sim_type))

pt[pred=="Equilibrium Fluctuating" & true=="Equilibrium Fluctuating", N]/50
pt[pred=="Long-Term Fluctuating" & true=="Long-Term Fluctuating", N]/50
pt[pred=="Equilibrium Balancing" & true=="Equilibrium Balancing", N]/50
pt[pred=="Long-Term Balancing" & true=="Long-Term Balancing", N]/50
pt[pred=="Neutral" & true=="Neutral", N]/50
pt[, Freq:=N/50, by=c("pred", "true")]

overall_mes=mean(predictions$class==test.data$sim_type)

mes01=ggplot(pt, aes(pred,true, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#F8766D") +
  labs(title="Middle and end of summer") + 
  theme(legend.position = "None", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank() )


## LDA using data from central window, adjacent window and window 250kb away from selected site
train.data <- data[ (n_win==250|n_win==275|n_win==251) & gen_year==10& run>50]
train.data=dcast(train.data, sel_type+sim_type+year+run+label ~ n_win, value.var = stats)
test.data <- data[ (n_win==250|n_win==275|n_win==251) & gen_year==10 &run<=50]
test.data=test.data=dcast(test.data, sel_type+sim_type+year+run+label ~ n_win, value.var = stats)

numberofCat=length(unique(train.data$sim_type))

formulaAll=sim_type~diversity_250+tajimas_d_site_250+theta_w_allele_250+ncd_5_250+ncd_4_250+ncd_3_250+ncd5_z_250+ncd4_z_250+ncd3_z_250+
  variance_250+skew_250+kurtosis_250+H1_250+H12_250+H123_250+H2H1_250+diversity_275+tajimas_d_site_275+theta_w_allele_275+ncd_5_275+ncd_4_275+
  ncd_3_275+ncd5_z_275+ncd4_z_275+ncd3_z_275+variance_275+skew_275+kurtosis_275+H1_275+H12_275+H123_275+H2H1_275+diversity_251+
  tajimas_d_site_251+theta_w_allele_251+ncd_5_251+ncd_4_251+
  ncd_3_251+ncd5_z_251+ncd4_z_251+ncd3_z_251+variance_251+skew_251+kurtosis_251+H1_251+H12_251+H123_251+H2H1_251

greedy.wilks(formulaAll,data=train.data, niveau = 0.05) 

formulaStepwise=sim_type ~ H2H1_250 + diversity_250 + kurtosis_250 + H1_250 + 
  H123_250 + tajimas_d_site_250 + ncd_3_250 + ncd_4_251 + variance_250

model <- lda(formulaStepwise,data=train.data, method="moment")

predictions <- model %>% predict(test.data, CV=TRUE)
pt=as.data.table(table(pred=predictions$class, true=test.data$sim_type))

pt[pred=="Equilibrium Fluctuating" & true=="Equilibrium Fluctuating", N]/50
pt[pred=="Long-Term Fluctuating" & true=="Long-Term Fluctuating", N]/50
pt[pred=="Equilibrium Balancing" & true=="Equilibrium Balancing", N]/50
pt[pred=="Long-Term Balancing" & true=="Long-Term Balancing", N]/50
pt[pred=="Neutral" & true=="Neutral", N]/50
pt[, Freq:=N/50, by=c("pred", "true")]

overall_mw=mean(predictions$class==test.data$sim_type)

mw01=ggplot(pt, aes(pred,true, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#F8766D") +
  labs(title="Multiple windows") + 
  theme(legend.position = "None", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank() )

## Combine figure panels to create full confusion matrix (Fig S16)
s01=ggarrange(cw01, mw01, esw01, mes01, nrow=2,ncol=2, widths = c(1.25,1,1,1))
s01=annotate_figure(s01, bottom=text_grob("Prediction", size=14), left=text_grob("True", size=14, rot=90))
ggsave("plots/S16.jpg", s01, width = 13, height = 10)

##s=0.5

# Central window and end of summer
train.data <- data[ n_win==250 & gen_year==10& run>50]
test.data <- data[ n_win==250 & gen_year==10 &run<=50]

numberofCat=length(unique(train.data$sim_type))

formulaAll=sim_type~diversity+tajimas_d_site+theta_w_allele+ncd_5+ncd_4+ncd_3+ncd5_z+ncd4_z+ncd3_z+variance+skew+kurtosis+H1+H12+H123+H2H1

greedy.wilks(formulaAll,data=train.data, niveau = 0.05) 

formulaStepwise=sim_type ~ H2H1 + variance + theta_w_allele + H1 + H123 + ncd_3 + 
  diversity + tajimas_d_site + kurtosis + skew + H12

model <- lda(formulaStepwise,data=train.data, method="moment")

predictions <- model %>% predict(test.data, CV=TRUE)
pt=as.data.table(table(pred=predictions$class, true=test.data$sim_type))

pt[pred=="Equilibrium\nFluctuating" & true=="Equilibrium\nFluctuating", N]/50
pt[pred=="Long-Term\nFluctuating" & true=="Long-Term\nFluctuating", N]/50
pt[pred=="Equilibrium Balancing" & true=="Equilibrium Balancing", N]/50
pt[pred=="Long-Term Balancing" & true=="Long-Term Balancing", N]/50
pt[pred=="Neutral" & true=="Neutral", N]/50

overall_cw5=mean(predictions$class==test.data$sim_type)

pt[, Freq:=N/50, by=c("pred", "true")]

cw5=ggplot(pt, aes(pred,true, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +theme_bw()+
  scale_fill_gradient(low="white", high="#00BA38") +
  labs(title="Central window")+
  theme(legend.position = "None", axis.title.x = element_blank(),axis.title.y = element_blank())

stats=c("diversity","tajimas_d_site","theta_w_allele","ncd_5","ncd_4","ncd_3","ncd5_z","ncd4_z","ncd3_z",
        "variance","skew","kurtosis","H1","H12","H123","H2H1")


## end of seasons
train.data <- data[ n_win==250 & (gen_year==10 |gen_year==0)& run>50]
train.data=dcast(train.data, sel_type+sim_type+year+run+label ~ gen_year, value.var = stats)
test.data <- data[ n_win==250 & (gen_year==10 |gen_year==0)&run<=50]
test.data=test.data=dcast(test.data, sel_type+sim_type+year+run+label ~ gen_year, value.var = stats)

numberofCat=length(unique(train.data$sim_type))

# Stepwise model was colinear, hence NCD (TF=0.5) - winter (ncd_5_0) was remove from formulaAll
formulaAll=sim_type~diversity_0+tajimas_d_site_0+theta_w_allele_0+ncd_4_0+ncd_3_0+ncd5_z_0+ncd4_z_0+ncd3_z_0+
  variance_0+skew_0+kurtosis_0+H1_0+H12_0+H123_0+H2H1_0+diversity_10+tajimas_d_site_10+theta_w_allele_10+ncd_5_10+ncd_4_10+
  ncd_3_10+ncd5_z_10+ncd4_z_10+ncd3_z_10+variance_10+skew_10+kurtosis_10+H1_10+H12_10+H123_10+H2H1_10

greedy.wilks(formulaAll,data=train.data, niveau = 0.05) 

formulaStepwise=sim_type ~ H2H1_0 + H2H1_10 + variance_10 + variance_0 + theta_w_allele_0 + 
  H1_10 + H123_10 + tajimas_d_site_10 + tajimas_d_site_0 + 
  theta_w_allele_10 + ncd_3_0 + diversity_10 + H123_0 + H1_0 + 
  ncd_4_0 + kurtosis_10

model <- lda(formulaStepwise,data=train.data, method="moment")

predictions <- model %>% predict(test.data, CV=TRUE)
pt=as.data.table(table(pred=predictions$class, true=test.data$sim_type))

pt[pred=="Equilibrium\nFluctuating" & true=="Equilibrium\nFluctuating", N]/50
pt[pred=="Long-Term Fluctuating" & true=="Long-Term Fluctuating", N]/50
pt[pred=="Equilibrium Balancing" & true=="Equilibrium Balancing", N]/50
pt[pred=="Long-Term Balancing" & true=="Long-Term Balancing", N]/50
pt[pred=="Neutral" & true=="Neutral", N]/50
pt[, Freq:=N/50, by=c("pred", "true")]

overall_esw5=mean(predictions$class==test.data$sim_type)

esw5=ggplot(pt, aes(pred,true, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +theme_bw()+
  scale_fill_gradient(low="white", high="#00BA38") +
  labs(title="Ends of seasons") + 
  theme(legend.position = "None", axis.title.x = element_blank(), axis.title.y = element_blank())


## mid and end of season
train.data <- data[ n_win==250 & (gen_year==10 |gen_year==5)& run>50]
train.data=dcast(train.data, sel_type+sim_type+year+run+label ~ gen_year, value.var = stats)
test.data <- data[ n_win==250 & (gen_year==10 |gen_year==5)&run<=50 ]
test.data=test.data=dcast(test.data, sel_type+sim_type+year+run+label ~ gen_year, value.var = stats)

numberofCat=length(unique(train.data$sim_type))

formulaAll=sim_type~diversity_5+tajimas_d_site_5+theta_w_allele_5+ncd_5_5+ncd_4_5+ncd_3_5+ncd5_z_5+ncd4_z_5+ncd3_z_5+
  variance_5+skew_5+kurtosis_5+H1_5+H12_5+H123_5+H2H1_5+diversity_10+tajimas_d_site_10+theta_w_allele_10+ncd_5_10+ncd_4_10+
  ncd_3_10+ncd5_z_10+ncd4_z_10+ncd3_z_10+variance_10+skew_10+kurtosis_10+H1_10+H12_10+H123_10+H2H1_10

greedy.wilks(formulaAll,data=train.data, niveau = 0.05) 

formulaStepwise=sim_type ~ H2H1_10 + variance_10 + variance_5 + theta_w_allele_10 + 
  H1_10 + H123_10 + tajimas_d_site_10 + diversity_5 + diversity_10 + 
  H2H1_5 + ncd_3_10 + H123_5 + H1_5 + theta_w_allele_5 + ncd_5_5 + 
  tajimas_d_site_5 + kurtosis_10 + H12_10

model <- lda(formulaStepwise,data=train.data, method="moment")

predictions <- model %>% predict(test.data, CV=TRUE)
pt=as.data.table(table(pred=predictions$class, true=test.data$sim_type))

pt[pred=="Equilibrium\nFluctuating" & true=="Equilibrium\nFluctuating", N]/50
pt[pred=="Long-Term Fluctuating" & true=="Long-Term Fluctuating", N]/50
pt[pred=="Equilibrium Balancing" & true=="Equilibrium Balancing", N]/50
pt[pred=="Long-Term Balancing" & true=="Long-Term Balancing", N]/50
pt[pred=="Neutral" & true=="Neutral", N]/50

pt[, Freq:=N/50, by=c("pred", "true")]

overall_mes5=mean(predictions$class==test.data$sim_type)

mes5=ggplot(pt, aes(pred,true, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +theme_bw()+
  scale_fill_gradient(low="white", high="#00BA38") +
  labs(title="Middle and end of summer") + 
  theme(legend.position = "None", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank() )

## windows central, flanking and 250kb away
train.data <- data[ (n_win==250|n_win==275|n_win==251) & gen_year==10& run>50]
train.data=dcast(train.data, sel_type+sim_type+year+run+label ~ n_win, value.var = stats)
test.data <- data[ (n_win==250|n_win==275|n_win==251) & gen_year==10 &run<=50]
test.data=test.data=dcast(test.data, sel_type+sim_type+year+run+label ~ n_win, value.var = stats)

numberofCat=length(unique(train.data$sim_type))

formulaAll=sim_type~diversity_250+tajimas_d_site_250+theta_w_allele_250+ncd_5_250+ncd_4_250+ncd_3_250+ncd5_z_250+ncd4_z_250+ncd3_z_250+
  variance_250+skew_250+kurtosis_250+H1_250+H12_250+H123_250+H2H1_250+diversity_275+tajimas_d_site_275+theta_w_allele_275+ncd_5_275+ncd_4_275+
  ncd_3_275+ncd5_z_275+ncd4_z_275+ncd3_z_275+variance_275+skew_275+kurtosis_275+H1_275+H12_275+H123_275+H2H1_275+diversity_251+
  tajimas_d_site_251+theta_w_allele_251+ncd_5_251+ncd_4_251+
  ncd_3_251+ncd5_z_251+ncd4_z_251+ncd3_z_251+variance_251+skew_251+kurtosis_251+H1_251+H12_251+H123_251+H2H1_251

greedy.wilks(formulaAll,data=train.data, niveau = 0.05) 

formulaStepwise=sim_type ~ H2H1_250 + variance_250 + theta_w_allele_250 + H1_250 + 
  H123_250 + theta_w_allele_251 + ncd_3_250 + H2H1_251 + theta_w_allele_275 + 
  diversity_250 + tajimas_d_site_250 + skew_250 + kurtosis_250 + 
  H123_251 + H1_251 + H12_251
model <- lda(formulaStepwise,data=train.data, method="moment")

predictions <- model %>% predict(test.data, CV=TRUE)
pt=as.data.table(table(pred=predictions$class, true=test.data$sim_type))

pt[pred=="Equilibrium\nFluctuating" & true=="Equilibrium\nFluctuating", N]/50
pt[pred=="Long-Term\nFluctuating" & true=="Long-Term\nFluctuating", N]/50
pt[pred=="Equilibrium Balancing" & true=="Equilibrium Balancing", N]/50
pt[pred=="Long-Term Balancing" & true=="Long-Term Balancing", N]/50
pt[pred=="Neutral" & true=="Neutral", N]/50
pt[, Freq:=N/50, by=c("pred", "true")]

overall_mw5=mean(predictions$class==test.data$sim_type)

mw5=ggplot(pt, aes(pred,true, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) + theme_bw()+
  scale_fill_gradient(low="white", high="#00BA38") +
  labs(title="Multiple windows") + 
  theme(legend.position = "None", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank() )

# Generate Figure 7
s05=ggarrange(cw5, mw5, esw5, mes5, nrow=2, ncol=2, widths = c(1.25,1,1,1))
s05=annotate_figure(s05, bottom=text_grob("Prediction", size=14), left=text_grob("True", size=14, rot=90))
ggsave("plots/Fig7.jpg", s05, width = 13, height = 10)

##s=1

# Central window and end of summer
train.data <- data[ n_win==250 & gen_year==10& run>50]
test.data <- data[ n_win==250 & gen_year==10&run<=50]

numberofCat=length(unique(train.data$sim_type))

# Stepwise model was colinear, hence NCD (TF=0.5) (ncd_5) was remove from formulaAll
formulaAll=sim_type~diversity+tajimas_d_site+theta_w_allele+ncd_4+ncd_3+ncd5_z+ncd4_z+ncd3_z+variance+skew+kurtosis+H1+H2H1+ H12 + H123

greedy.wilks(formulaAll,data=train.data, niveau = 0.05) 

formulaStepwise=sim_type ~ H2H1 + H123 + theta_w_allele + H1 + variance + tajimas_d_site + 
  diversity + ncd_3 + ncd_4 + skew + kurtosis + ncd3_z + ncd4_z + 
  ncd5_z

model <- lda(formulaStepwise,data=train.data, method="moment")

predictions <- model %>% predict(test.data, CV=TRUE)
pt=as.data.table(table(pred=predictions$class, true=test.data$sim_type))

pt[pred=="Equilibrium Fluctuating" & true=="Equilibrium Fluctuating", N]/50
pt[pred=="Long-Term Fluctuating" & true=="Long-Term Fluctuating", N]/50
pt[pred=="Equilibrium Balancing" & true=="Equilibrium Balancing", N]/50
pt[pred=="Long-Term Balancing" & true=="Long-Term Balancing", N]/50
pt[pred=="Neutral" & true=="Neutral", N]/50

pt[, Freq:=N/50, by=c("pred", "true")]

overall_cw1=mean(predictions$class==test.data$sim_type)


cw=ggplot(pt, aes(pred,true, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#619CFF") +
  labs( title="Central window")+
  theme(legend.position = "None", axis.title.x = element_blank(), axis.title.y=element_blank())

stats=c("diversity","tajimas_d_site","theta_w_allele","ncd_5","ncd_4","ncd_3","ncd5_z","ncd4_z","ncd3_z",
        "variance","skew","kurtosis","H1","H12","H123","H2H1")


## end of seasons
train.data <- data[ n_win==250 & (gen_year==10 |gen_year==0)& run>50]
train.data=dcast(train.data, sel_type+sim_type+year+run+label ~ gen_year, value.var = stats)
test.data <- data[ n_win==250 & (gen_year==10 |gen_year==0)&run<=50]
test.data=test.data=dcast(test.data, sel_type+sim_type+year+run+label ~ gen_year, value.var = stats)

numberofCat=length(unique(train.data$sim_type))

# Stepwise model was colinear, hence NCD (TF=0.5) - winter/summer (ncd_5_0 and ncd_5_10) was remove from formulaAll
formulaAll=sim_type~diversity_0+tajimas_d_site_0+theta_w_allele_0+ncd_4_0+ncd_3_0+ncd5_z_0+ncd4_z_0+ncd3_z_0+
  variance_0+skew_0+kurtosis_0+H1_0+H12_0+H123_0+H2H1_0+diversity_10+tajimas_d_site_10+theta_w_allele_10+ncd_4_10+
  ncd_3_10+ncd5_z_10+ncd4_z_10+ncd3_z_10+variance_10+skew_10+kurtosis_10+H1_10+H12_10+H123_10+H2H1_10
greedy.wilks(formulaAll,data=train.data, niveau = 0.05) 

formulaStepwise=sim_type ~ H12_10 + H2H1_0 + theta_w_allele_10 + H2H1_10 + H123_10 + 
  theta_w_allele_0 + tajimas_d_site_10 + tajimas_d_site_0 + 
  skew_0 + diversity_10 + H1_10 + H1_0 + H123_0 + kurtosis_10 + 
  kurtosis_0 + ncd_4_0 + ncd_3_0 + ncd5_z_10 + ncd3_z_10 + 
  ncd_4_10 + ncd_3_10 + ncd4_z_10

model <- lda(formulaStepwise,data=train.data, method="moment")

predictions <- model %>% predict(test.data, CV=TRUE)
pt=as.data.table(table(pred=predictions$class, true=test.data$sim_type))

pt[pred=="Equilibrium Fluctuating" & true=="Equilibrium Fluctuating", N]/50
pt[pred=="Long-Term Fluctuating" & true=="Long-Term Fluctuating", N]/50
pt[pred=="Equilibrium Balancing" & true=="Equilibrium Balancing", N]/50
pt[pred=="Long-Term Balancing" & true=="Long-Term Balancing", N]/50
pt[pred=="Neutral" & true=="Neutral", N]/50
pt[, Freq:=N/50, by=c("pred", "true")]

overall_esw1=mean(predictions$class==test.data$sim_type)


esw=ggplot(pt, aes(pred,true, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#619CFF") +
  labs( title="End of seasons")+
  theme(legend.position = "None", axis.title.x = element_blank(), axis.title.y = element_blank())


## mid and end of season
train.data <- data[ n_win==250 & (gen_year==10 |gen_year==5)& run>50]
train.data=dcast(train.data, sel_type+sim_type+year+run+label ~ gen_year, value.var = stats)
test.data <- data[ n_win==250 & (gen_year==10 |gen_year==5)&run<=50]
test.data=test.data=dcast(test.data, sel_type+sim_type+year+run+label ~ gen_year, value.var = stats)

numberofCat=length(unique(train.data$sim_type))

# Stepwise model was colinear, hence NCD (TF=0.5) - mid and end summer (ncd_5_5 and ncd_5_10) was remove from formulaAll
formulaAll=sim_type~diversity_5+tajimas_d_site_5+theta_w_allele_5+ncd_4_5+ncd_3_5+ ncd5_z_5+ncd4_z_5+ncd3_z_5+ ncd5_z_10+ncd4_z_10+ncd3_z_10+
  variance_5+skew_5+kurtosis_5+H1_5+H12_5+H123_5+H2H1_5+diversity_10+tajimas_d_site_10+theta_w_allele_10+ncd_4_10+
  ncd_3_10+variance_10+skew_10+kurtosis_10+H1_10+H12_10+H123_10+H2H1_10

greedy.wilks(formulaAll,data=train.data, niveau = 0.05) 

formulaStepwise=sim_type ~ H2H1_5 + H123_10 + theta_w_allele_10 + H1_10 + variance_10 + 
  variance_5 + H2H1_10 + tajimas_d_site_10 + diversity_5 + 
  diversity_10 + H1_5 + H123_5 + ncd_3_10 + ncd_4_10 + ncd3_z_10 + 
  ncd4_z_10 + ncd5_z_10 + H12_5 + kurtosis_10 + kurtosis_5 + 
  tajimas_d_site_5 + ncd_3_5 + ncd_4_5 + theta_w_allele_5

model <- lda(formulaStepwise,data=train.data, method="moment")

predictions <- model %>% predict(test.data, CV=TRUE)
pt=as.data.table(table(pred=predictions$class, true=test.data$sim_type))

pt[pred=="Equilibrium Fluctuating" & true=="Equilibrium Fluctuating", N]/50
pt[pred=="Long-Term Fluctuating" & true=="Long-Term Fluctuating", N]/50
pt[pred=="Equilibrium Balancing" & true=="Equilibrium Balancing", N]/50
pt[pred=="Long-Term Balancing" & true=="Long-Term Balancing", N]/50
pt[pred=="Neutral" & true=="Neutral", N]/50
pt[, Freq:=N/50, by=c("pred", "true")]

overall_mes1=mean(predictions$class==test.data$sim_type)

mes=ggplot(pt, aes(pred,true, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#619CFF") +
  labs( title="Middle and end of summer")+
  theme(legend.position = "None", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())


## central window, flanking and 250kb away
train.data <- data[ (n_win==250|n_win==275|n_win==251) & gen_year==10& run>50]
train.data=dcast(train.data, sel_type+sim_type+year+run+label ~ n_win, value.var = stats)
test.data <- data[ (n_win==250|n_win==275|n_win==251) & gen_year==10 &run<=50]
test.data=test.data=dcast(test.data, sel_type+sim_type+year+run+label ~ n_win, value.var = stats)

numberofCat=length(unique(train.data$sim_type))

formulaAll=sim_type~diversity_250+tajimas_d_site_250+theta_w_allele_250+ncd_4_250+ncd_3_250+ncd5_z_250+ncd4_z_250+ncd3_z_250+
  variance_250+skew_250+kurtosis_250+H1_250+H12_250+H123_250+H2H1_250+diversity_275+tajimas_d_site_275+theta_w_allele_275+ncd_4_275+
  ncd_3_275+ncd5_z_275+ncd4_z_275+ncd3_z_275+variance_275+skew_275+kurtosis_275+H1_275+H12_275+H123_275+H2H1_275+diversity_251+
  tajimas_d_site_251+theta_w_allele_251+ncd_4_251+
  ncd_3_251+ncd5_z_251+ncd4_z_251+ncd3_z_251+variance_251+skew_251+kurtosis_251+H1_251+H12_251+H123_251+H2H1_251

greedy.wilks(formulaAll,data=train.data, niveau = 0.05) 

formulaStepwise=sim_type ~ H12_250 + H2H1_250 + theta_w_allele_250 + theta_w_allele_251 + 
  H123_250 + theta_w_allele_275 + variance_250 + H2H1_251 + 
  tajimas_d_site_250 + diversity_250 + ncd_3_250 + ncd_4_250 + 
  H1_250 + skew_250 + kurtosis_250 + H123_251 + H1_275 + ncd_4_251 + 
  ncd_3_251 + H1_251 + ncd_4_275 + ncd_3_275 + skew_251

model <- lda(formulaStepwise,data=train.data, method="moment")

predictions <- model %>% predict(test.data, CV=TRUE)
pt=as.data.table(table(pred=predictions$class, true=test.data$sim_type))

pt[pred=="Equilibrium Fluctuating" & true=="Equilibrium Fluctuating", N]/50
pt[pred=="Long-Term Fluctuating" & true=="Long-Term Fluctuating", N]/50
pt[pred=="Equilibrium Balancing" & true=="Equilibrium Balancing", N]/50
pt[pred=="Long-Term Balancing" & true=="Long-Term Balancing", N]/50
pt[pred=="Neutral" & true=="Neutral", N]/50
pt[, Freq:=N/50, by=c("pred", "true")]

overall_mw1=mean(predictions$class==test.data$sim_type)

mw=ggplot(pt, aes(pred,true, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#619CFF") +
  labs( title="Multiple windows")+
  theme(legend.position = "None", axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank() )

# Generate Figure S17
s1=ggarrange(cw, mw, esw, mes, nrow=2, ncol=2, widths = c(1.25,1,1,1))
s1=annotate_figure(s1, bottom=text_grob("Prediction", size=14), left=text_grob("True", size=14, rot=90))
ggsave("plots/FigS17.pdf", s1, width = 13, height = 10)




