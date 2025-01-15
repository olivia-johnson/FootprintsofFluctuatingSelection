## Comparing varying strengths of fluctuating selection

## Load packages and install any missing from the user's R library
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, ggplot2, egg, ggpubr, geomtextpath, patchwork)

## Set path to directory
path='/Path/to/FootprintsofFluctuatingSelection/'
setwd(path)

## Load in data for neutral case to standardise NCD and generate neutral quantiles.
load(file=paste0(path, "single_locus_neutral.RData"))
sum_stats=sum_stats[label=="CP10k_EG_No Fitness"] 
ggplot(sum_stats)+ geom_density(aes(x=n_seg_sites, y=..density..)) # check distribution of neutral sites
bins=seq(from=0,to= 25*((max(sum_stats$n_seg_sites)%/%25)+1), by=25) # create bins of number of segregating sites
sum_stats[, bin:=cut(n_seg_sites, bins)] # bin windows  by number of segregating sites
#calculate mean and sd to use to standardise NCD based on number of segregating sites per window, as suggested in Bitarello et al. 2018.
ncd_stand=sum_stats[, .(m.ncd3=mean(ncd_3), m.ncd4=mean(ncd_4), m.ncd5=mean(ncd_5), 
                        sd.ncd3=sd(ncd_3), sd.ncd4=sd(ncd_4), sd.ncd5=sd(ncd_5)), by="bin"] 

sum_stats=merge(sum_stats, ncd_stand, by="bin") # add standardising values to summary data table

sum_stats[,ncd5_z:=(ncd_5-m.ncd5)/sd.ncd5, by=c("bin", "ncd_5")] # standardise NCD (TF=0.5)
sum_stats[,ncd4_z:=(ncd_4-m.ncd4)/sd.ncd4, by=c("bin", "ncd_4")] # standardise NCD (TF=0.4)
sum_stats[,ncd3_z:=(ncd_3-m.ncd3)/sd.ncd3, by=c("bin", "ncd_3")] # standardise NCD (TF=0.3)

## caculate 5%, 50% and 95% quartiles for neutral simulations to compare against fluctuating selection
div_qt=sum_stats[Gen==2560,quantile(diversity, probs=c(0.05, 0.5,0.95))]
thetaw_qt=sum_stats[Gen==2560,quantile(theta_w_allele, probs=c(0.05, 0.5,0.95))]
taj_qt=sum_stats[Gen==2560,quantile(tajimas_d_branch, probs=c(0.05, 0.5,0.95))]
ncd_qt=sum_stats[Gen==2560,quantile(ncd_5, probs=c(0.05, 0.5,0.95))]
ncdz_qt=sum_stats[Gen==2560,quantile(ncd5_z, probs=c(0.05, 0.5,0.95))]
var_qt=sum_stats[Gen==2560,quantile(variance, probs=c(0.05, 0.5,0.95))]
skew_qt=sum_stats[Gen==2560,quantile(skew, probs=c(0.05, 0.5,0.95))]
kurt_qt=sum_stats[Gen==2560,quantile(kurtosis, probs=c(0.05, 0.5,0.95))]
h1_qt=sum_stats[Gen==2560,quantile(H1, probs=c(0.05, 0.5,0.95))]
h12_qt=sum_stats[Gen==2560,quantile(H12, probs=c(0.05, 0.5,0.95))]
h123_qt=sum_stats[Gen==2560,quantile(H123, probs=c(0.05, 0.5,0.95))]
h2h1_qt=sum_stats[Gen==2560,quantile(H2H1, probs=c(0.05, 0.5,0.95))]


groups=c(1:6) # vector of unique parameter set IDs
load(file=paste0(path,"single_locus_fluctuating.RData"))
## Standardise NCD
bins=seq(from=0,to= 25*((max(sum_stats$n_seg_sites)%/%25)+1), by=25) # create bins
sum_stats[, bin:=cut(n_seg_sites, bins)] # bin windows by number of segregating sites
sum_stats=merge(sum_stats, ncd_stand, by="bin") # merge to dataset

sum_stats[,ncd5_z:=(ncd_5-m.ncd5)/sd.ncd5, by=c("bin", "ncd_5")] # Standardise NCD (TF = 0.5)
sum_stats[,ncd4_z:=(ncd_4-m.ncd4)/sd.ncd4, by=c("bin", "ncd_4")] # Standardise NCD (TF = 0.4)
sum_stats[,ncd3_z:=(ncd_3-m.ncd3)/sd.ncd3, by=c("bin", "ncd_3")] # Standardise NCD (TF = 0.3)

sum_stats[, year:=ifelse(gen_season=="EG",Gen%/%20, Gen%/%12), by="Gen"] # Determine the year/seasonal cycle of each generation


load(file=paste0(path,"af_single_locus_fluctuating.RData"))
freq_data[, year:=ifelse(gen_season=="EG",Gen%/%20, Gen%/%12), by="Gen"] # Determine the year/seasonal cycle of each generation

freq_data[, s_label:=paste0("s = ", s_s), by="s_s"] # Label for selection coefficient
freq_data[, h_label:=paste0("h = ", h_s), by="h_s"]  # Label for dominance coefficient
start_gens=c(cumprod(c(20, rep.int(2, 9))), 12000, 24000, 33000,48000,60000,72000,46000,960000) # Vector of standard times of sampling in full simulations

wittmann=sum_stats[Gen >=96040 & Gen<96060] # subset long-term time point
wittmann[, sim_type:="Long-Term Fluctuating"] # add label
wittmanne=NULL # subset data for early equilibrium timepoint, for loop that conditions of time equilibrium is reached for each parameter set and replicate
for (g in groups){
e.year=freq_data[group==g & gen_year==5|gen_year==15, mean(mut_freq), by=c("year", "run")][V1>0.45, min(year), by="run"]
for (j in 1:50){
  eq=e.year[ run==j, V1]
  values=sort(unique(sum_stats[group==g & run==j]$year))
  if (eq %in% values){
    sum_wit=sum_stats[group==g & run==j & year==eq &edges==F ]
    if (length(sum_wit[,unique(Gen)])!=20){
      val=values[values>eq]
      nearest=min(val)
      sum_wit=sum_stats[group==g & run==j & year==nearest &edges==F ]}}else{
        val=values[values>eq]
        nearest=min(val)
        sum_wit=sum_stats[group==g & run==j & year==nearest &edges==F ]}
  wittmanne=rbind(sum_wit, wittmanne)}
if (sum(wittmanne[group==g,length(unique(Gen)), by="run"]$V1==20)==50){print("Full season")}else{print("Not full season! - Wittmann")}}
wittmanne[, sim_type:="Early Equilibrium Fluctuating"]

wittmann=rbind(wittmann, wittmanne) # merge early equilibrium and long-term data

data=wittmann[edges==F & (gen_year==0 | gen_year==10|gen_year==5 | gen_year==15)] # remove windows that are not full or are over linkage breaks and only mid and ends of season
data[gen_year==0 | gen_year==10, fluctuation:=ifelse(gen_year==10, "End Summer", "End Winter"), by="gen_year"] # label time in season
data[gen_year==5 | gen_year==15, fluctuation:=ifelse(gen_year==5, "Mid Summer", "Mid Winter"), by="gen_year"]# label time in season
## calculate mean fo each statistic for each window, parameter set and seasonal timepoint
data[, mean_tajd:=mean(tajimas_d_branch), by=c( "midpoint", "sim_type","fluctuation", "label")]
data[, mean_div:=mean(diversity), by=c( "midpoint",  "sim_type","fluctuation", "label")]
data[, mean_thetaw:=mean(theta_w_allele), by=c( "midpoint",  "sim_type","fluctuation", "label")]
data[, mean_H1:=mean(H1), by=c( "midpoint",  "sim_type","fluctuation", "label")]
data[, mean_H12:=mean(H12), by=c( "midpoint",  "sim_type","fluctuation", "label")]
data[, mean_H123:=mean(H123), by=c( "midpoint",  "sim_type","fluctuation", "label")]
data[, mean_H2H1:=mean(H2H1), by=c( "midpoint",  "sim_type","fluctuation", "label")]
data[, mean_ncd_5:=mean(ncd_5), by=c( "midpoint",  "sim_type","fluctuation", "label")]
data[, mean_ncd5_z:=mean(ncd5_z), by=c( "midpoint",  "sim_type","fluctuation", "label")]
data[, mean_ncd_4:=mean(ncd_4), by=c( "midpoint",  "sim_type","fluctuation", "label")]
data[, mean_ncd_3:=mean(ncd_3), by=c( "midpoint",  "sim_type","fluctuation", "label")]
data[, mean_var:=mean(variance), by=c( "midpoint",  "sim_type","fluctuation", "label")]
data[, mean_skew:=mean(skew), by=c( "midpoint",  "sim_type","fluctuation", "label")]
data[, mean_kurtosis:=mean(kurtosis), by=c( "midpoint",  "sim_type","fluctuation", "label")]

data[, s_label:=paste0("s = ", s_s), by="s_s"] #label for seelction coefficient
data[, h_label:=paste0("h = ", h_s), by="h_s"] #label for dominance coefficient


## Plots ##

## Figure 2 - Footprints of fluctuating selection across time using SFS-based statistics.

div=ggplot(data[h_s==0.6&gen_year==10], aes(x=dist/100,  col=factor(s_s)))+
  geom_vline(data = subset(data),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=mean_div, group=block), alpha=0.8, size=.5)+
  theme_bw()+ 
  geom_texthline(aes(yintercept=div_qt[[1]], label = "5%"), lty = 2, colour= 'black',hjust = 0.95)+
  geom_texthline(aes(yintercept=div_qt[[2]], label = "50%"), lty = 2, colour= 'black',hjust = 0.95)+
  geom_texthline(aes(yintercept=div_qt[[3]], label = "95%"), lty = 2, colour= 'black',hjust = 0.95)+
  labs(x="Distance from selected site (Mb)", y="Nucleotide Diversity", col="Selection Coefficient")+
  facet_wrap(~sim_type, ncol=4)+
  theme(legend.position = "None")+
  theme(axis.title.x = element_blank())+
  coord_cartesian(x=c(-1,1))

theta=ggplot(data[h_s==0.6&gen_year==10], aes(x=dist/100,  col=factor(s_s)))+  
  geom_vline(data = subset(data),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=mean_thetaw, group=block), alpha=0.8, size=.5)+
  theme_bw()+
  geom_texthline(aes(yintercept=thetaw_qt[[1]], label = "5%"), lty = 2, colour= 'black',hjust = 0.95)+
  geom_texthline(aes(yintercept=thetaw_qt[[2]], label = "50%"), lty = 2, colour= 'black',hjust = 0.95)+
  geom_texthline(aes(yintercept=thetaw_qt[[3]], label = "95%"), lty = 2, colour= 'black',hjust = 0.95)+
  labs(x="Distance from selected site (Mb)", y="Watterson's Theta", col="Selection Coefficient")+
  geom_vline(data = subset(data),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  facet_wrap(~sim_type, ncol=4)+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "None")+
  coord_cartesian(x=c(-1,1))

taj=ggplot(data[h_s==0.6&gen_year==10], aes(x=dist/100,  col=factor(s_s)))+
  geom_vline(data = subset(data),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=mean_tajd, group=block), alpha=0.8, size=.5)+
  theme_bw()+
  geom_texthline(aes(yintercept=taj_qt[[1]], label = "5%"), lty = 2, colour= 'black',hjust = 0.95)+
  geom_texthline(aes(yintercept=taj_qt[[2]], label = "50%"), lty = 2, colour= 'black',hjust = 0.95)+
  geom_texthline(aes(yintercept=taj_qt[[3]], label = "95%"), lty = 2, colour= 'black',hjust = 0.95)+
  labs(x="Distance from selected site (Mb)", y="Tajima's D", col="Selection Coefficient")+
  facet_wrap(~sim_type, ncol=4)+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "None")+
  coord_cartesian(x=c(-1,1))

ncd5=ggplot(data[h_s==0.6&gen_year==10], aes(x=dist/100,  col=factor(s_s)))+
  geom_vline(data = subset(data),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=mean_ncd_5, group=block), alpha=0.8, size=.5)+
  theme_bw()+
  geom_texthline(aes(yintercept=ncd_qt[[1]], label = "5%"), lty = 2, colour= 'black',hjust = 0.95)+
  geom_texthline(aes(yintercept=ncd_qt[[2]], label = "50%"), lty = 2, colour= 'black',hjust = 0.95)+
  geom_texthline(aes(yintercept=ncd_qt[[3]], label = "95%"), lty = 2, colour= 'black',hjust = 0.95)+
  labs(x="Distance from selected site (Mb)", y="NCD (TF = 0.5)", col="Selection Coefficient")+
  facet_wrap(~sim_type, ncol=4)+ theme(axis.title.x.bottom = element_blank())+
  theme(legend.position="none")+
  coord_cartesian(x=c(-1,1))

ncd5z=ggplot(data[h_s==0.6&gen_year==10], aes(x=dist/100,  col=factor(s_s)))+
  geom_vline(data = subset(data),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=mean_ncd5_z, group=block), alpha=0.8, size=.5)+
  theme_bw()+
  geom_texthline(aes(yintercept=ncdz_qt[[1]], label = "5%"), lty = 2, colour= 'black',hjust = 0.95)+
  geom_texthline(aes(yintercept=ncdz_qt[[2]], label = "50%"), lty = 2, colour= 'black',hjust = 0.95)+
  geom_texthline(aes(yintercept=ncdz_qt[[3]], label = "95%"), lty = 2, colour= 'black',hjust = 0.95)+
  labs(x="Distance from selected site (Mb)", y="Standardised NCD (TF = 0.5)", col="Selection Coefficient")+
  facet_wrap(~sim_type, ncol=4)+ theme(axis.title.x.bottom = element_blank())+
  coord_cartesian(x=c(-1,1))

var=ggplot(data[h_s==0.6&gen_year==10], aes(x=dist/100,  col=factor(s_s)))+
  geom_vline(data = subset(data),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=mean_var, group=block), alpha=0.8, size=.5)+
  theme_bw()+
  geom_texthline(aes(yintercept=var_qt[[1]], label = "5%"), lty = 2, colour= 'black',hjust = 0.95)+
  geom_texthline(aes(yintercept=var_qt[[2]], label = "50%"), lty = 2, colour= 'black',hjust = 0.95)+
  geom_texthline(aes(yintercept=var_qt[[3]], label = "95%"), lty = 2, colour= 'black',hjust = 0.95)+
  labs(x="Distance from selected site (Mb)", y="Variance", col="Selection Coefficient")+
  facet_wrap(~sim_type, ncol=4)+
  theme(axis.title.x = element_blank())+
  coord_cartesian(x=c(-1,1))

skew=ggplot(data[h_s==0.6&gen_year==10], aes(x=dist/100,  col=factor(s_s)))+
  geom_vline(data = subset(data),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=mean_skew, group=block), alpha=0.8, size=.5)+
  theme_bw()+
  geom_texthline(aes(yintercept=skew_qt[[1]], label = "5%"), lty = 2, colour= 'black',hjust = 0.95)+
  geom_texthline(aes(yintercept=skew_qt[[2]], label = "50%"), lty = 2, colour= 'black',hjust = 0.95)+
  geom_texthline(aes(yintercept=skew_qt[[3]], label = "95%"), lty = 2, colour= 'black',hjust = 0.95)+
  labs(x="Distance from selected site (Mb)", y="Skew", col="Selection Coefficient")+
  facet_wrap(~sim_type, ncol=4)+
  theme(axis.title.x = element_blank())+
  coord_cartesian(x=c(-1,1))

kurtosis=ggplot(data[h_s==0.6&gen_year==10], aes(x=dist/100,  col=factor(s_s)))+
  geom_vline(data = subset(data),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=mean_kurtosis, group=block), alpha=0.8, size=.5)+
  theme_bw()+
  geom_texthline(aes(yintercept=kurt_qt[[1]], label = "5%"), lty = 2, colour= 'black',hjust = 0.95)+
  geom_texthline(aes(yintercept=kurt_qt[[2]], label = "50%"), lty = 2, colour= 'black',hjust = 0.95)+
  geom_texthline(aes(yintercept=kurt_qt[[3]], label = "95%"), lty = 2, colour= 'black',hjust = 0.95)+
  labs(x="Distance from selected site (Mb)", y="Kurtosis", col="Selection Coefficient")+
  facet_wrap(~sim_type, ncol=4)+
  theme(axis.title.x = element_blank())+
  coord_cartesian(x=c(-1,1))

fig2=(div/theta/taj/ncd5)|(var/skew/kurtosis/ncd5z)+ plot_layout(guides = "collect")&theme(legend.position="bottom")
fig2=wrap_elements(panel = fig2) +
  labs(tag = "Distance from selected site (Mb)") +
  theme(plot.tag = element_text(size = rel(1.25)),
        plot.tag.position = "bottom") 
ggexport(fig2, filename="plots/Fig2.pdf", width=10, height=11)
ggsave(fig2, filename="plots/Fig2.jpg", width=10, height=11)

## Figure 3 - Signatures of fluctuating selection in haplotype statistics.

H1_eq=ggplot(data[h_s==0.6& gen_year==10 &sim_type=="Early Equilibrium Fluctuating"], aes(x=dist/100,  col=factor(s_s)))+
  geom_vline(data = subset(data[sim_type=="Early Equilibrium Fluctuating"]),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=mean_H1, group=block), alpha=0.8, size=.5)+
  theme_bw()+
  labs(x="Distance from selected site (Mb)", y="Garud's H1", col="Selection Coefficient")+
  facet_wrap(~sim_type, ncol=4, scale="free_y")+
  theme(axis.title.x.bottom = element_blank(), legend.position = "None")+
  coord_cartesian(x=c(-1,1), y=c(0.005,0.145))
H1_lt=ggplot(data[h_s==0.6& gen_year==10&sim_type=="Long-Term Fluctuating"], aes(x=dist/100,  col=factor(s_s)))+
  geom_vline(data = subset(data[sim_type=="Long-Term Fluctuating"]),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=mean_H1, group=block), alpha=0.8, size=.5)+
  theme_bw()+
  labs(x="Distance from selected site (Mb)", y="Garud's H1", col="Selection Coefficient")+
  facet_wrap(~sim_type, ncol=4, scale="free_y")+
  theme(axis.title.x.bottom = element_blank(), axis.title.y = element_blank())+
  coord_cartesian(x=c(-1,1), y=c(0.007,0.014))
H1=ggarrange(H1_eq, H1_lt, nrow=1, widths=c(1,1.35))

H12_eq=ggplot(data[h_s==0.6& gen_year==10 &sim_type=="Early Equilibrium Fluctuating"], aes(x=dist/100,  col=factor(s_s)))+
  geom_vline(data = subset(data[sim_type=="Early Equilibrium Fluctuating"]),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=mean_H12, group=block), alpha=0.8, size=.5)+
  theme_bw()+
  labs(x="Distance from selected site (Mb)", y="Garud's H12", col="Selection Coefficient")+
  facet_wrap(~sim_type, ncol=4, scale="free_y")+
  theme(axis.title.x.bottom = element_blank(), legend.position = "None")+
  coord_cartesian(x=c(-1,1), y=c(0.005,0.145))
H12_lt=ggplot(data[h_s==0.6& gen_year==10&sim_type=="Long-Term Fluctuating"], aes(x=dist/100,  col=factor(s_s)))+
  geom_vline(data = subset(data[sim_type=="Long-Term Fluctuating"]),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=mean_H12, group=block), alpha=0.8, size=.5)+
  theme_bw()+
  labs(x="Distance from selected site (Mb)", y="Garud's H12", col="Selection Coefficient")+
  facet_wrap(~sim_type, ncol=4, scale="free_y")+
  theme(axis.title.x.bottom = element_blank(), axis.title.y = element_blank())+
  coord_cartesian(x=c(-1,1), y=c(0.007,0.014))
H12=ggarrange(H12_eq, H12_lt, nrow=1, widths=c(1,1.35))

H123_eq=ggplot(data[h_s==0.6& gen_year==10 &sim_type=="Early Equilibrium Fluctuating"], aes(x=dist/100,  col=factor(s_s)))+
  geom_vline(data = subset(data[sim_type=="Early Equilibrium Fluctuating"]),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=mean_H123, group=block), alpha=0.8, size=.5)+
  theme_bw()+
  labs(x="Distance from selected site (Mb)", y="Garud's H123", col="Selection Coefficient")+
  facet_wrap(~sim_type, ncol=4, scale="free_y")+
  theme(axis.title.x.bottom = element_blank(), legend.position = "None")+
  coord_cartesian(x=c(-1,1), y=c(0.005,0.145))
H123_lt=ggplot(data[h_s==0.6& gen_year==10&sim_type=="Long-Term Fluctuating"], aes(x=dist/100,  col=factor(s_s)))+
  geom_vline(data = subset(data[sim_type=="Long-Term Fluctuating"]),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=mean_H123, group=block), alpha=0.8, size=.5)+
  theme_bw()+
  labs(x="Distance from selected site (Mb)", y="Garud's H123", col="Selection Coefficient")+
  facet_wrap(~sim_type, ncol=4, scale="free_y")+
  theme(axis.title.x.bottom = element_blank(), axis.title.y = element_blank())+
  coord_cartesian(x=c(-1,1), y=c(0.007,0.014))
H123=ggarrange(H123_eq, H123_lt, nrow=1, widths=c(1,1.35))

H2H1_eq=ggplot(data[h_s==0.6& gen_year==10 &sim_type=="Early Equilibrium Fluctuating"], aes(x=dist/100,  col=factor(s_s)))+
  geom_vline(data = subset(data[sim_type=="Early Equilibrium Fluctuating"]),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=mean_H2H1, group=block),alpha=0.8, size=.5)+
  theme_bw()+
  labs(x="Distance from selected site (Mb)", y="Garud's H2/H1", col="Selection Coefficient")+
  facet_wrap(~sim_type, ncol=4, scale="free_y")+
  theme(axis.title.x.bottom = element_blank(), legend.position = "None")+
  coord_cartesian(x=c(-1,1), y=c(0.15,0.95))
H2H1_lt=ggplot(data[h_s==0.6& gen_year==10&sim_type=="Long-Term Fluctuating"], aes(x=dist/100,  col=factor(s_s)))+
  geom_vline(data = subset(data[sim_type=="Long-Term Fluctuating"]),aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=mean_H2H1, group=block), alpha=0.8, size=.5)+
  theme_bw()+
  labs(x="Distance from selected site (Mb)", y="Garud's H2/H1", col="Selection Coefficient")+
  facet_wrap(~sim_type, ncol=4, scale="free_y")+
  theme(axis.title.x.bottom = element_blank(), axis.title.y = element_blank())+
  coord_cartesian(x=c(-1,1), y=c(0.89,0.955))
H2H1=ggarrange(H2H1_eq, H2H1_lt, nrow=1, widths=c(1,1.35))

fig3=(H1_eq|H1_lt)/(H12_eq|H12_lt)/(H123_eq|H123_lt)/(H2H1_eq|H2H1_lt)+plot_layout(guides = "collect")&theme(legend.position="bottom", legend.justification = c(1,0),legend.margin = margin(0, 0, 0, 0))
fig3=wrap_elements(panel = fig3) +
  labs(tag = "Distance from selected site (Mb)") +
  theme(plot.tag = element_text(size = rel(1.25)),
        plot.tag.position = "bottom")
ggexport(fig3, filename="plots/Fig3.pdf", width=8, height=10)
ggsave(fig3, filename="plots/Fig3.jpg", width=8, height=10)

## Figure 4 - Signatures of fluctuating selection differ at different points of the seasonal cycle.

fH_eq = ggplot(data[h_s==0.6 & s_s==1 & sim_type=="Early Equilibrium Fluctuating"], aes(x=dist/100))+
  geom_vline(aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=H1,group=block, col=factor(s_s)), alpha=0.3, col="#619CFF")+
  geom_line(aes(y=mean_H1), col="black", size=0.8)+
  theme_bw()+
  labs(x="Distance from selected site (Mb)", y="Garud's H1")+
  facet_wrap(~sim_type+fluctuation, ncol=4)+
  theme( legend.position ="none", axis.title.x = element_blank(), axis.title.y = element_blank())+
  coord_cartesian(x=c(-0.25,0.25), y=c(0.0065, 0.2))

fH_lt = ggplot(data[h_s==0.6 & s_s==1 & sim_type=="Long-Term Fluctuating"], aes(x=dist/100))+
  geom_vline(aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=H1,group=block, col=factor(s_s)), alpha=0.3, col="#619CFF")+
  geom_line(aes(y=mean_H1), col="black", size=0.8)+
  theme_bw()+
  labs(x="Distance from selected site (Mb)", y="Garud's H1")+
  facet_wrap(~sim_type+fluctuation, ncol=4)+
  theme( legend.position ="none", axis.title.x = element_blank(), axis.title.y = element_blank())+
  coord_cartesian(x=c(-0.25,0.25), y=c(0.0065, 0.0135))
fH=ggarrange(fH_eq,fH_lt, ncol=1)
fH=annotate_figure(fH, left ="Garud's H1")

ft_eq = ggplot(data[h_s==0.6 & s_s==1 & sim_type=="Early Equilibrium Fluctuating"], aes(x=dist/100))+
  geom_vline(aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=diversity,group=block, col=factor(s_s)), alpha=0.3, col="#619CFF")+
  geom_line(aes(y=mean_div), col="black", size=0.8)+
  theme_bw()+
  labs(x="Distance from selected site (Mb)", y="Garud's H1")+
  facet_wrap(~sim_type+fluctuation, ncol=4)+
  theme( legend.position ="none", axis.title.x = element_blank(), axis.title.y = element_blank())+
  coord_cartesian(x=c(-0.25,0.25))

ft_lt = ggplot(data[h_s==0.6 & s_s==1 & sim_type=="Long-Term Fluctuating"], aes(x=dist/100))+
  geom_vline(aes(xintercept = 0), col= "red",linetype="dotted", size=.7)+
  geom_line(aes(y=diversity,group=block, col=factor(s_s)), alpha=0.3, col="#619CFF")+
  geom_line(aes(y=mean_div), col="black", size=0.8)+
  theme_bw()+
  labs(x="Distance from selected site (Mb)", y="Garud's H1")+
  facet_wrap(~sim_type+fluctuation, ncol=4)+
  theme( legend.position ="none", axis.title.x = element_blank(), axis.title.y = element_blank())+
  coord_cartesian(x=c(-0.7,0.7))
ft=ggarrange(ft_eq,ft_lt, ncol=1)
ft=annotate_figure(ft, left ="Nucleotide Diversity")
fig4=ggarrange(fH,ft, ncol=1)
fig4=annotate_figure(fig4, bottom = "Distance from selected site (Mb)")
ggexport(fig4, filename="plots/Fig4.pdf", width=8, height=10)
ggsave(plot=fig4, filename="plots/Fig4.jpg", width=8, height=10)

