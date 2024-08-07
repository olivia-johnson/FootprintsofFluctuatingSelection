## Collate and analyse SFS output ##

if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, ggplot2, ggpubr)

## set path and working directory
path='/Path/to/FootprintsofFluctuatingSelection/'
setwd(path)

## forms of selection in sfs labels
types=c("hard", "soft", "neutral", "fluctuating", "250kb_fluctuating", "balancing")

## collate sfs data
alldata=NULL
for (sim_type in types){ # loop through sfs labels
  if (sim_type=="250kb_fluctuating"){ # set sel_type to extract files from correct directory
    sel_type="fluctuating"
  }else{
    sel_type=sim_type
    }
  # list of sfs files
  sfs_list  <- list.files(path =paste0(path, sel_type, "/"),pattern =paste0("sfs_", sim_type))
  # if sfs is only calculated at one timepoint (neutral/hard sweeps/soft sweep) 
  if (length(sfs_list)<51){
    # create vempty vector to add frequency spectrums to
    rep_folded=vector(mode="numeric", length=199)
    rep_unfolded=vector(mode="numeric", length=199)
    for (i in 1:50){
      # for each replciate add sfs at each frequency
      rep=fread(sfs_list[i])
      rep_folded=rep$Folded+rep_folded
      rep_unfolded=rep$Unfolded+rep_unfolded
    }
    # create data table with frequency as x and count and value, add label for folded or unfolded spectrum and type of seelction
    data_1=data.table(x=1:199, value=rep_folded, sfs="Folded", sim_type=sim_type) # folded dt
    data_2=data.table(x=1:199, value=rep_unfolded, sfs="Unfolded", sim_type=sim_type) #unfolded dt
    data=rbind(data_1, data_2)
    alldata=rbind(alldata, data) ## merge datasets to full dataset fo all replicates
  }
  else{ ## if sfs is calculated at two timepoints (balancing/fluctuating)
    # files from long-term timepoint
    lt_list = sfs_list[grep("_96050.txt", sfs_list)] 
    # files from early equilibrium timepoint
    eq_list = sfs_list[-grep("_96050.txt", sfs_list)]
    print(length(lt_list))
    print(length(eq_list))
    
    #create empty vectors
    eqrep_folded=vector(mode="numeric", length=199)
    eqrep_unfolded=vector(mode="numeric", length=199)
    for (i in 1:50){ # collate equilibrium sfs
      rep=fread(eq_list[i])
      eqrep_folded=rep$Folded+eqrep_folded
      eqrep_unfolded=rep$Unfolded+eqrep_unfolded
    }
    #create empty vectors
    ltrep_folded=vector(mode="numeric", length=199)
    ltrep_unfolded=vector(mode="numeric", length=199)
    for (i in 1:50){ # collate long-term sfs
      rep=fread(lt_list[i])
      ltrep_folded=rep$Folded+ltrep_folded
      ltrep_unfolded=rep$Unfolded+ltrep_unfolded
    }
    # create tables for each timepoitn and folded and unfolded spectrums
    eqdata_1=data.table(x=1:199, value=eqrep_folded, sfs="Folded", sim_type=paste0("eq_",sim_type))
    eqdata_2=data.table(x=1:199, value=eqrep_unfolded, sfs="Unfolded", sim_type=paste0("eq_",sim_type))
    ltdata_1=data.table(x=1:199, value=ltrep_folded, sfs="Folded", sim_type=paste0("lt_",sim_type))
    ltdata_2=data.table(x=1:199, value=ltrep_unfolded, sfs="Unfolded", sim_type=paste0("lt_",sim_type))
    data=rbind(eqdata_1, eqdata_2, ltdata_1, ltdata_2)
    alldata=rbind(alldata, data) # merge datasets
    
}
  
}
# Add labels to datasets
alldata[sim_type=="hard", label:="Hard Sweep"]
alldata[sim_type=="soft", label:="Soft Sweep"]
alldata[sim_type=="neutral", label:="Neutral"]
alldata[sim_type=="eq_wittmann_unlinked", label:="Early Eq.\nFluctuating\nCentral"]
alldata[sim_type=="lt_wittmann_unlinked", label:="Long-Term\nFluctuating\nCentral"]
alldata[sim_type=="eq_250kb_wittmann_unlinked", label:="Early Eq.\nFluctuating\n250kb"]
alldata[sim_type=="lt_250kb_wittmann_unlinked", label:="Long-Term\nFluctuating\n250kb"]
alldata[sim_type=="eq_balancing", label:="Early Eq.\nBalancing\n(s=0.1)"]
alldata[sim_type=="lt_balancing", label:="Long-Term\nBalancing\n(s=0.1)"]

# Create the foled and unfolded neutral expectations
folded_exp=data.table(x=1:100, y=(1/1:100)*alldata[sim_type=="neutral", max(value)]/20)
unfolded_exp=data.table(x=1:200, y=(1/1:200)*alldata[sim_type=="neutral", max(value)]/20)

# Plot of the folded spectrums
folded=ggplot(alldata[sfs=="Folded" & sim_type!="lt_250kb_wittmann_unlinked"], aes(x=x, y= value/20, fill=sim_type))+ 
  geom_line(data=folded_exp, aes(x=x, y=y), alpha=0.5, linewidth=0.3)+
  geom_col()+ 
  theme_bw()+
  facet_wrap(~factor(label, c("Neutral","Early Eq.\nFluctuating\nCentral","Early Eq.\nFluctuating\n250kb", "Long-Term\nFluctuating\nCentral",
                                 "Early Eq.\nBalancing\n(s=0.1)", "Long-Term\nBalancing\n(s=0.1)", "Hard Sweep", "Soft Sweep")), ncol=1, switch = "both") +
  theme(legend.position = 'none',axis.title.y = element_blank()) +labs(x="Folded", y="") +coord_cartesian(x=c(0,100))
  
# plot of the unfolded spectrums
unfolded=ggplot(alldata[sfs=="Unfolded" & sim_type!="lt_250kb_wittmann_unlinked"], aes(x=x, y= value/20, fill=sim_type))+ 
  geom_line(data=unfolded_exp, aes(x=x, y=y), alpha=0.5, linewidth=0.3)+
  geom_col()+ 
  theme_bw()+
  facet_wrap(~factor(label, c("Neutral","Early Eq.\nFluctuating\nCentral","Early Eq.\nFluctuating\n250kb", "Long-Term\nFluctuating\nCentral",
                                 "Early Eq.\nBalancing\n(s=0.1)", "Long-Term\nBalancing\n(s=0.1)", "Hard Sweep", "Soft Sweep")), ncol=1) +
  theme(legend.position = 'none',strip.text.x = element_blank(), axis.title.y = element_blank()) +labs(x="Unfolded", y="")

# add folded and unfolded panels
plot=ggarrange(folded, unfolded)

#export SFS plot (Fig S2)
ggsave("FigS2.pdf", plot=plot, width=8, height=9.5)
ggsave("FigS2.jpg", plot=plot, width=8, height=9.5)



