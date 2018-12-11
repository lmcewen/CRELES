---
title: "NICOYA_DNAm_Variability_analyses_toupload"
output: html_document
---

Please email: lmcewen@cmmt.ubc.ca if you have any questions about this code. 

### Publication details (accepted at Epigenetics and Chromatin on April 12, 2017):
__Title:__ Differential DNA methylation and lymphocyte proportions in a Costa Rican high longevity region

__Authors:__ Lisa M. McEwen, BSc; Alexander M. Morin, BSc.; Rachel D. Edgar, MSc.; Julia L. MacIsaac, PhD; Meaghan J. Jones, PhD.; William H. Dow, PhD.; Luis Rosero-Bixby, PhD.; Michael S. Kobor, PhD.; David H. Rehkopf, ScD., MPH.

__Date:__ April 2017


Overview
======
This R script is to preprocess and normalize data generated with the Infinium HumanMethylation450 Array (450K).

Whole blood was collected from 95 individuals from the [Costa Rican Berkeley CRELES cohort](http://www.creles.berkeley.edu/). We sampled from the Nicoyan population ('Nicoyans',n=47) and from other regions of Costa Rica('Controls' or 'non-Nicoyans',n=47). The Nicoya peninsula hosts a population of Costa Rica that exceed expected life expectancy, see [Luis Rosero-Bixby et al., 2013](http://pubmedcentralcanada.ca/pmcc/articles/PMC4241350/). The goal of this study was to explore the epigenetics (DNA methylation) of Nicoyan individuals as compared to non-Nicoyans, considering that DNA methylation is closely related to biological aging. 

Study Design: 
95 samples + 1 technical replicate, randomly distributed across 12 arrays ("chips", as each chip has 12 samples) run in one batch.

All R objects and files are located:
1. Kobor Server: big_data/lmcewen/Rehkopf
2. Stanford Server: [LINK](to be determined)

# DNAm Variability
```{r}
###FIGURE 3A 
library(lumi)
load("~/FinalObjects/CRELES_celltypecorrected.RData")
load("~/FinalObjects/meta_norep.RData")
CRELES_betas <- betas(CRELES_celltypecorrected)
nic <- subset(pdat_correctedjune2016, nicoya.x == "1")
con <- subset(pdat_correctedjune2016, nicoya.x == "0")

nicoya_betas <- as.data.frame(CRELES_betas[,rownames(nic)])
control_betas<- as.data.frame(CRELES_betas[,rownames(con)])

Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}

nicoya_IQR <-sapply(1:nrow(nicoya_betas), function(y)
    Variation(as.numeric(nicoya_betas[y,]))) 
control_IQR <-sapply(1:nrow(control_betas), function(y)
    Variation(as.numeric(control_betas[y,]))) 


IQR_nic_con <- data.frame("CpG "= rownames(CRELES_betas), 
                             "IQR_Nicoya" = nicoya_IQR,
                             "IQR_Control" = control_IQR)


IQR_nic_con$sd_negative<-IQR_nic_con$IQR_Nicoya-(0.20)*IQR_nic_con$IQR_Nicoya
IQR_nic_con$sd_positive<-IQR_nic_con$IQR_Nicoya+(0.20)*IQR_nic_con$IQR_Nicoya


IQR_nic_con$ind<-NA
for(i in 1:nrow(IQR_nic_con)){
  if(IQR_nic_con$IQR_Control[i]<IQR_nic_con$sd_negative[i]){ 
    IQR_nic_con$ind[i]<-"MOREVARNICOYA"
  } else if(IQR_nic_con$IQR_Control[i]>IQR_nic_con$sd_positive[i]){ 
    IQR_nic_con$ind[i]<-"MOREVARCONTROL"
  } else {
    IQR_nic_con$ind[i] <- "NA"
}}

#save(IQR_nic_con, file = "IQR_nic_con_2016Dec18.RData")
dim(IQR_nic_con[IQR_nic_con$ind == "MOREVARNICOYA",]) # 41695 CpGs are more variable in nicoya
dim(IQR_nic_con[IQR_nic_con$ind == "MOREVARCONTROL",]) #98073 CpGs are more variable in controls


###FIGURE 3B
fullcor <- betas(CRELES_fullcorrected)
nicoya <- subset(pdat_correctedjune2016, nicoya == "1")
nonnicoya <-  subset(pdat_correctedjune2016, nicoya == "0")
nic_betas <- fullcor[,rownames(nicoya)]
nonnic_betas <- fullcor[,rownames(nonnicoya)]

young_nic <- subset(nicoya, age.x <= 80)
old_nic <- subset(nicoya, age.x >= 80)
young_nic_betas <- nic_betas[,rownames(young_nic)]
old_nic_betas <- nic_betas[,rownames(old_nic)]

young_nonnic <- subset(nonnicoya, age.x <= 80)
old_nonnic <- subset(nonnicoya, age.x >= 80)
young_nonnic_betas <- nonnic_betas[,rownames(young_nonnic)]
old_nonnic_betas <- nonnic_betas[,rownames(old_nonnic)]

Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}
young_nic_betas_IQR <-sapply(1:nrow(young_nic_betas), function(y)
  Variation(as.numeric(young_nic_betas[y,]))) 
old_nic_betas_IQR <-sapply(1:nrow(old_nic_betas), function(y)
  Variation(as.numeric(old_nic_betas[y,]))) 
young_nonnic_betas_IQR <-sapply(1:nrow(young_nonnic_betas), function(y)
  Variation(as.numeric(young_nonnic_betas[y,]))) 
old_nonnic_betas_IQR <-sapply(1:nrow(old_nonnic_betas), function(y)
  Variation(as.numeric(old_nonnic_betas[y,]))) 

#save(young_nic_betas_IQR, young_nonnic_betas_IQR, old_nic_betas_IQR, old_nonnic_betas_IQR, file = "nicoya_nonnicoya_IQR_2016Nov23.RData")

young_nicoya_top <- young_nic_betas_IQR[young_nic_betas_IQR >= 0.05] #115,552
old_nicoya_top <- old_nic_betas_IQR[old_nic_betas_IQR >= 0.05] #125,327
young_nonnic_top <- young_nonnic_betas_IQR[young_nonnic_betas_IQR >= 0.05] #120,017
old_nonnic_top <- old_nonnic_betas_IQR[old_nonnic_betas_IQR >= 0.05] #141,922
