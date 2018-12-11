---
title: "NICOYA_EWAS_analyses_toupload"
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

# EWAS 

```{r}
load("~/FinalObjects/CRELES_notcelltypecorrected.RData")
load("~/FinalObjects/meta_norep.RData")
pdat <- meta_norep
mvals_uncor <- exprs(CRELES_notcelltypecorrected)

mod.nicoya <- model.matrix(~nicoya.x+ male.x + age+CD8.naive+CD8pCD28nCD45RAn+CD4T+Gran+NK+Mono+PlasmaBlast, data=pdat)
fit.nicoya <- lmFit(mvals_uncor, mod.nicoya_EWASA)  
eFit.nicoya <- eBayes(fit.nicoya_EWASA)
topt.nicoya <- topTable(eFit.nicoya_EWASA, coef = "nicoya.x1", n=Inf)
nicoya.hits <- subset(topt.nicoya_EWASA, P.Value <= 0.0000005) #9 

#calculate delta beta
betas_uncor <- betas(CRELES_notcelltypecorrected)
nicoya <- subset(pdat_correctedjune2016, nicoya =="1")
control <- subset(pdat_correctedjune2016, nicoya =="0")
delbeta_nicoya <- as.data.frame(rowMeans(betas_uncor[,rownames(nicoya)]) -rowMeans(betas_uncor[,rownames(control)]))
nicoya.hits$delB <- delbeta_nicoya[rownames(nicoya.hits),]

#table of sig hits
EWAS_resultsTable <- cbind(rownames(nicoya.hits), nicoya.hits[,c("P.Value","adj.P.Val", "delB")], fdat[rownames(nicoya.hits_EWAS),])
write.csv(EWAS_resultsTable, file = "~/FinalObjects/EWAS_resultsTable.csv")
```


