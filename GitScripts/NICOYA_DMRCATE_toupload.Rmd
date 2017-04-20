---
title: "NICOYA_DMRCATE"
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

# DMRcate Analyses

```{r}
load("~/FinalObjects/meta_norep.RData")
load("~/FinalObjects/CRELES_notcelltypecorrected.RData")
mvals_uncorrected <- exprs(CRELES_notcelltypecorrected)
library(DMRcate)
nic.design <- model.matrix(~as.factor(nicoya)+ age+CD8.naive+CD8pCD28nCD45RAn+CD4T+Gran+NK+Mono+PlasmaBlast, data=meta_norep)
nic.myannotation <- cpg.annotate("array", object =mvals_uncorrected, analysis.type="differential", design=nic.design, coef=2) 
nic.dmrcoutput <- dmrcate(nic.myannotation, lambda=1000, C=2, pcutoff =0.05) #here, you can change the FDR and beta cutoffs
nic.results.ranges <- extractRanges(nic.dmrcoutput, genome = "hg19")
length(nic.results.ranges) #3328 @ p< 0.05
two.results <- nic.results.ranges[which(nic.results.ranges$no.cpgs>=3),] # looks like there are a bunch with 1 probe only... removing them 
length(two.results) #2666 with 3 or more probes 
two.results.df <- as.data.frame(two.results)
head(two.results.df)
sig_dmr_results <- subset(two.results.df, meanbetafc >= 0.05 | meanbetafc <= -0.05 & minfdr <= 0.05) #8
dim(sig_dmr_results) #20
```

