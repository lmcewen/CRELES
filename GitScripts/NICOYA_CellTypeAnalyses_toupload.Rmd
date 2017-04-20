---
title: "Nicoya_CellTypeAnalyses"
output: html_document
---


# CRELES Study

### Publication details (currently submitted to Epigenetics and Chromatin):
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
2. Stanford Server: [LINK](tbd)

# Estimated cell type analyses with chronological age and nicoya group

```{r eval = FALSE}
load("~/FinalObjects/meta_norep.RData")

#hannum estimated cell proportions (fig 1a)
kruskal.test(meta_norep$CD8T~meta_norep$nicoya.x) # 0.00383
kruskal.test(meta_norep$Gran~meta_norep$nicoya.x) #0.04855
kruskal.test(meta_norep$NK~meta_norep$nicoya.x) #0.4657
kruskal.test(meta_norep$CD4T~meta_norep$nicoya.x) #0.2626
kruskal.test(meta_norep$Bcell~meta_norep$nicoya.x) #0.2871
kruskal.test(meta_norep$Mono~meta_norep$nicoya.x) #0.7887

#horvath estimated proportions (fig 1b)
kruskal.test(meta_norep$CD8.naive~meta_norep$nicoya.x) #0.008789
kruskal.test(meta_norep$CD8pCD28nCD45RAn~meta_norep$nicoya.x) #0.02187

#correlation between estimated horvath proportions and chronological age
cor.test(subset(meta_norep, nicoya.x =="1")[,"CD8.naive"], log(subset(meta_norep, nicoya.x =="1")[,"age.x"])) #r=-0.547933
cor.test(subset(meta_norep, nicoya.x =="1")[,"CD8pCD28nCD45RAn"], log(subset(meta_norep, nicoya.x =="1")[,"age.x"])) #r=0.6165151 
cor.test(subset(meta_norep, nicoya.x =="0")[,"CD8.naive"], log(subset(meta_norep, nicoya.x =="0")[,"age.x"])) #r=-0.6037924 
cor.test(subset(meta_norep, nicoya.x =="0")[,"CD8pCD28nCD45RAn"], log(subset(meta_norep, nicoya.x =="0")[,"age.x"])) #r=0.6113757 
```

