---
title: "Nicoya - Preprocessing & Normalization"
author: "Lisa M. McEwen"
date: "January 9, 2017"
output: html_document
---

# CRELES Study

Please email: lmcewen@bcchr.ca if you have any questions about this code. 

### Publication details (accepted at Epigenetics and Chromatin on April 12, 2017):
__Title:__ Differential DNA methylation and lymphocyte proportions in a Costa Rican high longevity region

__Authors:__ Lisa M. McEwen, BSc; Alexander M. Morin, BSc.; Rachel D. Edgar, MSc.; Julia L. MacIsaac, PhD; Meaghan J. Jones, PhD.; William H. Dow, PhD.; Luis Rosero-Bixby, PhD.; Michael S. Kobor, PhD.; David H. Rehkopf, ScD., MPH.

__Date:__ April 2017

__Paper Link:__ [Click here](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-017-0128-2)


Overview
======
This R script is to preprocess and normalize data generated with the Infinium HumanMethylation450 Array (450K).

Whole blood was collected from 95 individuals from the [Costa Rican Berkeley CRELES cohort](http://www.creles.berkeley.edu/). We sampled from the Nicoyan population ('Nicoyans',n=47) and from other regions of Costa Rica('Controls' or 'non-Nicoyans',n=47). The Nicoya peninsula hosts a population of Costa Rica that exceed expected life expectancy, see [Luis Rosero-Bixby et al., 2013](http://pubmedcentralcanada.ca/pmcc/articles/PMC4241350/). The goal of this study was to explore the epigenetics (DNA methylation) of Nicoyan individuals as compared to non-Nicoyans, considering that DNA methylation is closely related to biological aging. 

Study Design: 
95 samples + 1 technical replicate, randomly distributed across 12 arrays ("chips", as each chip has 12 samples) run in one batch.

All R objects and files are located:
1. Kobor Server: big_data/lmcewen/Rehkopf (locally saved)
2. Stanford Server: please email drehkopf@stanford.edu or lmcewen@bcchr.ca to request access

LINKS TO CODE:

[Preprocessing & Normalization](https://github.com/lmcewen/CRELES/blob/master/PreprocessingNormalization.Md)

[DNAm Age Analysis](https://github.com/lmcewen/CRELES/blob/master/GitScripts/NICOYA_DNAmAgeAnalysis_toupload.md)

[Cell Type Analysis](https://github.com/lmcewen/CRELES/blob/master/GitScripts/NICOYA_CellTypeAnalyses_toupload.md)

[DMRcate Analysis](https://github.com/lmcewen/CRELES/blob/master/GitScripts/NICOYA_DMRCATE_toupload.md)

[DNAm Variability Analysis](https://github.com/lmcewen/CRELES/blob/master/GitScripts/NICOYA_DNAm_Variability_analyses_toupload.md)

[GLINT Analysis](https://github.com/lmcewen/CRELES/blob/master/GitScripts/Nicoya_glint.md)

[Pyro verification](https://github.com/lmcewen/CRELES/blob/master/GitScripts/Nicoya_pyro.md)
