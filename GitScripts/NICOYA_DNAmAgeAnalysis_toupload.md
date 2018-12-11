---
title: "NICOYA_DNAmAgeAnalysis"
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


1. Prepare data to upload to Horvath epigenetic clock...
```{r eval = FALSE}
**You will need the following files for this**
- [datMiniAnnotation27k.csv](https://github.com/kobor-lab/450k-Analysis/blob/master/analysis%20scripts/DNAm%20Age%20-%20Horvath/datMiniAnnotation27k.csv)

datMiniAnnotation=read.csv("~/R Scripts/datMiniAnnotation.csv") #get from horvath website above

#make function to generate data for online tool 
horvathPrep <- function(dat){
  dat <- as.data.frame(dat)
dat0= cbind(rownames(dat), dat) #need col 1 to be CpG probe ID
datMiniAnnotation=read.csv("~/R Scripts/datMiniAnnotation.csv") #get from horvath website above
match1=match(datMiniAnnotation[,1], dat0[,1] )
dat0Reduced=dat0[match1,]
dat0Reduced[,1]=as.character(dat0Reduced[,1])
dat0Reduced[is.na(match1),1]= as.character(datMiniAnnotation[is.na(match1),1])
datout=data.frame(dat0Reduced)
# make sure you output numeric variables...
for (i in 2:dim(datout)[[2]]  ){datout[,i]= as.numeric(as.character(gsub(x=datout[,i],pattern="\"",replacement=""))) }
colnames(datout)[1] <- "Probe"                         
return(datout)}

#run data prep function (ran on not cell type corrected data because this is a pan-tissue clock)
betas <- betas(CRELES_notcelltypecorrected)
inputdata <- horvathDNAmAge(CRELES_notcelltypecorrected)

write.table(betas_horvath,"betas_horvath.csv", row.names=F, sep=",") #added this data to meta file
```

```{r}
load("~/FinalObjects/meta_norep.RData")

#correlation between chronological age and DNAm ages...
cor.test(meta_norep$age.x, meta_norep$DNAmAge, method = "spearman") #0.8687636 #horvath output "DNAmAge"
cor.test(meta_norep$age.x, meta_norep$Weid99, method = "spearman") #0.84949 #this was calculated separately according to  Q Lin - 2016 method #(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4789590/)
cor.test(meta_norep$age.x, meta_norep$BioAge1HA, method = "spearman") #0.8524417, calibrated version that is from horvath output

mean(meta$AgeAccelerationDiff) #-6.9 years
```


