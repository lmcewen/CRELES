---
title: "NICOYA_GLINT"
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

# GLINT
reference: [http://safrabio.cs.tau.ac.il/download/Papers/324.pdf](http://safrabio.cs.tau.ac.il/download/Papers/324.pdf)

```{r}
#creating required files
load("/big_data/lmcewen/Rehkopf/pdat_correctedjune2016.RData")
pdat <- pdat_correctedjune2016
mycovar  <- as.data.frame(cbind(ID = rownames(pdat), Group = pdat$nicoya.y, Sex = pdat$male.x, CD8Naive = as.numeric(pdat$CD8.naive), CD8memory = as.numeric(pdat$CD8pCD28nCD45RAn), CD4 = as.numeric(pdat$CD4T), Gran= as.numeric(pdat$Gran), NK= as.numeric(pdat$NK), Mono=as.numeric(pdat$Mono), PBlast=as.numeric(pdat$PlasmaBlast),Age= as.numeric(pdat$age.x)))
mypheno <- as.data.frame(cbind(ID = rownames(pdat),Group = pdat$nicoya.x, Age =as.numeric(pdat$age)))

#this is how I saved my files (hardest part was trying to figure out the right format!):
write.table(mydatafile_nov10, file = "~/rehkopf glint/datafile_nov10.txt", sep = "\t", row.names=F, col.names=T, quote = F)
write.table(mycovar, file = "~/rehkopf glint/mycovar_nov10.txt", sep = "\t", row.names=F, col.names=T,quote = F)
write.table(mypheno, file = "~/rehkopf glint/mypheno_nov10.txt", sep = "\t", row.names=F, col.names=T,quote = F)

#final files:
write.table(mydatafile_nov10, file = "~/rehkopf glint/datafile_nov10.txt", sep = "\t", row.names=F, col.names=T, quote = F)
write.table(mycovar, file = "~/rehkopf glint/mycovar_nov10.txt", sep = "\t", row.names=F, col.names=T,quote = F)
write.table(mypheno, file = "~/rehkopf glint/mypheno_nov10.txt", sep = "\t", row.names=F, col.names=T,quote = F)

#in terminal python
python glint.py --datafile datafile_nov10.txt --covarfile mycovar_nov10.txt --phenofile mypheno_nov28.txt --gsave
 
python glint.py --datafile datafile_nov10.glint --plot --plotpcs --numpcs 2 --out pcs_plot_nov10

python glint.py --datafile datafile_nov10.glint --maxpcstd 1 4 --gsave --out mydata_cleaned_nov10

python glint.py --datafile datafile_nov10.glint --epi --covar Sex CD8Naive CD8memory CD4 Gran NK Mono PBlast Age --savepcs 6 --gsave --out mydata_final_nov10

#back in R
mydata_final_nov10.samples <- read.delim("~/rehkopf glint/mydata_final_nov10.samples.txt", header=F, comment.char="#")
load("/big_data/lmcewen/Rehkopf/pdat_correctedjune2016.RData")
pdat_correctedjune2016$Epi1 <- mydata_final_nov10.samples$V14
pdat_correctedjune2016$Epi2 <- mydata_final_nov10.samples$V15

ggplot(pdat_correctedjune2016, aes(x=EPI1, y=EPI2, colour = as.factor(nicoya.x)))+
  scale_colour_manual(values=c("black", "dodgerblue"))+
  geom_point()+
  theme_bw()+
    theme(axis.text=element_text(colour="black", size = 14), 
        axis.title=element_text(face = "bold", colour="black", size = 14),
        title=element_text(face = "bold", colour="black", size = 14), 
        legend.key.size = unit(1, "cm"), legend.title = element_text(size = 14), 
        legend.text = element_text(size = 14))+
 labs(colour="Group")+
  ylab("PC2")+
  xlab("PC1")
  ```



