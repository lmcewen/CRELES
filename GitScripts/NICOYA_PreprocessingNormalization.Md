---
title: "Nicoya - Preprocessing & Normalization"
author: "Lisa M. McEwen"
date: "January 9, 2017"
output: html_document
---

Please email: lmcewen@bcchr.ca if you have any questions about this code. 

### Publication details:
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
1. (locally) Kobor Server: big_data/lmcewen/Rehkopf
2. (by request) Stanford Server: please email drehkopf@stanford.edu or lmcewen@bcchr.ca

450K Data preprocessing
======

## Load libraries & files

```{r eval = FALSE}
#libraries
setwd("/big_data/lmcewen/Rehkopf")
library(minfi)
library(FlowSorted.Blood.450k)
library(IlluminaHumanMethylation450kmanifest)
library(quadprog)
library(limma)
library(methylumi)
library(lumi)
library(gridExtra)
library(gplots)
library(ggdendro)
library(marray)
library(lumi)
library(lattice)
library("RColorBrewer")
library(knitr)
library(xtable)
library(plyr)
library(compositions)
install.packages('dendextend')
library(dendextend)
source("http://bioconductor.org/biocLite.R")
biocLite("FlowSorted.Blood.450k")
biocLite("quadprog")
biocLite("IlluminaHumanMethylation450kmanifest")
source('~/GitScripts/mypreprocessingfunctions.R')

#objects
allDat <- ("/home/lmcewen/FinalObjects/CRELES_allDat.txt") #exported from GenomeStudio software, already colour corrected and background subtracted)
avgBeta <- ("OriginalFiles/CRELES_avgBeta.txt")
qcDat <- ("OriginalFiles/CRELES_qcDat.txt")
```

Create methylumi objects...
```{r eval = FALSE}
dim(CRELES.raw <- lumiMethyR(allDat)) # probes = 485,577, n = 96
dim(CRELES.2.raw <- methylumiR(avgBeta, qcfile = qcDat)) # probes = 485,577, n = 96
sampleNames(CRELES.raw) <- gsub("_R2","", sampleNames(CRELES.raw))
sampleNames(CRELES.raw) <- gsub(" ","",sampleNames(CRELES.raw))
meta_withrep<-read.csv("/home/lmcewen/FinalObjects/meta_withrep.csv",row.names=1)
meta_withrep <- meta_withrep[sampleNames(CRELES.raw),]
identical(rownames(meta_withrep), sampleNames(CRELES.raw))
pData(CRELES.raw) <- meta_withrep

save(CRELES.raw, file ="/home/lmcewen/FinalObjects/CRELES.raw.RData")
```


## Initial QC
1. Denodogram of samples based on SNP probes
2. MDS plot of samples to illustrate outliers
3. Check beta and m-value distributions


### 1. Denodogram of samples based on SNP probes

```{r eval = FALSE}
CRELES.rs <- CRELES.raw[substring(featureNames(CRELES.raw), 1, 2) == "rs", ] #SNP probes = 65, n = 96
dend <- t(exprs(CRELES.rs)) %>%  scale %>% 
        dist %>% hclust %>% as.dendrogram
dend %>% set("leaves_pch", 19) %>%  # node point type
  set("leaves_cex", 1) %>%  # node point size
  set("leaves_col", "purple") %>% 
  set("labels_cex", .75) %>%
  plot()
```
![SNP Dendrogram](https://github.com/lmcewen/CRELES/blob/master/GitFigs/CRELES_preprocessing_SNPs.png)

We do this to double check that our one technical replicate ("1659_Rep1" & "1659_Rep2") clusters together. We see here that this is the case, so based on the SNPs there are no obvious outliers at this point.


### 2. MDS plot to identify outliers

First we colour by sex to ensure no sample mixups...
```{r eval = FALSE}
sexCol <- as.factor(pdat2$male.y)
plotSampleRelation(CRELES.raw, method = "mds", col=sexCol)
hist(test[,"1485"]) #clearly this is actually a male and NOT a female. Females would have a tridodal beta value distribution because of their X chromosome. For example:
fdat <- fData(CRELES.raw)
hist(betas.raw[rownames(subset(fdat, CHR == "X")),"1485"], main = "Incorrectly labelled sample - should be male", xlab= "beta value distribution of X chromosome")
```
![MDS Plot](https://github.com/kobor-lab/Publication-Scripts/blob/master/Lisa%20McEwen/Nicoya%20Paper/GitFigs/CRELES_preprocessingMDS.png)

One sample is mislabelled (sample "1485"), seen on the right there is sample that is black clustering with the red (red = males, black = females). There seems to be a male mislabelled as a female in this crowd. We can investigate this more given females have X-inactivation of their X chromosome, which can be visualized by plotting the methylation of the X chromosome (females should be trimodal, males should be bimodal):

![Incorrectly Labelled Sample](https://github.com/kobor-lab/Publication-Scripts/blob/master/Lisa%20McEwen/Nicoya%20Paper/GitFigs/CRELES_preprocessing_mislabelledSample.png)

^^This should actually be a male given the x chromosome distribution. This was updated in the metadata seperately.

For reference, this is properly labelled female...
```{r eval = FALSE}
hist(betas.raw[rownames(subset(fdat, CHR == "X")),"2055"], main = "Correctly labelled female", xlab= "beta value distribution of X chromosome") #SEE
```
![Correctly Labelled Sample](https://github.com/kobor-lab/Publication-Scripts/blob/master/Lisa%20McEwen/Nicoya%20Paper/GitFigs/CRELES_preprocessing_correctlyLabelled.png)

### 3. Check distributions of betas & m-values
```{r eval=FALSE}
#subset a random 20000 probes to plot
pdat_withRep$SampleID <- rownames(pdat_withRep) 

#BETAS
betas_sub <-betas.raw[sample(1:nrow(betas.raw), 200000),]
betas_sub_m <- melt(betas_sub)
plot_beta<-merge(betas_sub_m,pdat_withRep, by.x="Var2", by.y="SampleID")

ggplot(plot_beta, aes(value, group=Var2, color=as.factor(nicoya.x)))+
  geom_density()+theme_bw()+
  theme(axis.text = element_text(size =10, color="black"),
            axis.text.x = element_text(),
          axis.title = element_text(size =15))+
  scale_color_aaas(guide = guide_legend(title = "Nicoya Group"))
```

![raw beta distribution](https://github.com/kobor-lab/Publication-Scripts/blob/master/Lisa%20McEwen/Nicoya%20Paper/GitFigs/CRELES_preprocessing_rawbeta_dist.png)

```{r eval = FALSE}
#MVALS
rawMvals <- exprs(CRELES.raw)
mvals_sub <-rawMvals[sample(1:nrow(rawMvals), 200000),]
mvals_sub_m <- melt(mvals_sub)
plot_mvals<-merge(mvals_sub_m,pdat_withRep, by.x="Var2", by.y="SampleID")

ggplot(plot_mvals, aes(value, group=Var2, color=as.factor(nicoya.x)))+
  geom_density()+theme_bw()+
  theme(axis.text = element_text(size =10, color="black"),
            axis.text.x = element_text(),
          axis.title = element_text(size =15))+
  scale_color_manual(values = c("darkorchid4","forestgreen"), guide = guide_legend(title = "Nicoya Group"))
```

![raw m-value distribution](https://github.com/kobor-lab/Publication-Scripts/blob/master/Lisa%20McEwen/Nicoya%20Paper/GitFigs/CRELES_preprocessing_rawmvas_dist.png)




## Probe Filtering
Based on [EM Price et al.,2013](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/1756-8935-6-4)
1. Remove SNP rs probes 
2. Remove probes targeted to XY chromosomes
3. Remove poor performing probes: bad detection pvalue & bad avg beta probes
4. Remove probes targeted to polymorphic CpGs
5. Remove probes that are predicted to cross hybridize to other areas of the genome

### 1. Remove SNP rs probes
The 450K array has 65 rs probes designed for control purposes. 
```{r eval = FALSE}
dim(CRELES.preprocessed <- CRELES.raw[substring(featureNames(CRELES.raw),1,2) != "rs", ]) # probes = 485512, samples = 96 
```

### 2. Remove probes targeted to XY chromosomes
```{r eval = FALSE}
dim(CRELES.preprocessed <- CRELES.preprocessed[!fData(CRELES.preprocessed)$CHR%in%c("X", "Y"), ]) #probes = 473864, samples = 96 
```

### 3. Remove poor performing probes: bad detection pvalue & bad avg beta probes
```{r eval = FALSE}
sampleNames(CRELES.2.raw) <- sampleNames(CRELES.raw) 
avgbeta <- betas(CRELES.2.raw)[featureNames(CRELES.preprocessed),sampleNames(CRELES.preprocessed)]
badAvgbeta <- is.na(avgbeta)
nbadAvgbeta <- print(sum(rowSums(badAvgbeta)>=5))  #792
dotchart(colSums(badAvgbeta), cex=0.4, main = paste("Bad Avg Beta: probes = 473,864"))
badDetP <- detection(CRELES.preprocessed)>0.01
nbadDetP <- print(sum(rowSums(badDetP)>=5)) #1235
dotchart(colSums(badDetP),main=paste("Bad Det P, probes = 473,684"), cex = 0.4)
sum(badProbes <- rowSums(badAvgbeta)>=5|rowSums(badDetP)>=5)# badProbes = 1,976
dim(CRELES.preprocessed <- CRELES.preprocessed[!badProbes,]) # probes = 471,888, n = 96
```

### 4. Remove probes targeted to polymorphic CpGs
```{r eval = FALSE}
snp_index <- which(is.na(fData(CRELES.preprocessed)$n_SNPCpG)) #removed 20134 polymorphic probes
dim(CRELES.preprocessed <- CRELES.preprocessed[snp_index, ]) #451,754, 96
```

### 5. Remove probes that are predicted to cross hybridize to other areas of the genome
```{r eval = FALSE}
xy_hit_index <- which(fData(CRELES.preprocessed)$XY_Hits == "XY_NO") #10645 cross hybridizing probes 
dim(CRELES.preprocessed <- CRELES.preprocessed[xy_hit_index, ]) #441,109 probes, 96 samples 
```

### Save preprocessed objects
```{r eval = FALSE}
save(CRELES.preprocessed, file = "/home/lmcewen/FinalObjects/CRELES.preprocessed.RData")
save(CRELES.2.raw, file = "/home/lmcewen/FinalObjects/CRELES.2.raw.RData")

```


## Data normalization
1. Look at heat scree PCA plot before norm
2. Probe type normalization - SWAN
3. Correct for technical variation (sentrid ID/row) - ComBat
4. Estimate cell type proportions - look at PCA plot - remove variance


### 1. Look at heat scree PCA plot before norm
```{r eval = FALSE}
load("/home/lmcewen/FinalObjects/CRELES.preprocessed.RData")
source('~/GitScripts/mypreprocessingfunctions.R')

## PCA
CRELES_processed_betas <- betas(CRELES.preprocessed)
PCA_full<-princomp(na.omit(CRELES_processed_betas))
Loadings<-as.data.frame(unclass(PCA_full$loadings))
vars <- PCA_full$sdev^2
Importance<-vars/sum(vars)

#Specify which covariates are categorical and/or categorical
meta <- pdat_withRep
meta$Sentrix_ID <- as.factor(meta$Sentrix_ID)
meta$Row <- as.factor(meta$Row)
meta$nicoya.x <- as.factor(meta$nicoya.x)
meta$male.x <- as.factor(meta$male.x)
meta$lowEdu <- as.factor(meta$lowEdu)
meta$lowWealth<- as.factor(meta$lowWealth)

meta_categorical<-meta[,c("Sentrix_ID", "Row", "nicoya.x","male.x", "lowEdu", "lowWealth")]
meta_continuous<-meta[,c("bmi.x","age","crp","telobp.x","kneeAvg")]



#Apply function on PCA results, pulls in the meta data and beta values from above
heat_scree_plot(Loadings, Importance)

```
![PCA heat scree plot before SWAN](https://github.com/kobor-lab/Publication-Scripts/blob/master/Lisa%20McEwen/Nicoya%20Paper/GitFigs/CRELES_normalization_heatscreeprenorm.png)


### 2a. Probe type normalization - SWAN
[Jovana Maksimovic et al., 2012](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-6-r44)
```{r eval = FALSE}
load("/home/lmcewen/GitScripts/nCpG all.RData")
load("/home/lmcewen/FinalObjects/CRELES.preprocessed.RData")
load("/home/lmcewen/FinalObjects/CRELES.2.raw.RData")
dim(CRELES.preprocessed) # Should have 441,109 probes and 96 samples
dim(CRELES.2.raw) # Should have 485,577 probes and 96 samples

target <- fData(QCdata(CRELES.2.raw))$TargetID =="NEGATIVE"
grn <- methylated(QCdata(CRELES.2.raw))[target,]
red <- unmethylated(QCdata(CRELES.2.raw))[target,]
bkg <- rowMeans(cbind(colMeans(grn),colMeans(red)))
all(names(bkg)%in%sampleNames(CRELES.preprocessed)) #true
bkg <- bkg[sampleNames(CRELES.preprocessed)]

all(featureNames(CRELES.preprocessed)%in%names(nCpG.all))
fData(CRELES.preprocessed)$nCpG <- nCpG.all[featureNames(CRELES.preprocessed)]
CRELES.swan <- preprocessMSWAN(CRELES.preprocessed,bg=bkg)
sum(is.na(exprs(CRELES.swan))) # should be 0
sum(is.na(unmethylated(CRELES.swan))) # should be 0
sum(is.na(methylated(CRELES.swan))) # should be 0
save(CRELES.swan,file = "/home/lmcewen/FinalObjects/CRELES.swan.RData")
```


### 2b. Look at heat scree PCA plot after SWAN
```{r eval = FALSE}
load("/home/lmcewen/FinalObjects/CRELES.swan.RData")
meta <- pdat_withRep[colnames(CRELES.swan),]
meta_categorical<-meta[,c("Sentrix_ID", "Row", "nicoya.x","male.x", "lowEdu", "lowWealth")]
meta_continuous<-meta[,c("bmi.x","age","crp","telobp.x","kneeAvg")]
ord <- 1:11
CRELES_swan_betas <- betas(CRELES.swan)
PCA_full<-princomp(na.omit(CRELES_swan_betas))
Loadings<-as.data.frame(unclass(PCA_full$loadings))
vars <- PCA_full$sdev^2
Importance<-vars/sum(vars)
heat_scree_plot(Loadings, Importance)
```
![PCA heat scree plot after SWAN](https://github.com/kobor-lab/Publication-Scripts/blob/master/Lisa%20McEwen/Nicoya%20Paper/GitFigs/CRELES_normalization_heatscree_afterSWAN.png)


### 3a. Correct for technical variation (sentrid ID) - ComBat
```{r eval = FALSE}
CRELES.swan.bat <- CRELES.swan
mod <- model.matrix(~as.factor(male) + as.factor(age) + as.factor(nicoya), data=metadata) #variables want to keep
test.bat<- ComBat(exprs(CRELES.swan), metadata$Sentrix_ID) #numCovs is the columns in which to be treated as continuous variables
exprs(CRELES.swan.bat) <- test.bat ##changing expression data to now be combat corrected for sentrix ID
save(CRELES.swan.bat, file = "/home/lmcewen/FinalObjects/CRELES.swan.bat.RData")
```


### 3b. Look at heat scree PCA plot after ComBat (sentrix ID)
```{r eval = FALSE}
load("/home/lmcewen/FinalObjects/CRELES.swan.bat.RData")
meta <- pdat_withRep[colnames(CRELES_swan.bat_betas),]
meta_categorical<-meta[,c("Sentrix_ID", "Row", "nicoya.x","male.x", "lowEdu", "lowWealth")]
meta_continuous<-meta[,c("bmi.x","age","crp","telobp.x","kneeAvg")]
ord <- 1:11

CRELES_swan.bat_betas <- betas(CRELES.swan.bat)
PCA_full<-princomp(na.omit(CRELES_swan.bat_betas))
Loadings<-as.data.frame(unclass(PCA_full$loadings))
vars <- PCA_full$sdev^2
Importance<-vars/sum(vars)
heat_scree_plot(Loadings, Importance)
```
![PCA heat scree plot after SWAN and combat for sentrix ID](https://github.com/kobor-lab/Publication-Scripts/blob/master/Lisa%20McEwen/Nicoya%20Paper/GitFigs/CRELES_normalization_heatscree_afterswanbat.png)

### 3c. Correct for technical variation (sentrix row) - ComBat
```{r eval = FALSE}
#ComBat for row
head(betas(CRELES.swan.bat))
CRELES.swan.bat.2 <- CRELES.swan.bat
test.bat.2 <- ComBat(exprs(CRELES.swan.bat), meta$Row) 

exprs(CRELES.swan.bat.2) <- test.bat.2 
save(CRELES.swan.bat.2, file = "/home/lmcewen/FinalObjects/CRELES.swan.bat.2.RData")

#remove technical replicate
replicate <- c("1659_Rep2")
CRELES.swan.bat.2<- CRELES.swan.bat.2[, !sampleNames(CRELES.swan.bat.2)%in%replicates]
sampleNames(CRELES.swan.bat.2) <- gsub("_Rep1",sampleNames(CRELES.swan.bat.2))
sampleNames(CRELES.swan.bat.2) <- gsub("_R2",sampleNames(CRELES.swan.bat.2))
CRELES_notcelltypecorrected <- CRELES.swan.bat.2

save(CRELES_notcelltypecorrected, file = "/home/lmcewen/FinalObjects/CRELES_notcelltypecorrected.RData")
```

### 3d. Look at heat scree PCA plot after ComBat (sentrix ID and row)
```{r}
load("/home/lmcewen/FinalObjects/CRELES_notcelltypecorrected.RData")

meta_categorical<-meta[,c("Sentrix_ID", "Row", "nicoya.x","male.x", "lowEdu", "lowWealth")]
meta_continuous<-meta[,c("bmi.x","age","crp","telobp.x","kneeAvg")]

CRELES_notcelltypecorrected_betas <- betas(CRELES_notcelltypecorrected)
PCA_full<-princomp(na.omit(CRELES_notcelltypecorrected_betas))
Loadings<-as.data.frame(unclass(PCA_full$loadings))
vars <- PCA_full$sdev^2
Importance<-vars/sum(vars)
heat_scree_plot(Loadings, Importance)

##remove technical replicate
replicate <- c("1659_Rep2")
CRELES.swan.bat.2<- CRELES.swan.bat.2[, !sampleNames(CRELES.swan.bat.2)%in%replicates]
sampleNames(CRELES.swan.bat.2) <- gsub("_Rep1",sampleNames(CRELES.swan.bat.2))
sampleNames(CRELES.swan.bat.2) <- gsub("_R2",sampleNames(CRELES.swan.bat.2))
CRELES_notcelltypecorrected <- CRELES.swan.bat.2
save(CRELES_notcelltypecorrected, file = "/home/lmcewen/FinalObjects/CRELES_notcelltypecorrected.RData")

meta <- meta[!rownames(meta)%in%replicate,]
rownames(meta) <- gsub("_Rep1","", rownames(meta))
```
![PCA heat scree plot after SWAN and combat for both sentrix ID and row](https://github.com/kobor-lab/Publication-Scripts/blob/master/Lisa%20McEwen/Nicoya%20Paper/GitFigs/CRELES_normalization_heatscree_afterswanbat2.png)

### 4a. Estimate cell type proportions - look at PCA plot - remove variance
```{r eval = FALSE}
CRELES_celltypecorrected <- CRELES_notcelltypecorrected 

#load all libraries
library(minfi)
library(FlowSorted.Blood.450k)
library(IlluminaHumanMethylation450kmanifest)
library(quadprog)
library(limma)

#deconvolution script
RGset = as.data.frame(betas(CRELES_celltypecorrected))
source("/home/lmcewen/R Scripts/ECC2_original.R")
detach("package:dplyr", unload=TRUE)
BiocGenerics::combine
cellprop <- ECC2(RGset)
save(cellprop,file="CRELES_predicted_cellprops.RData")
cellpropcounts <- cellprop$counts
write.csv(cellpropcounts, file = "cellpropcounts.csv")
         
diff <- as.data.frame(cellprop$counts)
identical(colnames(betas(CRELES.swan.bat.2.decon)), rownames(diff)) #TRUE
```

### 4b. look at PCA plot with cell type proportions
```{r eval = FALSE}
#merge estimated cell types with meta

meta <- cbind(meta, cellprop$counts)
colnames(meta)
meta_categorical<-meta[,c("Sentrix_ID", "Row", "nicoya.x","male.x", "lowEdu", "lowWealth")]
meta_continuous<-meta[,c("bmi.x","age","crp","telobp.x","kneeAvg","CD8T","Bcell","CD4T", "NK", "Gran", "Mono")]

ord <-1:17
CRELES_notcelltypecorrected_betas <- betas(CRELES_notcelltypecorrected)
PCA_full<-princomp(na.omit(CRELES_notcelltypecorrected_betas))
Loadings<-as.data.frame(unclass(PCA_full$loadings))
vars <- PCA_full$sdev^2
Importance<-vars/sum(vars)
heat_scree_plot(Loadings, Importance)
```
![PCA heat scree plot after SWAN and combat for both sentrix ID and row but include estimated cell types](https://github.com/kobor-lab/Publication-Scripts/blob/master/Lisa%20McEwen/Nicoya%20Paper/GitFigs/CRELES_normalization_heatscree_afterswanbat2_withcelltypes.png)

### 4c. Remove cell type associated variance
```{r eval = FALSE}
betas <- betas(CRELES_celltypecorrected)

# fit linear model for cell type effects on each probe in betas
betas.lm <- apply(betas, 1, function(x){
  blood <- diff[colnames(betas),]
  lm(x~0+CD8T+CD4T+NK+Bcell+Mono+Gran,data=blood)
})

# extract matrix of residuals from resulting linear models
residuals <- t(sapply(betas.lm, function(x)residuals(summary(x))))
colnames(residuals) <- colnames(betas) # re-name residuals columns with sample names

# add the residuals of each regression model to the mean methylation value of each probe (mean across all samples) to obtain the “adjusted” methylation data.
adj.residuals <- residuals+matrix(apply(betas, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))
# To make sure we do not induce any NAs into the dataset when we convert the beta values back M-values (by log2 transformation), we need to ensure we do not have any corrected beta values that are greater or equal to zero or any beta values that are greater than 1.

adj.residuals[adj.residuals<=0] <- 0.001 # convert any values that are less than or equal to zero to 0.001
adj.residuals[adj.residuals>1] <- 0.999 # convert any values that are greater than 1 to 0.999
adj.M.values<-beta2m(adj.residuals)
any(is.na(adj.M.values)) # should be FALSE indicating there are no NAs
summary(adj.residuals)

# replace adjusted data into a cloned methylumi object
identical(rownames(adj.M.values), rownames(betas(CRELES.swan.bat.2.decon))) #should be TRUE
identical(colnames(adj.M.values), colnames(betas(CRELES.swan.bat.2.decon))) #should be TRUE
exprs(CRELES_celltypecorrected) <- adj.M.values
all.equal(betas(CRELES_celltypecorrected), adj.residuals)

save(CRELES_celltypecorrected,file = "/home/lmcewen/FinalObjects/CRELES_celltypecorrected.RData")
```

### 4d. Look at heat scree PCA plot after cell type correction
```{r eval = FALSE}
load("/home/lmcewen/FinalObjects/CRELES_celltypecorrected.RData")

CRELES_celltypecorrected_betas <- betas(CRELES_celltypecorrected)
PCA_full<-princomp(na.omit(CRELES_celltypecorrected_betas))
Loadings<-as.data.frame(unclass(PCA_full$loadings))
vars <- PCA_full$sdev^2
Importance<-vars/sum(vars)
heat_scree_plot(Loadings, Importance)
```
![PCA heat scree plot after SWAN and combat for both sentrix ID and row but include estimated cell types](https://github.com/kobor-lab/Publication-Scripts/blob/master/Lisa%20McEwen/Nicoya%20Paper/GitFigs/CRELES_normalization_heatscree_afterswanbat2decon.png)


