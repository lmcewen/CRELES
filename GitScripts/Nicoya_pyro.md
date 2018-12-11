2016 Oct 27

Pyrosequencing Verification for CRELES - chose 3 CpGs from 4 significant sites

```{r}

#1 - cg02853387 

cg02853387_pyro <- read.csv("~/Rehkopf/pyro oct 2016/cg02853387_pyro1_final.csv", row.names=1)

cg02853387_pyro$betas450k <- betas_450k["cg02853387",rownames(cg02853387_pyro)] #1
cor.test(as.numeric(cg02853387_pyro$betas450k), as.numeric(cg02853387_pyro$cg02853387)) #0.8724391
cg02853387_pyro$cg02853387<-cg02853387_pyro$cg02853387/100
cg02853387_pyro$diff <- cg02853387_pyro$betas450k - cg02853387_pyro$cg02853387
cg02853387_pyro$avg <- sapply(1:nrow(cg02853387_pyro), function(x) mean(cg02853387_pyro[x,"betas450k"],cg02853387_pyro[x,"cg02853387"]))

m1<-mean(cg02853387_pyro$diff)
upper1<-m1+2*(sd(cg02853387_pyro$diff))
lower1<-m1-2*(sd(cg02853387_pyro$diff))

ggplot(cg02853387_pyro, aes(avg, diff, color =rownames(cg02853387_pyro), label = rownames(cg02853387_pyro))) +
  geom_text(size=3.5)+
  geom_hline(aes(yintercept=(m1))) + geom_hline(aes(yintercept=(upper1)),linetype="dashed") +
  geom_hline(aes(yintercept=(lower1)),linetype="dashed") + 
  theme_bw() + ylab("Diff (450K - Pyro)") + xlab("Average of 450K and Pyro Measures")+
  ggtitle("cg02853387")+
  theme(legend.position="none")
ggsave("~/Rehkopf/pyro oct 2016/cg02853387_blandaltman_final.pdf",scale =1)

ggplot(cg02853387_pyro, aes(x=cg02853387, y=betas450k, color =rownames(cg02853387_pyro), label = rownames(cg02853387_pyro))) +
  geom_text(size=3.5)+
   theme_bw()+
  stat_smooth(method = "lm", color = "grey",se=F)+
  theme(legend.position="none")
ggsave("~/Rehkopf/pyro oct 2016/cg02853387_scatter_final.pdf",scale =1)

cg02853387_test<- cbind(cg02853387_pyro,pdat_correctedjune2016[rownames(cg02853387_pyro),c("nicoya.y","male.x","CD8.naive","CD8pCD28nCD45RAn", "CD4T", "Gran", "NK", "Mono","PlasmaBlast","age.x")])#3
summary(aov(cg02853387 ~ nicoya.y+ male.x + age.x+CD8.naive+CD8pCD28nCD45RAn+CD4T+Gran+NK+Mono+PlasmaBlast, data = cg02853387_test)) #7.9e-05 

ggplot(cg02853387_test, aes(x=as.factor(nicoya.y),y=cg02853387))+
  geom_jitter(width = 0.25, size=1)+
  geom_boxplot(outlier.colour = NA, fill = c("lightgrey","cornflowerblue"))+
  theme_bw()+
  ylim(c(0,1))
  ggsave("~/Rehkopf/pyro oct 2016/cg02853387_boxplot_final.pdf",scale =1)
  
mean(cg02853387_pyro[rownames(subset(pdat_correctedjune2016, nicoya.x=="1")),"cg02853387"]) -
mean(cg02853387_pyro[rownames(subset(pdat_correctedjune2016, nicoya.x=="0")),"cg02853387"]) 
  
#2 - cg02438481
cg02438481_pyro <- read.csv("~/Rehkopf/pyro oct 2016/cg02438481_pyro2_final.csv", row.names=1)
cg02438481_pyro$betas450k <- betas_450k["cg02438481",rownames(cg02438481_pyro)]#2
cor.test(as.numeric(cg02438481_pyro$betas450k), as.numeric(cg02438481_pyro$cg02438481)) #0.92


cg02438481_pyro$cg02438481 <- cg02438481_pyro$cg02438481/100
cg02438481_pyro$diff <- cg02438481_pyro$betas450k - cg02438481_pyro$cg02438481
cg02438481_pyro$avg <- sapply(1:nrow(cg02438481_pyro), function(x) mean(cg02438481_pyro[x,"betas450k"], cg02438481_pyro[x,"cg02438481"]))

m2<-mean(cg02438481_pyro$diff)
upper2<-m2+2*(sd(cg02438481_pyro$diff))
lower2<-m2-2*(sd(cg02438481_pyro$diff))

ggplot(cg02438481_pyro, aes(avg, diff, color =rownames(cg02438481_pyro),  label = rownames(cg02438481_pyro))) +
  geom_text(size=3.5)+
  geom_hline(aes(yintercept=(m2))) + geom_hline(aes(yintercept=(upper2)),linetype="dashed") +
  geom_hline(aes(yintercept=(lower2)),linetype="dashed") + 
  theme_bw() + ylab("Diff (450K - Pyro)") + xlab("Average of 450K and Pyro Measures")+
  ggtitle("cg02438481")+
  theme(legend.position="none")# lines represent mean difference + or - 2sd
ggsave("~/Rehkopf/pyro oct 2016/cg02438481_blandaltman_final.pdf",scale =1)

ggplot(cg02438481_pyro, aes(x=cg02438481, y=betas450k, color =rownames(cg02438481_pyro), label = rownames(cg02438481_pyro))) +
  geom_text(size=3.5)+
   theme_bw()+
  stat_smooth(method = "lm", color = "grey",se=F)+
  theme(legend.position="none")
ggsave("~/Rehkopf/pyro oct 2016/cg02438481_scatter_final.pdf",scale =1)

cg02438481_test<- cbind(cg02438481_pyro,pdat_correctedjune2016[rownames(cg02438481_pyro),c("nicoya.y","male.x","CD8.naive","CD8pCD28nCD45RAn", "CD4T", "Gran", "NK", "Mono","PlasmaBlast","age.x")])#3
summary(aov(cg02438481 ~ nicoya.y+ male.x + age.x+CD8.naive+CD8pCD28nCD45RAn+CD4T+Gran+NK+Mono+PlasmaBlast, data = cg02438481_test)) #1.94e-05, 

ggplot(cg02438481_test, aes(x=as.factor(nicoya.y),y=cg02438481))+
  geom_jitter(width = 0.25, size=1)+
  geom_boxplot(outlier.colour = NA, fill = c("lightgrey","cornflowerblue"))+
  theme_bw()+
  ylim(c(0,1))
  ggsave("~/Rehkopf/pyro oct 2016/cg02438481_boxplot_final.pdf",scale =1)

 mean(cg02438481_pyro[rownames(subset(pdat_correctedjune2016, nicoya.x=="1")),"cg02438481"]) -
mean(cg02438481_pyro[rownames(subset(pdat_correctedjune2016, nicoya.x=="0")),"cg02438481"]) 


#cg13979274_pyro 3
cg13979274_pyro <- read.csv("~/Rehkopf/pyro oct 2016/cg13979274_pyro_final.csv", row.names=1)
head(cg13979274_pyro)
cg13979274_pyro$betas450k <- betas_450k["cg13979274",rownames(cg13979274_pyro)]#3
cor.test(cg13979274_pyro$cg13979274, cg13979274_pyro$betas450k) #.8801523

cg13979274_pyro$cg13979274 <- cg13979274_pyro$cg13979274/100
cg13979274_pyro$diff <-  cg13979274_pyro[,2] - cg13979274_pyro[,1]
cg13979274_pyro$avg <- sapply(1:nrow(cg13979274_pyro), function(x) mean(cg13979274_pyro[x,1], cg13979274_pyro[x,2]))

m3<-mean(cg13979274_pyro$diff)
upper3<-m3+2*(sd(cg13979274_pyro$diff))
lower3<-m3-2*(sd(cg13979274_pyro$diff))

ggplot(cg13979274_pyro, aes(avg, diff, color =rownames(cg13979274_pyro), label = rownames(cg13979274_pyro))) +
  geom_text(size=3.5)+
  geom_hline(aes(yintercept=(m3))) + geom_hline(aes(yintercept=(upper3)),linetype="dashed") +
  geom_hline(aes(yintercept=(lower3)),linetype="dashed") + 
  theme_bw() + ylab("Diff (450K - Pyro)") + xlab("Average of 450K and Pyro Measures")+
  ggtitle("cg13979274")+
  theme(legend.position="none")# lines represent mean difference + or - 2sd
ggsave("~/Rehkopf/pyro oct 2016/cg13979274_blandaltman_final.pdf",scale =1)

ggplot(cg13979274_pyro, aes(x=cg13979274, y=betas450k, color =rownames(cg13979274_pyro), label = rownames(cg13979274_pyro))) +
  geom_text(size=3.5)+
   theme_bw()+
  stat_smooth(method = "lm", color = "grey",se=F)+
  theme(legend.position="none")
  ggsave("~/Rehkopf/pyro oct 2016/cg13979274_scatter_final.pdf",scale =1)

  
cg13979274_test<- cbind(cg13979274_pyro,pdat_correctedjune2016[rownames(cg13979274_pyro),c("nicoya.y","male.x","CD8.naive","CD8pCD28nCD45RAn", "CD4T", "Gran", "NK", "Mono","PlasmaBlast","age.x")])#3
summary(aov(cg13979274 ~ nicoya.y+ male.x + age.x+CD8.naive+CD8pCD28nCD45RAn+CD4T+Gran+NK+Mono+PlasmaBlast, data = cg13979274_test)) #9.1e-07
mean(cg13979274_pyro[rownames(subset(pdat_correctedjune2016, nicoya.x=="1")),"cg13979274"]) -
mean(cg13979274_pyro[rownames(subset(pdat_correctedjune2016, nicoya.x=="0")),"cg13979274"])   

ggplot(cg13979274_test, aes(x=as.factor(nicoya.y),y=cg13979274))+
  geom_jitter(width = 0.25, size=1)+
  geom_boxplot(outlier.colour = NA, fill = c("lightgrey","cornflowerblue"))+
  theme_bw()+
  ylim(c(0,1))
  ggsave("~/Rehkopf/pyro oct 2016/cg13979274_boxplot_final.pdf",scale =1)
```




* PYRO SEQUENCING PRIMERS *

Can also be found in supp: https://static-content.springer.com/esm/art%3A10.1186%2Fs13072-017-0128-2/MediaObjects/13072_2017_128_MOESM6_ESM.pdf

*cg02853387*

- Forward:  AGGGAAGAAAAGTTATTAAGTTGT

- Reverse (5’biotinylated): ACAAATACAAAACCCATATTCTCAA

- Sequencing: GTGTAGGTTTTTAGTTTATAGT

*cg02438481*

- Forward: GTTTTGGGTTTGGTGATTTGGTTTTA

- Reverse(5’biotinylated): ATTTCTTAATCAATACCACCTTCTTCTATA

- Sequencing: TTTATTTTAGGTGGGAGT

*cg13979274*

- Forward(5’biotinylated): AGGGGAGTATTTTAGTTTAGTGTATAG

- Reverse: CCAACTTAAAAAAACCAAACTTCAATATC

- Sequencing: AAAACAATTACAACCCTC
