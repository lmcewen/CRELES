#MDS PLOT
mymdsplot <- function(x, subset=NULL, cv.Th=0.1, standardize=TRUE, dimension=c(1,2), color=NULL, main=NULL, pch=NULL, addLegend=TRUE, ...) {
  
  if (is(x, 'ExpressionSet')) {
    dataMatrix <- exprs(x)
  } else if (is.matrix(x)) {
    dataMatrix <- x
  } else {
    stop('The class of "x" should be matrix or LumiBatch!')
  }
  
  if (standardize) 
    dataMatrix <- scale(dataMatrix)
  if (is.null(subset)) {
    probeList <- rownames(dataMatrix)
    if (is.null(probeList)) probeList <- 1:nrow(dataMatrix)
    if (cv.Th > 0) {
      cv.gene <- apply(dataMatrix, 1, function(x) sd(x)/mean(x))
      subset <- probeList[abs(cv.gene) > cv.Th]
      if (is.null(main)) main <- paste('Sample relations based on', length(subset), 'probes with sd/mean >', cv.Th)
    } else {
      subset <- probeList
      if (is.null(main)) main <- paste('Sample relations based on', length(subset), 'probes')
    }
  } else {
    if (length(subset) == 1 && is.numeric(subset)) {
      subset <- sample(1:nrow(dataMatrix), min(subset, nrow(dataMatrix)))
    }
    if (is.null(main)) main <- paste('Sample relations based on', length(subset), 'selected probes')
  }
  
  dd <- dist(t(dataMatrix[subset,]))
  
  mds.result <- cmdscale(dd, k=max(dimension), eig=TRUE)
  ppoints <- mds.result$points
  eig <- mds.result$eig
  percent <- round(eig/sum(eig) * 100, 1)
  
  colorLegend <- NULL
  if (is.null(color)) {
    color <- 1
  } else {
    if (!is.numeric(color)) {
      allColor <- colors()
      if (!all(is.element(color, allColor))) {
        colorLegend <- unique(color)
        color <- as.numeric(factor(color, levels=colorLegend))
      } 
    }
  }
  if (missing(pch)) {
    plot(ppoints[,dimension[1]], ppoints[,dimension[2]], type='n', 
         xlab=paste('Principal Component ', dimension[1], " (", percent[dimension[1]], "%)", sep=""),
         ylab=paste('Principal Component ', dimension[2], " (", percent[dimension[2]], "%)", sep=""), 
         main=main, ...)
    text(ppoints[,dimension[1]], ppoints[,dimension[2]], col=color, labels=colnames(dataMatrix), cex=0.75)
  } else {
    plot(ppoints[,dimension[1]], ppoints[,dimension[2]],  
         xlab=paste('Principal Component ', dimension[1], " (", percent[dimension[1]], "%)", sep=""),
         ylab=paste('Principal Component ', dimension[2], " (", percent[dimension[2]], "%)", sep=""), 
         main=main, col=color, pch=pch, ...)
  }
  attr(ppoints, 'geneNum') <- length(subset)
  attr(ppoints, 'threshold') <- cv.Th
  
  ## add legend if color is a factor
  if (!is.null(colorLegend) && addLegend) {
    if (!missing(pch)) {
      legend('topleft', legend=colorLegend, col=unique(color), pch=unique(pch))
    } else {
      legend('topleft', legend=colorLegend, text.col=1:length(colorLegend))
    }
  }
  
  return(invisible(ppoints))	
}

#HEAT SCREE PLOT PCA
heat_scree_plot<-function(Loadings, Importance){
  adjust<-1-Importance[1]
  pca_adjusted<-Importance[2:length(Importance)]/adjust
  pca_df<-data.frame(adjusted_variance=pca_adjusted, PC=seq(1:length(pca_adjusted)))
  
  scree<-ggplot(pca_df[which(pca_df$PC<25),],aes(PC,adjusted_variance))+geom_bar(stat = "identity",color="black",fill="grey")+theme_bw()+
    theme(axis.text = element_text(size =12),
          axis.title = element_text(size =15),
          plot.margin=unit(c(1,1.5,0.2,2.25),"cm"))+ylab("Adjusted Variance")+
    scale_x_continuous(breaks = seq(1,24,1))
  
  #### Heat
  ## correlate meta with PCS
  ## Run anova of each PC on each meta data variable
  
  aov_PC_meta<-lapply(1:ncol(meta_categorical), function(covar) sapply(1:ncol(Loadings), function(PC) summary(aov(Loadings[,PC]~meta_categorical[,covar]))[[1]]$"Pr(>F)"[1]))
  cor_PC_meta<-lapply(1:ncol(meta_continuous), function(covar) sapply(1:ncol(Loadings), function(PC) suppressWarnings(cor.test(Loadings[,PC],as.numeric(meta_continuous[,covar]),alternative = "two.sided", method="spearman", na.action=na.omit)$p.value)))
  names(aov_PC_meta)<-colnames(meta_categorical)
  names(cor_PC_meta)<-colnames(meta_continuous)
  aov_PC_meta<-do.call(rbind, aov_PC_meta)
  cor_PC_meta<-do.call(rbind, cor_PC_meta)
  aov_PC_meta<-rbind(aov_PC_meta, cor_PC_meta)
  aov_PC_meta<-as.data.frame(aov_PC_meta)
  #adjust
  aov_PC_meta_adjust<-aov_PC_meta[,2:ncol(aov_PC_meta)]
  
  #reshape
  avo<-aov_PC_meta_adjust[,1:24]
  avo_heat_num<-apply(avo,2, as.numeric)
  avo_heat<-as.data.frame(avo_heat_num)
  colnames(avo_heat)<-sapply(1:24, function(x) paste("PC",x, sep=""))
  avo_heat$meta<-rownames(avo)
  avo_heat_melt<-melt(avo_heat, id=c("meta"))
  
  # cluster meta data
  meta_var_order<-unique(avo_heat_melt$meta)[rev(ord)]
  avo_heat_melt$meta <- factor(avo_heat_melt$meta, levels = meta_var_order)
  
  # color if sig
  # avo_heat_melt$Pvalue<-sapply(1:nrow(avo_heat_melt), function(x) if(avo_heat_melt$value[x]>=0.9){">=0.9"}else{
  # if(avo_heat_melt$value[x]>=0.5){">=0.5"}else{
  # if(avo_heat_melt$value[x]>=0.1){">=0.1"}else{"<0.1"}}})
  avo_heat_melt$Pvalue<-sapply(1:nrow(avo_heat_melt), function(x) if(avo_heat_melt$value[x]<=0.001){"<=0.001"}else{
    if(avo_heat_melt$value[x]<=0.01){"<=0.01"}else{
      if(avo_heat_melt$value[x]<=0.05){"<=0.05"}else{">0.05"}}})
  
  heat<-ggplot(avo_heat_melt, aes(variable,meta, fill = Pvalue)) +
    geom_tile(color = "black",size=0.5) +
    theme_gray(8)+scale_fill_manual(values=c("#084594","#4292c6","#9ecae1","#deebf7"))+
    theme(axis.text = element_text(size =10, color="black"),
          axis.text.x = element_text(),
          axis.title = element_text(size =15),
          legend.text = element_text(size =14),
          legend.title = element_text(size =12),
          legend.position = c(1, 0), legend.justification = c(1,0),
          plot.margin=unit(c(0,2.25,1,1),"cm"))+
    xlab("Adjusted Principal Component")+ylab(NULL)
  
  grid.arrange(scree, heat, ncol=1, heights = c(2, 4))
}




#SWAN
## generated 20120703 temporary since bg is faked instead of from real data
## preprocessMSWAN - adapted preprocessSWAN for methylumiM

preprocessMSWAN <- function(MethyLumiM,bg=NULL){
  if(class(MethyLumiM)!="MethyLumiM")stop("Object must be MethyLumiM class")
  typeI <- fData(MethyLumiM)[fData(MethyLumiM)$INFINIUM_DESIGN_TYPE=="I",c("NAME","nCpG")]
  typeII <- fData(MethyLumiM)[fData(MethyLumiM)$INFINIUM_DESIGN_TYPE=="II",c("NAME","nCpG")]
  subset <- min(table(typeI$nCpG[typeI$nCpG <= 3 & typeI$nCpG > 0]),
                table(typeII$nCpG[typeII$nCpG <= 3 & typeII$nCpG >0]))
  CpG.counts <- rbind(typeI,typeII)
  CpG.counts$NAME <- as.character(CpG.counts$NAME)
  CpG.counts$Type <- rep(c("I", "II"), times = c(nrow(typeI), nrow(typeII)))
  names(CpG.counts)[1:2] <- c("Name","CpGs")
  counts <- CpG.counts[CpG.counts$Name %in% featureNames(MethyLumiM), ]
  ######### fake generated background!
  if(is.null(bg))bg <- rep(200.00,ncol(MethyLumiM))
  bg <- bg##    bg <- bgIntensitySwan(
  #########
  methData <- methylated(MethyLumiM)
  unmethData <- unmethylated(MethyLumiM)
  normMethData <- NULL
  normUnmethData <- NULL
  xNormSet <- vector("list", 2)
  xNormSet[[1]] <- getSubset(counts$CpGs[counts$Type == "I"],
                             subset)
  xNormSet[[2]] <- getSubset(counts$CpGs[counts$Type == "II"],
                             subset)
  for (i in 1:ncol(MethyLumiM)) {
    cat(sprintf("Normalizing array %d of %d\n", i, ncol(MethyLumiM)))
    normMethData <- cbind(normMethData,
                          normaliseChannel(methData[rownames(methData) %in% counts$Name[counts$Type == "I"], i],
                                           methData[rownames(methData) %in% counts$Name[counts$Type == "II"], i],
                                           xNormSet, bg[i]))
    normUnmethData <- cbind(normUnmethData,
                            normaliseChannel(unmethData[rownames(unmethData) %in%
                                                          counts$Name[counts$Type == "I"], i],
                                             unmethData[rownames(unmethData) %in%
                                                          counts$Name[counts$Type == "II"], i],
                                             xNormSet, bg[i]))
  }
  colnames(normMethData) <- sampleNames(MethyLumiM)
  colnames(normUnmethData) <- sampleNames(MethyLumiM)
  normMethyLumiM <- MethyLumiM
  stopifnot(all(featureNames(normMethyLumiM)%in%rownames(normMethData)))
  assayDataElement(normMethyLumiM, "methylated") <- normMethData[featureNames(normMethyLumiM),]
  stopifnot(all(featureNames(normMethyLumiM)%in%rownames(normUnmethData)))
  assayDataElement(normMethyLumiM, "unmethylated") <- normUnmethData[featureNames(normMethyLumiM),]
  normMethyLumiM@annotation <- c(sprintf("SWAN (based on a MethylSet preprocesses as in minfi)",
                                         MethyLumiM@annotation[1]))
  normMethyLumiM <- estimateM(normMethyLumiM)
  normMethyLumiM
}

subsetQuantileNorm <- function (x, xNormSet, xTarget, bg)
{
  for (i in 1:length(x)) {
    n <- length(x[[i]])
    nTarget <- length(xTarget)
    nNormSet <- sum(xNormSet[[i]])
    if (nNormSet != nTarget) {
      targetQuantiles <- (0:(nTarget - 1))/(nTarget - 1)
      r <- rank(x[xNormSet[, i], i])
      xNew <- (r - 1)/(nNormSet - 1)
      xNew <- xNew[order(xNew)]
      xNorm <- approx(x = targetQuantiles, y = xTarget,
                      xout = xNew, ties = "ordered", rule = 2)$y
    }
    else {
      xNorm <- xTarget
    }
    r <- rank(x[[i]])
    xNew <- (r - 1)/(n - 1)
    quantiles <- xNew[xNormSet[[i]]]
    quantiles <- quantiles[order(quantiles)]
    xmin <- min(x[[i]][xNormSet[[i]]])
    xmax <- max(x[[i]][xNormSet[[i]]])
    kmax <- which(xNew > max(quantiles))
    kmin <- which(xNew < min(quantiles))
    offsets.max <- x[[i]][kmax] - xmax
    offsets.min <- x[[i]][kmin] - xmin
    x[[i]] <- approx(x = quantiles, y = xNorm, xout = xNew,
                     ties = "ordered")$y
    x[[i]][kmax] <- max(xNorm) + offsets.max
    x[[i]][kmin] <- min(xNorm) + offsets.min
    x[[i]] = ifelse(x[[i]] <= 0, bg, x[[i]])
  }
  x
}

aveQuantile <- function (X)
{
  nbrOfChannels <- length(X)
  if (nbrOfChannels == 1) {
    return(X)
  }
  nbrOfObservations <- unlist(lapply(X, FUN = length), use.names = FALSE)
  maxNbrOfObservations <- max(nbrOfObservations)
  if (maxNbrOfObservations == 1) {
    return(X)
  }
  quantiles <- (0:(maxNbrOfObservations - 1))/(maxNbrOfObservations -
                                                 1)
  xTarget <- vector("double", maxNbrOfObservations)
  for (cc in 1:nbrOfChannels) {
    Xcc <- X[[cc]]
    Scc <- sort(Xcc)
    nobs <- length(Scc)
    if (nobs < maxNbrOfObservations) {
      bins <- (0:(nobs - 1))/(nobs - 1)
      Scc <- approx(x = bins, y = Scc, xout = quantiles,
                    ties = "ordered")$y
    }
    xTarget <- xTarget + Scc
  }
  rm(Scc, Xcc)
  xTarget <- xTarget/nbrOfChannels
  xTarget
}

normaliseChannel <- function (intensityI, intensityII, xNormSet, bg)
{
  xTarget <- aveQuantile(list(intensityI[xNormSet[[1]]], intensityII[xNormSet[[2]]]))
  xNorm <- unlist(subsetQuantileNorm(list(intensityI, intensityII),
                                     xNormSet, xTarget, bg))
  names(xNorm) <- names(c(intensityI, intensityII))
  xNorm
}

getSubset <- function (counts, subset)
{
  x <- numeric(0)
  for (i in 1:3) {
    x <- c(x, sample(seq(1, length(counts), by = 1)[counts ==
                                                      i], subset))
  }
  return(seq(1, length(counts)) %in% x)
}

## ## Original scripts from minfi package
## bgIntensitySwan <- function (rgSet)
## {
##     grnMed <- colMedians(getGreen(rgSet)[getControlAddress(rgSet,
##         controlType = "NEGATIVE"), ])
##     redMed <- colMedians(getRed(rgSet)[getControlAddress(rgSet,
##         controlType = "NEGATIVE"), ])
##     return(rowMeans(cbind(grnMed, redMed)))
## }

## ## Original scripts from minfi package
## preprocessSWAN <- function (rgSet, mSet = NULL)
## {
##     if (is.null(mSet))
##         mSet <- preprocessRaw(rgSet)
##     typeI <- getProbeInfo(rgSet, type = "I")[, c("Name", "nCpG")]
##     typeII <- getProbeInfo(rgSet, type = "II")[, c("Name", "nCpG")]
##     subset <- min(table(typeI$nCpG[typeI$nCpG <= 3 & typeI$nCpG >
##         0]), table(typeII$nCpG[typeII$nCpG <= 3 & typeII$nCpG >
##         0]))
##     CpG.counts <- rbind(typeI, typeII)
##     CpG.counts$Name <- as.character(CpG.counts$Name)
##     CpG.counts$Type <- rep(c("I", "II"), times = c(nrow(typeI),
##         nrow(typeII)))
##     names(CpG.counts)[2] <- "CpGs"
##     counts <- CpG.counts[CpG.counts$Name %in% featureNames(mSet),
##         ]
##     bg <- bgIntensitySwan(rgSet)
##     methData <- getMeth(mSet)
##     unmethData <- getUnmeth(mSet)
##     normMethData <- NULL
##     normUnmethData <- NULL
##     xNormSet <- vector("list", 2)
##     xNormSet[[1]] <- getSubset(counts$CpGs[counts$Type == "I"],
##         subset)
##     xNormSet[[2]] <- getSubset(counts$CpGs[counts$Type == "II"],
##         subset)
##     for (i in 1:ncol(mSet)) {
##         cat(sprintf("Normalizing array %d of %d\n", i, ncol(mSet)))
##         normMethData <- cbind(normMethData, normaliseChannel(methData[rownames(methData) %in%
##             counts$Name[counts$Type == "I"], i], methData[rownames(methData) %in%
##             counts$Name[counts$Type == "II"], i], xNormSet, bg[i]))
##         normUnmethData <- cbind(normUnmethData, normaliseChannel(unmethData[rownames(unmethData) %in%
##             counts$Name[counts$Type == "I"], i], unmethData[rownames(unmethData) %in%
##             counts$Name[counts$Type == "II"], i], xNormSet, bg[i]))
##     }
##     colnames(normMethData) <- sampleNames(mSet)
##     colnames(normUnmethData) <- sampleNames(mSet)
##     normSet <- mSet
##     assayDataElement(normSet, "Meth") <- normMethData
##     assayDataElement(normSet, "Unmeth") <- normUnmethData
##     normSet@preprocessMethod <- c(sprintf("SWAN (based on a MethylSet preprocesses as '%s'",
##         mSet@preprocessMethod[1]), as.character(packageVersion("minfi")),
##         as.character(packageVersion("IlluminaHumanMethylation450kmanifest")))
##     normSet
## }



ECC5 <- function (RGset, meanPlot = TRUE, cellType = "Blood", verbose = TRUE, 
                  ...) 
{
  require(quadprog)
  library(genefilter)
  library(matrixStats)
  platform <- "450k"
  referencePkg <- sprintf("FlowSorted.%s.%s", "Blood", platform)
  if (!require(referencePkg, character.only = TRUE)) 
    stop(sprintf("Could not find reference data package for cellType '%s' and platform '%s' (inferred package name is '%s')", 
                 cellType, platform, referencePkg))
  data(list = referencePkg)
  referenceRGset <- get(referencePkg)
  if (verbose) 
    cat("[estimateCellCounts] Combining Data with Flow Sorted Data.\n")
  #RGsetComb <- combine(RGset, referenceRGset)
  RGsetComb <-   BiocGenerics::combine(RGset, as.data.frame(getBeta(referenceRGset))) #something aint right
  newpd <- data.frame(sampleNames = c(names(RGset), sampleNames(referenceRGset)), 
                      studyIndex = rep(c("user", "reference"), times = c(ncol(RGset), 
                                                                         ncol(referenceRGset))), stringsAsFactors = FALSE)
  #pData(RGsetComb) <- newpd
  if (verbose) 
    cat("[estimateCellCounts] Normalizing Data Together.\n")
  #Mset <- preprocessQuantile(RGsetComb, removeBadSamples = FALSE,...)
  #referenceMset <- Mset[, Mset$studyIndex == "reference"]
  library(limma)
  Mset <- normalizeQuantiles(RGsetComb)
  referenceMset <- Mset[, newpd[,2] == "reference"]
  refpd <- cbind(newpd[newpd[,2] == "reference",],referenceRGset$CellType)
  #pData(referenceMset) <- as(pData(referenceRGset), "DataFrame")
  #Mset <- Mset[, Mset$studyIndex == "user"]
  Mset <- Mset[, newpd[,2] == "user"]
  #pData(Mset) <- as(pData(RGset), "DataFrame")
  if (verbose) 
    cat("[estimateCellCounts] Picking Probes for Composition Estimation.\n")
  #compData <- pickCompProbes(referenceMset)
  compData <- pickCompProbes2(referenceMset,refpd)
  coefs <- compData$coefEsts
  if (verbose) 
    cat("[estimateCellCounts] Estimating Composition.\n")
  #counts <- projectWBC(getBeta(Mset)[rownames(coefs), ], coefs)
  counts <- projectWBC(Mset[rownames(coefs), ], coefs)
  #rownames(counts) <- names(RGset)
  print(table(rownames(counts) == names(RGset)))
  if (meanPlot) {
    smeans <- compData$sampleMeans
    smeans <- smeans[order(names(smeans))]
    #sampleMeans <- c(colMeans(getBeta(Mset)[rownames(coefs),]), smeans)
    sampleMeans <- c(colMeans(Mset[rownames(coefs),]), smeans)
    sampleColors <- c(rep(1, ncol(Mset)), 1 + as.numeric(factor(names(smeans))))
    plot(sampleMeans, pch = 21, bg = sampleColors)
    legend("bottomleft", c("blood", levels(factor(names(smeans)))), 
           col = 1:7, pch = 15)
  }
  list(counts = counts, compTable = compData$compTable, sortedData = referenceMset)
}


pickCompProbes2 <- function(Mset, pdref, numProbes = 50) {
  splitit <- function(x) {
    split(seq(along=x), x)
  }
  
  #p <- getBeta(Mset)
  p <- as.matrix(Mset)
  #pd <- as.data.frame(pData(Mset))
  pd <- as.data.frame(pdref)
  
  ## only keep 6 components from kere
  keep <- which(pd[,3] %in% c("Mono", "Bcell", 
                              "Gran", "CD4T", "CD8T", "NK"))
  pd <- pd[keep,]
  p <- p[,keep]
  ?rowRanges
  ## make cell type a factor 
  pd$CellType <- factor(pd[,3], 
                        levels = c("CD8T","CD4T", "NK","Bcell","Mono","Gran"))
  
  # get fstats
  ffComp <- rowFtests(p, pd$CellType)
  prof <- sapply(splitit(pd$CellType), function(i) rowMeans(p[,i]))
  r <- cbind(apply(p, 1, min), apply(p, 1, max))
  compTable <- cbind(ffComp, prof, r, abs(r[,1] - r[,2]))
  
  
  names(compTable)[c(1, 9:11)] <- c("Fstat", "low", "high", "range")
  
  # t-test by cell type
  tIndexes <- splitit(pd$CellType)
  str(tIndexes)
  tstatList <- lapply(tIndexes, function(i) {
    x <- rep(0,ncol(p))
    x[i] <- 1
    return(rowttests(p, factor(x)))
  })
  
  ## take N up and N down
  probeList <- lapply(tstatList, function(x) {
    y <- x[x[,"p.value"] < 1e-8,]
    yUp <- y[order(y[,"dm"], decreasing=TRUE),]
    yDown <- y[order(y[,"dm"], decreasing=FALSE),]
    c(rownames(yUp)[1:numProbes], rownames(yDown)[1:numProbes])
  })
  write.csv(probeList, file="probes_comp.csv")
  trainingProbes <- unlist(probeList)
  p <- p[trainingProbes,]
  
  
  pMeans <- colMeans(p)
  names(pMeans) <- pd$CellType
  
  mod <- model.matrix(~pd$CellType-1)
  colnames(mod) <- levels(pd$CellType)
  #form <- as.formula(sprintf("y ~ %s - 1", colnames(mod), collapse="+"))
  form <- as.formula("y ~ CD8T + CD4T + NK + Bcell + Mono + Gran -1")
  
  tmp <- validationWBC(p,data.frame(mod),form)
  coefEsts <- tmp$coefEsts
  
  out <- list(coefEsts = coefEsts, compTable = compTable,
              sampleMeans = pMeans)
  return(out)
}

projectWBC <- function(Y, coefWBC, contrastWBC=NULL, nonnegative=TRUE, lessThanOne=FALSE){ 
  if(is.null(contrastWBC))
    Xmat <- coefWBC
  else
    Xmat <- coefWBC %*% t(contrastWBC) 
  
  nCol <- dim(Xmat)[2]
  nSubj <- dim(Y)[2]
  
  mixCoef <- matrix(0, nSubj, nCol)
  rownames(mixCoef) <- colnames(Y)
  colnames(mixCoef) <- colnames(Xmat)
  
  if(nonnegative){
    library(quadprog)
    
    if(lessThanOne){
      Amat <- cbind(rep(-1,nCol), diag(nCol))
      b0vec <- c(-1,rep(0,nCol))
    } else {
      Amat <- diag(nCol)
      b0vec <- rep(0,nCol)
    }
    
    for(i in 1:nSubj) {
      obs <- which(!is.na(Y[,i])) 
      Dmat <- t(Xmat[obs,]) %*% Xmat[obs,]
      mixCoef[i,] <- solve.QP(Dmat, t(Xmat[obs,]) %*% Y[obs,i], Amat, b0vec)$sol
    }
  } else {
    for(i in 1:nSubj) {
      obs <- which(!is.na(Y[,i])) 
      Dmat <- t(Xmat[obs,]) %*% Xmat[obs,]
      mixCoef[i,] <- solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
    }
  }
  return(mixCoef)
}

validationWBC <- function(Y, pheno, modelFix, modelBatch=NULL, L.forFstat = NULL){
  N <- dim(pheno)[1]
  pheno$y <- rep(0, N)
  xTest <- model.matrix(modelFix, pheno)
  sizeModel <- dim(xTest)[2]
  M <- dim(Y)[1]
  
  if(is.null(L.forFstat)) {
    L.forFstat <- diag(sizeModel)[-1,]  #All non-intercept coefficients
    colnames(L.forFstat) <- colnames(xTest) 
    rownames(L.forFstat) <- colnames(xTest)[-1] 
  }
  
  ## Initialize various containers
  sigmaResid <- sigmaIcept <- nObserved <- nClusters <- Fstat <- rep(NA, M)
  coefEsts <- matrix(NA, M, sizeModel)
  coefVcovs <- list()
  
  for(j in 1:M) { # For each CpG
    ## Remove missing methylation values
    ii <- !is.na(Y[j,])
    nObserved[j] <- sum(ii)
    pheno$y <- Y[j,]
    
    if(j%%round(M/10)==0)
      cat(j,"\n") # Report progress
    
    try({ # Try to fit a mixed model to adjust for plate
      if(!is.null(modelBatch)) {
        fit <- try(lme(modelFix, random=modelBatch, data=pheno[ii,]))
        OLS <- inherits(fit,"try-error") # If LME can't be fit, just use OLS
      } else
        OLS <- TRUE
      
      if(OLS) {
        fit <- lm(modelFix, data=pheno[ii,])
        fitCoef <- fit$coef
        sigmaResid[j] <- summary(fit)$sigma
        sigmaIcept[j] <- 0
        nClusters[j] <- 0
      } else { 
        fitCoef <- fit$coef$fixed
        sigmaResid[j] <- fit$sigma
        sigmaIcept[j] <- sqrt(getVarCov(fit)[1])
        nClusters[j] <- length(fit$coef$random[[1]])
      }
      coefEsts[j,] <- fitCoef
      coefVcovs[[j]] <- vcov(fit)
      
      useCoef <- L.forFstat %*% fitCoef
      useV <- L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
      Fstat[j] <- (t(useCoef) %*% solve(useV, useCoef))/sizeModel
    })
  }
  ## Name the rows so that they can be easily matched to the target data set
  rownames(coefEsts) <- rownames(Y)
  colnames(coefEsts) <- names(fitCoef)
  degFree <- nObserved - nClusters - sizeModel + 1
  
  ## Get P values corresponding to F statistics
  Pval <- 1-pf(Fstat, sizeModel, degFree)
  
  out <- list(coefEsts=coefEsts, coefVcovs=coefVcovs, modelFix=modelFix, modelBatch=modelBatch,
              sigmaIcept=sigmaIcept, sigmaResid=sigmaResid, L.forFstat=L.forFstat, Pval=Pval,
              orderFstat=order(-Fstat), Fstat=Fstat, nClusters=nClusters, nObserved=nObserved,
              degFree=degFree)
  
  out
}


