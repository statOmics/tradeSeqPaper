---
title: 'Case study: OE dataset - three versions of DE'
author: "Koen Van den Berge"
date: "`r Sys.Date()`"
---

# pca data

```{r loadPCAData}
load("~/PhD_Data/singleCell/fletcher/ShareWithKelly/E4c2b_slingshot_wsforkelly.RData")
library(slingshot)
library(rgl)
library(rafalib) ; mypar()
library(RColorBrewer)
library(mgcv)
library(tradeR)

origData <- t(pcax$x %*% t(pcax$rotation))
colpal=cc
rd <- X[,1:5]
lin <- getLineages(rd, clusterLabels = clus.labels, start.clus = "1", end.clus="4")
crv <- getCurves(lin)
plot(X[,1:2], col=colpal[as.factor(clus.labels)], pch=16)
lines(crv)
#3D plot
# rgl::plot3d(X[,1:3], t='p', col=colpal[as.factor(clus.labels)],alpha=0.3, pch = 19, cex = 2, size=8, xlab="PC 1", ylab="PC 2", zlab="PC 3", aspect="iso", box=FALSE, axes=FALSE)
# #rgl::plot3d(X[,1:3], t='p', col=c("white","black")[hlp+1],alpha=0.3, pch = 19, cex = 2, size=8, xlab="PC 1", ylab="PC 2", zlab="PC 3", aspect="iso", box=FALSE, axes=FALSE)
# rgl::axes3d(tick=FALSE)
# rgl::par3d(windowRect = c(20, 30, 800, 800))
# for (i in seq_along(curves)){
#   rgl::plot3d(crv@curves[[i]]$s[order(crv@curves[[i]]$lambda),1:3], type='l',add=TRUE, lwd=4,col=colpal[which.max(tail(lin@lineages[[i]],1)==levels(clus.labels))])
# }
# rgl.postscript("~/fletcher3d.pdf", fmt="pdf")

```

In the plot below, trajectory 1 (neuronal trajectory) is the green curve, trajectory 2 is the yellow curve, and trajectory 3 (sustentacular trajectory) is the brown curve.

```{r plot1, echo=TRUE, fig.cap="Fletcher 3D PCA", include=identical(knitr:::pandoc_to(), 'html')}
knitr::include_graphics("/Users/koenvandenberge/fletcher3d.pdf")
#![Fletcher 3D PCA](/Users/koenvandenberge/fletcher3d.pdf)
```

```{r cufflinksCountData}
suppressPackageStartupMessages(library(SummarizedExperiment))
load("~/Downloads/GSE95601_oeHBCdiff_Cufflinks_eSet.rda")
counts <- assayData(Cufflinks_eSet)$counts_table
counts = counts[!apply(counts,1,function(row) any(is.na(row))),]
counts <- counts[,colnames(counts)%in%colnames(origData)]
keep <- rowSums(edgeR::cpm(counts)>5)>=15
countsFiltered <- counts[keep,]
```

```{r EDA}
# pdf("~/batchOE.pdf")
# for(i in 1:24) plotMDS(d, pch=16, col=droplevels(phenoData(Cufflinks_eSet)[[i]]))
# dev.off()

# batch from Github
library(scone) ; library(edgeR)
load("~/p63-HBC-diff/output/clust/oeHBCdiff/oeHBCdiff_scone.Rda")
batch <- droplevels(colData(scone_out)$batch[match(colnames(countsFiltered),rownames(colData(scone_out)))])
run <- droplevels(phenoData(Cufflinks_eSet)[[3]])[colnames(exprs(Cufflinks_eSet))%in%colnames(origData)]
#table(batch,run)
```

# fit smoothers on unaligned data

## get ZI weights from zinbwave

```{r zinbwave}
library(zinbwave)
library(BiocParallel)
library(doParallel)
NCORES <- 2
registerDoParallel(NCORES)
register(DoparParam())

core <- SummarizedExperiment(countsFiltered, colData = data.frame(clusLabel = clus.labels, batch=batch))
#zinb_c <- zinbFit(core, X = '~ clusLabel + batch', commondispersion = TRUE)
#save(zinb_c,file="~/zinbFletcher.rda")
load("~/zinbFletcher.rda")
weights <- computeObservationalWeights(zinb_c, countsFiltered)
```


```{r fitGAM}
#gamList <- tradeR::fitGAM(countsFiltered, X=model.matrix(~batch), pseudotime=slingPseudotime(crv, na=FALSE), cellWeights=slingCurveWeights(crv), weights=weights)
#save(gamList,file="~/gamListOE_tradeR.rda")
load("~/gamListOE_tradeR.rda")
```


#####################################
### BETWEEN LINEAGE COMPARISONS #####
#####################################

# Comparing the expression at the differentiated cell type

```{r}
### ERROR: model 2548 gives negative variance for the lineage 1 - lineage 2 contrast. how come?

endTestGam <- endPointTest(gamList, omnibus=TRUE, pairwise=TRUE)
endOmnibusPval <- endTestGam$pvalue
sum(is.na(endOmnibusPval)) #genes we could not fit or test.
endOmnibusPadj <- p.adjust(endOmnibusPval,"fdr")
sum(endOmnibusPadj<=0.05, na.rm=TRUE) ; mean(endOmnibusPadj<=0.05, na.rm=TRUE)
deGenesEndGam <- which(endOmnibusPadj<=0.05)
hist(endOmnibusPval)
```

# edgeR analysis

In the edgeR analysis, I compare the gene expression of the differentiated cell types (i.e. final cluster of a trajectory) to each other, between the different trajectories.

```{r}
library(edgeR)
d <- DGEList(countsFiltered)
d <- calcNormFactors(d)
design <- model.matrix(~clus.labels+batch)
d$weights <- weights
d <- estimateDisp(d,design)
fit <- glmFit(d,design)
L <- matrix(0,nrow=ncol(fit$coefficients),ncol=3)
rownames(L) <- colnames(fit$coefficients)
# trajectory 1 vs. 2
L[c("clus.labels12","clus.labels15"),1] <- c(1,-1)
# trajectory 3 vs 1
L[c("clus.labels4","clus.labels12"),2] <- c(1,-1)
# trajectory 3 vs 2
L[c("clus.labels4","clus.labels15"),3] <- c(1,-1)
lrt <- zinbwave::glmWeightedF(fit,contrast=L)
deGenesEdgeR <- which(p.adjust(lrt$table$PValue,"fdr")<=0.05)
length(deGenesEdgeR)
mean(deGenesEdgeR%in%deGenesEndGam)
```

# compare analyses

We retrieve 85% of the genes that edgeR finds, but also obtain a bunch of other genes.

```{r}
mypar(mfrow=c(1,1))

## compare
edgeRDE <- p.adjust(lrt$table$PValue,"fdr")<=0.05
gamDE <- endOmnibusPadj<=0.05
vennC <- cbind(edgeRDE,gamDE)
vennDiagram(vennC, main="end point comparison across all trajectories")
```

# plot unique GAM genes

```{r}
uniqueGamEndId <- which(gamDE==TRUE & edgeRDE==FALSE)

i=0
while(i<10){
  i=i+1
  plotSmoothers(gamList[[uniqueGamEndId[i]]])
}
```

# plot unique edgeR genes

```{r}
uniqueEdgeREndId <- which(gamDE==FALSE & edgeRDE==TRUE)

i=0
while(i<10){
  i=i+1
  plotSmoothers(gamList[[uniqueEdgeREndId[i]]])
}
```



```{r patternTest}
resPat <- patternTest(gamList, omnibus=TRUE, pairwise=TRUE)
o <- order(resPat$waldStat, decreasing=TRUE)

i=0
while(i<10){
  i=i+1
  plotSmoothers(gamList[[o[i]]], main=i)
}

# for paper
mypar(mfrow=c(3,2))
i=0
while(i<6){
  i=i+1
  plotSmoothers(gamList[[o[i]]], main=rownames(resPat)[o][i])
}

# Stage-wise testing
library(stageR)
pScreen <-  resPat$pvalue
names(pScreen) <- rownames(resPat)
pConfirmation <- cbind(resPat$pvalue_1vs2, resPat$pvalue_1vs3, resPat$pvalue_2vs3)
dimnames(pConfirmation) <- list(rownames(resPat),c("1v2","1v3","2v3"))
stageObj <- stageR(pScreen, pConfirmation, pScreenAdjusted=FALSE)
stageObj <- stageWiseAdjustment(stageObj, alpha=.05, method="holm", allowNA=TRUE)
res <- getResults(stageObj)
colSums(res)
allCompGenes <- names(which(rowSums(res)==4))
write.table(allCompGenes, file="~/allCompGenes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
# submit these in http://software.broadinstitute.org/gsea/msigdb/compute_overlaps.jsp

 overlap=readLines("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/overlapCellCycleOE")
overlapSets <- overlap[10:30]
overlapSetsSplit <- sapply(overlapSets, function(x) strsplit(x,split="\t"))
gsNames <- unname(unlist(lapply(overlapSetsSplit,"[[",1)))
genesInSet <- unname(unlist(lapply(overlapSetsSplit,"[[",2)))
genesInOverlap <- unname(unlist(lapply(overlapSetsSplit,"[[",4)))
qval <- unname(unlist(lapply(overlapSetsSplit,"[[",7)))
gsNames <- tolower(gsNames)
gsNames <- gsub(x=gsNames,pattern="_",replacement=" ")
gsNames[-1] <- unname(sapply(gsNames[-1],function(x) substr(x,4,nchar(x))))
tab <- data.frame(geneSet=gsNames[-1],
                    overlap=genesInOverlap[-1],
                  genesInSet=genesInSet[-1],
                  qvalue=qval[-1])
#library(xtable)
#xtable(tab)

# cluster fitted values of genes significant for pattern test and create a heatmap.
# or plot gene clusters their smoothing profiles
log(fitted(gamList[[1]])+.1)

```





###### within lineage comparisons

```{r cellCycleWithin}
plot(x=levels(clus.labels), pch=16, col=colpal)
# colors 5,7,4 are last three clusters of SUS lineage
# colors 2,14,10,9,11 are last 5 clusters of NEUR lineage
# colors 1, 8, 5, 3, 2,14,10,9,11 are all clusters of NEUR lineage

L <- matrix(0,nrow=ncol(fit$coefficients),ncol=9)
rownames(L) <- colnames(fit$coefficients)
colnames(L) <- crv@lineages[[1]]
L[c("clus.labels8","clus.labels5","clus.labels3","clus.labels2","clus.labels14","clus.labels10","clus.labels9","clus.labels11"),1] <- rep(-1/8,8)
L[c("clus.labels8","clus.labels5","clus.labels3","clus.labels2","clus.labels14","clus.labels10","clus.labels9","clus.labels11"),2] <- c(1,rep(-1/8,7))
L[c("clus.labels8","clus.labels5","clus.labels3","clus.labels2","clus.labels14","clus.labels10","clus.labels9","clus.labels11"),3] <- c(-1,1,rep(-1/8,6))
L[c("clus.labels8","clus.labels5","clus.labels3","clus.labels2","clus.labels14","clus.labels10","clus.labels9","clus.labels11"),4] <- c(rep(-1/8,2),1,rep(-1/8,5))
L[c("clus.labels8","clus.labels5","clus.labels3","clus.labels2","clus.labels14","clus.labels10","clus.labels9","clus.labels11"),5] <- c(rep(-1/8,3),1,rep(-1/8,4))
L[c("clus.labels8","clus.labels5","clus.labels3","clus.labels2","clus.labels14","clus.labels10","clus.labels9","clus.labels11"),6] <- c(rep(-1/8,4),1,rep(-1/8,3))
L[c("clus.labels8","clus.labels5","clus.labels3","clus.labels2","clus.labels14","clus.labels10","clus.labels9","clus.labels11"),7] <- c(rep(-1/8,5),1,rep(-1/8,2))
L[c("clus.labels8","clus.labels5","clus.labels3","clus.labels2","clus.labels14","clus.labels10","clus.labels9","clus.labels11"),8] <- c(rep(-1/8,6),1,rep(-1/8,1))
L[c("clus.labels8","clus.labels5","clus.labels3","clus.labels2","clus.labels14","clus.labels10","clus.labels9","clus.labels11"),9] <- c(rep(-1/8,7),1)
lrtList <- list()
for(i in 1:ncol(L)) lrtList[[i]] <- zinbwave::glmWeightedF(fit,contrast=L[,i])


cellGenes=read.table("~/Downloads/GO_term_summary_20181001_091947.txt", header=TRUE, sep="\t", row.names=NULL)
cellGenes$MGI.Gene.Marker.ID <- as.character(cellGenes$MGI.Gene.Marker.ID)
mean(cellGenes$MGI.Gene.Marker.ID%in%rownames(counts))
mean(cellGenes$MGI.Gene.Marker.ID%in%rownames(countsFiltered))
sum(cellGenes$MGI.Gene.Marker.ID%in%rownames(countsFiltered))

sumFoundMat <- matrix(NA, nrow=nrow(countsFiltered), ncol=9)
for(i in 1:length(lrtList)){
  pvals <- lrtList[[i]]$table$PValue
  o <- order(pvals)
  pvalsO <- pvals[o]
  namesO <- rownames(lrtList[[i]])[o]
  sumFound <- sapply(1:nrow(countsFiltered), function(ii){
    sum(cellGenes$MGI.Gene.Marker.ID %in% namesO[1:ii])
  })
  sumFoundMat[,i] <- sumFound
}

# aggregate cluster comparisons
library(aggregation)
pvalMatEdgeR <- do.call(cbind,lapply(lrtList,function(x) x$table$PValue))
fisherP <- apply(pvalMatEdgeR,1,fisher)
names(fisherP) <- rownames(lrtList[[1]]$table)
o <- order(fisherP)
sumFoundFisher <- sapply(1:nrow(countsFiltered), function(ii){
  sum(cellGenes$MGI.Gene.Marker.ID %in% names(fisherP[o])[1:ii])
})

smootherStats <- getSmootherTestStats(gamList)
o <- order(smootherStats[,1], decreasing=TRUE)
sumFoundPat <- sapply(1:nrow(countsFiltered), function(ii){
  sum(cellGenes$MGI.Gene.Marker.ID %in% rownames(smootherStats[o,])[1:ii])
})


#in top 1000 genes we find 394 cell cycle genes
sum(cellGenes$MGI.Gene.Marker.ID %in% rownames(smootherStats[o,])[1:1000])
# out of a total of 2889 genes present in filtered dataset
sum(cellGenes$MGI.Gene.Marker.ID %in% rownames(smootherStats))

# make enrichment plot
# absolute
plot(x=1:nrow(smootherStats), y=sumFoundPat, type="l")
lines(x=1:nrow(smootherStats), y=(2889/nrow(smootherStats))*(1:nrow(smootherStats)), type="l", col="steelblue")
expRandom <- (2889/nrow(smootherStats))*(1:nrow(smootherStats))
# relative
plot(x=1:nrow(smootherStats), y=sumFoundPat/expRandom, type="l", ylab="# genes found with tradeR / # genes found under random selection")
abline(h=1, lty=2, col="red")



#library(scales)
#plot(x=1:nrow(countsFiltered), y=sumFoundPat, type="l", log="y", col="darkseagreen3", lwd=2)
#lines(x=1:nrow(countsFiltered), y=sumFoundFisher, col="red", lwd=2)
#for(i in 1:9) lines(x=1:nrow(countsFiltered), y=sumFoundMat[,i], col=alpha(colorRampPalette(c("lightblue","darkblue"))(9)[i],.5))
```

# Plot for paper


```{r}
plotSmoothersIk <- function(m, nPoints = 100, ...){

  data <- m$model
  y <- data$y

  #construct time variable based on cell assignments.
  nCurves <- length(m$smooth)
  timeAll <- c()
  col <- c()
  for (jj in seq_len(nCurves)) {
    for (ii in 1:nrow(data)) {
      if (data[ii, paste0("l", jj)] == 1) {
        timeAll[ii] <- data[ii, paste0("t", jj)]
        col[ii] <- jj
      } else {
        next
      }
    }
  }

  # plot raw data
    cols <- c("#B3DE69", "#FFED6F", "#A6761D")
  plot(x = timeAll, y = log(y + 1), col = cols[col], pch = 16, cex = 2 / 3,
       ylab = "Expression + 1 (log-scale)", xlab = "Pseudotime", ...)

  #predict and plot smoothers across the range
  for (jj in seq_len(nCurves)) {
    df <- .getPredictRangeDf(m, jj, nPoints = nPoints)
  yhat <- predict(m, newdata = df, type = "response")
  lines(x = df[, paste0("t", jj)], y = log(yhat + 1), col = cols[jj], lwd = 2)
  }
  legend("bottomleft", c("Neuronal", "Microvillous", "Sustentacular"),col = cols,
         lty = 1, lwd = 2, bty = "n", cex = 2 / 3)
}



#par(mfrow=c(1,2))
rgl::plot3d(X[,1:3], t='p', col=colpal[as.factor(clus.labels)],alpha=0.3, pch = 19, cex = 2, size=8, xlab="", ylab="", zlab="", aspect="iso", box=FALSE, axes=TRUE)
for (i in seq_along(curves)){
  rgl::plot3d(crv@curves[[i]]$s[order(crv@curves[[i]]$lambda),1:3], type='l',add=TRUE, lwd=4,col=colpal[which.max(tail(lin@lineages[[i]],1)==levels(clus.labels))])
}
#rgl::text3d(x=0,y=-100,z=0,"PC1")
rgl.postscript("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/plots/fletcher3d.pdf", fmt="pdf")

png("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/plots/figCaseFletcher2.png", width=9, height=6.5, units="in", res=200)
mypar(bty="l")
layout(matrix(c(1,1,2,3),nrow=2,ncol=2,byrow=TRUE))
# important gene plot
plotSmoothersIk(gamList[[o[6]]])
mtext("b", at=-30, font=2, cex=2)

# barplot DE
nrDEPat <- colSums(res)[-1]
barplot(nrDEPat, ylab="Frequency", names=c("Neur vs. Micro", "Neur vs. Sus", "Micro vs. Sus"))
mtext("c", at=-.2, font=2, cex=2)

# relative enrichment plot
par(mar=c(2.5,4,1.6,1.1), bty="l")
plot(x=1:nrow(smootherStats), y=sumFoundPat/expRandom, type="l", ylab="Ratio of cell cycle genes found \n with tradeR vs. random selection", xlab="Rank")
abline(h=1, lty=2, col="red")
mtext("d", at=-2300, font=2, cex=2)
dev.off()


```



# Check for HBC markers

```{r}
startRes <- startPointTest(gamList)
top250Genes <- rownames(head(startRes[order(startRes$waldStat, decreasing=TRUE),],250))
write.table(top250Genes, file="~/startTestOETop250.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
# submit to http://software.broadinstitute.org/gsea/msigdb/compute_overlaps.jsp for top 20 gene sets


overlap=readLines("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/overlapStartTestOE")
overlapSets <- overlap[10:30]
overlapSetsSplit <- sapply(overlapSets, function(x) strsplit(x,split="\t"))
gsNames <- unname(unlist(lapply(overlapSetsSplit,"[[",1)))
genesInSet <- unname(unlist(lapply(overlapSetsSplit,"[[",2)))
genesInOverlap <- unname(unlist(lapply(overlapSetsSplit,"[[",4)))
qval <- unname(unlist(lapply(overlapSetsSplit,"[[",7)))
gsNames <- tolower(gsNames)
gsNames <- gsub(x=gsNames,pattern="_",replacement=" ")
gsNames[-1] <- unname(sapply(gsNames[-1],function(x) substr(x,4,nchar(x))))
tab <- data.frame(geneSet=gsNames[-1],
                    overlap=genesInOverlap[-1],
                  genesInSet=genesInSet[-1],
                  qvalue=qval[-1])
#library(xtable)
#xtable(tab)
```







# Cell cycle functional annotation


# pattern contrast between neuronal and sustentacular lineage should reveal cell cycle genes, since the cells in neuronal lineage undergo cell division, while the sustentacular lineage cells do not. Also this occurs in intermediary stages, and not at end stage, so we need the pattern test for this.
# cell cycle gene set: http://www.informatics.jax.org/go/term/GO:0007049
# downloaded on October 1st 2018


```{r cc}
cellGenes=read.table("~/Downloads/GO_term_summary_20181001_091947.txt", header=TRUE, sep="\t", row.names=NULL)
cellGenes$MGI.Gene.Marker.ID <- as.character(cellGenes$MGI.Gene.Marker.ID)
mean(cellGenes$MGI.Gene.Marker.ID%in%rownames(counts))


# we expect the cell cycle genes to have different patterns but to not be DE at end point.
plot(x=log(endTestGam$waldStat_1vs3), y=log(resPat$waldStat_1vs3))
# note that Wald stats are not directly comparable due to df.
plot(x=endTestGam$waldStat_1vs3, y=resPat$waldStat_1vs3/resPat$df_1vs3,log="xy") ; abline(0,1,col="red")
# compare p-values:
logit=function(x) log(x/(1-x))
plot(x=logit(endTestGam$pvalue_1vs3), y=logit(resPat$pvalue1vs3), pch=16, cex=2/3)

# code from vignette, written by Hector:
library(dplyr)
compare <- inner_join(resPat %>% mutate(Gene = rownames(resPat),
                                            pattern = waldStat_1vs3) %>%
                                     select(Gene, pattern),
                      endTestGam %>% mutate(Gene = rownames(endTestGam),
                                        end = waldStat_1vs3) %>%
                                 select(Gene, end),
                      by = c("Gene" = "Gene")) %>%
           mutate(transientScore = (min_rank(desc(end)))^2 +
                                   (min_rank(pattern) - 1)^2)
library(ggplot2)
ggplot(compare, aes(x = log(pattern), y = log(end))) +
  geom_point(aes(col = transientScore)) +
  labs(x = "endPointTest Wald Statistic (log scale)",
       y = "patternTest Wald Statistic (log scale)") +
  scale_color_continuous(low = "yellow", high = "red") +
  theme_classic()

topTransient <- (compare %>% arrange(desc(transientScore)))



sumFoundTradeR <- sapply(1:nrow(countsFiltered), function(ii){
  sum(cellGenes$MGI.Gene.Marker.ID %in% topTransient$Gene[1:ii])
})


## get range of cluster comparisons to check if the tradeR analysis improves upon ANY cluster comparison
#in PCA plot, col=colpal[as.factor(clus.labels)]
plot(x=levels(clus.labels), pch=16, col=colpal)
# colors 5,7,4 are last three clusters of SUS lineage
# colors 2,14,10,9,11 are last 5 clusters of NEUR lineage
L <- matrix(0,nrow=ncol(fit$coefficients),ncol=15)
rownames(L) <- colnames(fit$coefficients)
L[c("clus.labels5","clus.labels2"),1] <- c(1,-1)
L[c("clus.labels5","clus.labels14"),2] <- c(1,-1)
L[c("clus.labels5","clus.labels10"),3] <- c(1,-1)
L[c("clus.labels5","clus.labels9"),4] <- c(1,-1)
L[c("clus.labels5","clus.labels11"),5] <- c(1,-1)
L[c("clus.labels7","clus.labels2"),6] <- c(1,-1)
L[c("clus.labels7","clus.labels14"),7] <- c(1,-1)
L[c("clus.labels7","clus.labels10"),8] <- c(1,-1)
L[c("clus.labels7","clus.labels9"),9] <- c(1,-1)
L[c("clus.labels7","clus.labels11"),10] <- c(1,-1)
L[c("clus.labels4","clus.labels2"),11] <- c(1,-1)
L[c("clus.labels4","clus.labels14"),12] <- c(1,-1)
L[c("clus.labels4","clus.labels10"),13] <- c(1,-1)
L[c("clus.labels4","clus.labels9"),14] <- c(1,-1)
L[c("clus.labels4","clus.labels11"),15] <- c(1,-1)
lrtList <- list()
for(i in 1:ncol(L)) lrtList[[i]] <- zinbwave::glmWeightedF(fit,contrast=L[,i])



sumFoundMat <- matrix(NA, nrow=nrow(countsFiltered), ncol=15)
for(i in 1:length(lrtList)){
  pvals <- lrtList[[i]]$table$PValue
  o <- order(pvals)
  pvalsO <- pvals[o]
  namesO <- rownames(lrtList[[i]])[o]
  sumFound <- sapply(1:nrow(countsFiltered), function(ii){
    sum(cellGenes$MGI.Gene.Marker.ID %in% namesO[1:ii])
  })
  sumFoundMat[,i] <- sumFound
}


o <- order(endTestGam$waldStat_1vs3, decreasing=TRUE)
sumFoundEnd <- sapply(1:nrow(countsFiltered), function(ii){
  sum(cellGenes$MGI.Gene.Marker.ID %in% rownames(endTestGam[o,])[1:ii])
})

o <- order(resPat$waldStat_1vs3, decreasing=TRUE)
sumFoundPat <- sapply(1:nrow(countsFiltered), function(ii){
  sum(cellGenes$MGI.Gene.Marker.ID %in% rownames(resPat[o,])[1:ii])
})


library(aggregation)
pvalMatEdgeR <- do.call(cbind,lapply(lrtList,function(x) x$table$PValue))
fisherP <- apply(pvalMatEdgeR,1,fisher)
names(fisherP) <- rownames(lrtList[[1]]$table)
o <- order(fisherP)
sumFoundFisher <- sapply(1:nrow(countsFiltered), function(ii){
  sum(cellGenes$MGI.Gene.Marker.ID %in% names(fisherP[o])[1:ii])
})

library(scales)
plot(x=1:nrow(countsFiltered), y=sumFoundPat, type="l", log="y", col="darkseagreen3", lwd=2)
#lines(x=1:nrow(countsFiltered), y=sumFoundEnd, col="red")
lines(x=1:nrow(countsFiltered), y=sumFoundFisher, col="red", lwd=2)
for(i in 1:15) lines(x=1:nrow(countsFiltered), y=sumFoundMat[,i], col=alpha(colorRampPalette(c("lightblue","darkblue"))(15)[i],.5))


# we can also do a pattern test on a fixed region...
# first check which cluster comparisons are performing that good.
# 1, 6, 11 are performing best. Note that all of these involve a comparison of the INP1 cluster with the three other SUS clusters. Hence, it is this cluster that matters; the differences are local and this is thus probably not a good example.
```










































# Junk

```{R}
#Look for neuronal gene sets?
olfactoryBehavior=read.table("~/Downloads/GO_term_summary_20181002_051017.txt", header=TRUE, sep="\t", row.names=NULL)
olfactoryBehavior$MGI.Gene.Marker.ID <- as.character(olfactoryBehavior$MGI.Gene.Marker.ID)
mean(olfactoryBehavior$MGI.Gene.Marker.ID%in%rownames(counts))


# microvilous
microvillusGenes=read.table("~/Downloads/GO_term_summary_20181002_080501.txt", header=TRUE, sep="\t", row.names=NULL)
microvillusGenes$MGI.Gene.Marker.ID <- as.character(microvillusGenes$MGI.Gene.Marker.ID)
mean(microvillusGenes$MGI.Gene.Marker.ID%in%rownames(counts))

# edgeR
#11 and 15 are unique lineage 2 clusters.
#7 and 4 are unique lineage 3 clusters.
L <- matrix(0,nrow=ncol(fit$coefficients),ncol=4)
rownames(L) <- colnames(fit$coefficients)
L[c("clus.labels11","clus.labels7"),1] <- c(1,-1)
L[c("clus.labels11","clus.labels4"),2] <- c(1,-1)
L[c("clus.labels15","clus.labels7"),3] <- c(1,-1)
L[c("clus.labels15","clus.labels4"),4] <- c(1,-1)
lrtListMicro <- list()
for(i in 1:ncol(L)) lrtListMicro[[i]] <- zinbwave::glmWeightedF(fit,contrast=L[,i])

sumFoundMatMicro <- matrix(NA, nrow=nrow(countsFiltered), ncol=15)
for(i in 1:length(lrtList)){
  pvals <- lrtList[[i]]$table$PValue
  o <- order(pvals)
  pvalsO <- pvals[o]
  namesO <- rownames(lrtList[[i]])[o]
  sumFound <- sapply(1:nrow(countsFiltered), function(ii){
    sum(microvillusGenes$MGI.Gene.Marker.ID %in% namesO[1:ii])
  })
  sumFoundMatMicro[,i] <- sumFound
}

o <- order(resPat$waldStat_2vs3, decreasing=TRUE)
sumFoundPatMicro <- sapply(1:nrow(countsFiltered), function(ii){
  sum(cellGenes$MGI.Gene.Marker.ID %in% rownames(resPat[o,])[1:ii])
})


plot(x=1:nrow(countsFiltered), y=sumFoundPatMicro, type="l", log="y", col="darkseagreen3", lwd=2)
for(i in 1:15) lines(x=1:nrow(countsFiltered), y=sumFoundMatMicro[,i], col=alpha(colorRampPalette(c("lightblue","darkblue"))(4)[i],.5))


```













# make heatmap of unique GAM genes compared to genes shared between GAM and edgeR

```{r}
countsGAMOnly <- countsFiltered[uniqueGamEndId,]
# lineage 1
countsGAMOnly1 <- countsGAMOnly[,gamList[[1]]$model$l1==1]
o1 <- order(gamList[[1]]$model$t1[gamList[[1]]$model$l1==1])
countsGAMOnly1 <- countsGAMOnly1[,o1]
countsGAMOnly1 <- countsGAMOnly[!rowSums(countsGAMOnly1)==0,]
countsGAMOnly1 <- t(scale(t(countsGAMOnly1)))
## cluster gene expression
library(clusterExperiment)
geneClConsensusMerge <- function(geneClusters, cutoff, lineageName) {
  ceg <- combineMany(geneClusters, proportion=0.7, minSize=10,propUnassigned=0.9, clusterFunction="hierarchical01")
  plotClusters(ceg)
  plotCoClustering(ceg)
  ceg <- makeDendrogram(ceg,reduceMethod="none", filterIgnoresUnassigned=TRUE)
  ceg <- mergeClusters(ceg, mergeMethod="locfdr", plotType="colorblock", cutoff=cutoff, DEMethod="limma")
  assign(paste0("ceg", lineageName), ceg)
  #save(list=paste0("ceg", lineageName), file=file.path(gClust_dir,paste0(expt_str, "_", lineageName, "_geneCl_final.Rda")))
  return(ceg)
}
countsCluster1 <- t(countsGAMOnly1)
clGeneLin1 <- clusterMany(countsCluster1, reduceMethod="none", clusterFunction=c("pam"), minSizes=1, subsample=F, sequential=F, ncores=detectCores()-1, random.seed=88, run=TRUE, ks = 4:15, verbose=TRUE, isCount=FALSE)
consLin1 <- geneClConsensusMerge(clGeneLin1, 0.04, "neuronal")






```


# Session information

```{r}
sessionInfo()
```
