library(monocle) # Load Monocle
RNGversion("3.5.0")
library(rafalib)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(clusterExperiment)
library(cluster)
library(scales)
gcolpal <- c(brewer.pal(8, "Dark2")[-c(2, 3, 5)],
             brewer.pal(12, "Paired")[c(1, 2, 8, 10, 9)],
             brewer.pal(12, "Set3")[c(7, 8, 12)],
             brewer.pal(8, "Pastel2")[8], brewer.pal(11, "BrBG")[11],
             brewer.pal(11, "PiYG")[1], "cyan", "darkblue", "darkorchid2",
             "brown1", "springgreen1", "deepskyblue4", "darkolivegreen",
             "antiquewhite2")
cell_type_color <- c("Basophils" = "#E088B8",
                     "Dendritic cells" = "#46C7EF",
                     "Eosinophls" = "#EFAD1E",
                     "Erythrocyte" = "#8CB3DF",
                     "Monocytes" = "#53C0AD",
                     "Multipotent progenitors" = "#4EB859",
                     "GMP" = "#D097C4",
                     "Megakaryocytes" = "#ACC436",
                     "Neutrophils" = "#F5918A",
                     'NA' = '#000080')


download.file(
  "https://github.com/statOmics/tradeSeqPaper/raw/master/data/se_paul.rda",
  destfile = "./se_paul.rda")
load("./se_paul.rda")
rd <- reducedDim(se)
set.seed(97)
cl <- kmeans(rd, centers = 7)$cluster
plot(rd, col = brewer.pal(9, "Set1")[cl], pch = 16, asp = 1)
lin <- getLineages(rd, clusterLabels = cl, start.clus = 4)
plot(rd, col = brewer.pal(9, "Set1")[cl], xlab = "UMAP1", ylab = "UMAP2")
lines(lin, lwd = 2)
crv <- getCurves(lin)
plot(rd, col = brewer.pal(9, "Set1")[cl], main = "color by cluster", xlab = "UMAP1", ylab = "UMAP2")
lines(crv, lwd = 2)
plot(rd, col = cell_type_color[colData(se)$cell_type2], main = "color by cell type", xlab = "UMAP1", ylab = "UMAP2", pch = 16)
lines(crv, lwd = 2)
counts <- as.matrix(assays(se)$counts)


gamListPaul <- fitGAM(counts, pseudotime=slingPseudotime(crv,na=FALSE), cellWeights=slingCurveWeights(crv), verbose=TRUE, nknots=6)

# pattern test
patternResPaul <- patternTest(gamListPaul)
sum(p.adjust(patternResPaul$pvalue, "fdr") <= 0.05, na.rm = TRUE)
patternGenes <- rownames(counts)[which(p.adjust(patternResPaul$pvalue, "fdr") <= 0.05)]


#### fitted values
#RSEC on PCA
resRSEC <- clusterExpressionPatterns(gamListPaul, nPoints=100, genes=patternGenes)
clRSEC <- primaryCluster(resRSEC$rsec)

#RSEC directly
resRSECNoDR <- clusterExpressionPatterns(gamListPaul, nPoints=100, genes=patternGenes,
                                         reduceMethod="none")

## same number of clusters with PAM.
set.seed(882)
yhatScaled <- resRSEC$yhatScaled
resPam <- cluster::pam(x=yhatScaled, k=length(unique(clRSEC))-1)
clPam <- resPam$clustering

## the silhouette value is biased to k-means / PAM since it's based on the same distance.
## a non-parametric approach would be to use bootstrapping.
pt <- slingPseudotime(crv,na=FALSE)
cw <- slingCurveWeights(crv)

ariRSEC <- c()
ariRSECNoDR <- c()
ariPam <- c()
for(ii in 1:3){
  set.seed(ii)
  ## bootstrap cells
  bootCells <- sample(1:ncol(counts), replace=TRUE)
  ptBoot <- pt[bootCells,]
  cwBoot <- cw[bootCells,]
  countsBoot <- counts[,bootCells]
  
  ## refit tradeSeq
  glBoot <- fitGAM(countsBoot[patternGenes,], pseudotime=ptBoot, cellWeights=cwBoot, 
                   verbose=TRUE, nknots=6)
  
  ## cluster using RSEC
  rsecBoot <- clusterExpressionPatterns(glBoot, nPoints=100, genes=patternGenes)
  
  ## cluster using RSEC without PCA
  rsecBootNoDR <- clusterExpressionPatterns(glBoot, nPoints=100, genes=patternGenes,
                                            reduceMethod="none")
  
  ## cluster using PAM
  yhatBoot <- rsecBoot$yhatScaled
  pamBoot <- cluster::pam(yhatBoot, k=length(unique(clRSEC))-1)
  
  ## calculate ARI with full data clustering of respective clustering method.
  ariRSEC[ii] <- mclust::adjustedRandIndex(primaryCluster(rsecBoot$rsec),
                            primaryCluster(resRSEC$rsec))
  
  ariRSECNoDR[ii] <- mclust::adjustedRandIndex(primaryCluster(rsecBootNoDR$rsec),
                            primaryCluster(resRSECNoDR$rsec))
  
  ariPam[ii] <- mclust::adjustedRandIndex(resPam$clustering,
                            pamBoot$clustering)
}

ariRSEC
ariRSECNoDR
ariPam




