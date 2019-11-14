library(slingshot)
library(RColorBrewer)
library(mgcv)
library(tradeSeq)
library(edgeR)
library(rafalib)
library(wesanderson)
library(BiocParallel)
library(doParallel)
library(iCOBRA)
library(scales)
NCORES <- 2

  FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
  }

dataAll <- readRDS("~/Dropbox/research/PhD/research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyngen_cycle_72/datasets/datasets_for_koen.rds")

RNGversion("3.5.0")
for(datasetIter in 1:10){

  pdf(paste0("~/Dropbox/research/PhD/research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyngen_cycle_72/dataset",datasetIter,"_v2IncludingEdgeR.pdf"))

  data <- dataAll[[datasetIter]]
  counts <- as.matrix(t(data$counts))
  falseGenes <- data$feature_info$gene_id[!data$feature_info$housekeeping]
  nullGenes <- data$feature_info$gene_id[data$feature_info$housekeeping]

  pal <- wes_palette("Zissou1", 12, type = "continuous")
  truePseudotime <- data$prior_information$timecourse_continuous
  g <- Hmisc::cut2(truePseudotime,g=12)

  # quantile normalization
  normCounts <- FQnorm(counts)

  ## dim red
  pca <- prcomp(log1p(t(normCounts)), scale. = FALSE)
  rd <- pca$x[,1:2]
  plot(rd, pch=16, asp = 1, col=pal[g])

  library(princurve)
  pcc <- principal_curve(rd, smoother="periodic_lowess")
  lines(x=pcc$s[order(pcc$lambda),1], y=pcc$s[order(pcc$lambda),2], col="red", lwd=2)

  # fit smoothers on raw data
  cWeights <- rep(1,ncol(counts))
  pseudoT <- matrix(pcc$lambda,nrow=ncol(counts),ncol=1)
  gamList <- fitGAM(counts, pseudotime=pseudoT, cellWeights=cWeights, nknots=5)

  # Test for association of expression with the trajectory
  assocTestRes <- associationTest(gamList)
  #hist(assocTestRes$pvalue)
  
  # edgeR in PCA space
  X <- predict(gamList[[1]], type="lpmatrix")
  d <- DGEList(counts)
  d <- calcNormFactors(d)
  design <- model.matrix(~ -1 + X)
  d <- estimateDisp(d, design)
  fit <- glmFit(d, design)
  lrt <- glmLRT(fit, coef=2:ncol(fit$coefficients))

  # Monocle 3
  library(monocle)
  fd <- data.frame(gene_short_name=rownames(counts))
  fd <- new("AnnotatedDataFrame",fd)
  rownames(fd) <- rownames(counts)
  pd <- data.frame(cellid=colnames(counts))
  pd <- new("AnnotatedDataFrame",pd)
  rownames(pd) <- colnames(counts)
  cds <- newCellDataSet(counts, featureData=fd, phenoData=pd)
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  cds <- preprocessCDS(cds, num_dim = 20)
  cds <- reduceDimension(cds, reduction_method = 'UMAP')
  cds <- partitionCells(cds)
  cds <- learnGraph(cds,  RGE_method = 'SimplePPT')
  pr_graph_test <- try(principalGraphTest(cds, k=3, cores=1))
  if(class(pr_graph_test) == "try-error") {
    pr_graph_test <- data.frame(pval = rep(1, nrow(counts)))
  }
  print(plot_cell_trajectory(cds,color_by="louvain_component") + ggtitle(paste0("Dataset ",datasetIter))) # fails to discover cyclic pattern.

  # slingshot in UMAP space
  umapRD <- t(cds@reducedDimS)
  plot(umapRD, pch=16, asp = 1, col=pal[g], main="UMAP")
  pccUMAP <- principal_curve(umapRD, smoother="periodic_lowess")
  lines(x=pccUMAP$s[order(pccUMAP$lambda),1], y=pccUMAP$s[order(pccUMAP$lambda),2], col="red", lwd=2)

  # tradeSeq on UMAP slingshot
  pseudoT <- matrix(pccUMAP$lambda,nrow=ncol(counts),ncol=1)
  gamListUMAP <- fitGAM(counts, pseudotime=pseudoT, cellWeights=cWeights, nknots=5)
  assocTestResUMAP <- associationTest(gamListUMAP)



  ########
  # FDP-TPR
  ########
  truth <- as.data.frame(matrix(rep(0,nrow(counts)), dimnames=list(rownames(counts),"status")))
  truth[falseGenes,"status"] <- 1

  ### estimated pseudotime
  pval <- data.frame( tradeSeq_slingshot_assoc=assocTestRes$pval,
                      tradeSeq_slingshot_UMAP_assoc=assocTestResUMAP$pval,
                      Monocle3=pr_graph_test$pval,
                      edgeR=lrt$table$PValue,
                        row.names=rownames(counts))
  cobra <- COBRAData(pval=pval, truth=truth)
  saveRDS(cobra, file=paste0("~/Dropbox/research/PhD/research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyngen_cycle_72/datasets/cobra",datasetIter,"_v2IncludingEdgeR.rds"))
  dev.off()
}
