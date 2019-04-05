library(slingshot)
library(RColorBrewer)
library(mgcv)
library(tradeR)
library(edgeR)
library(rafalib)
library(wesanderson)

  FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
  }

dataAll <- readRDS("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyngen_cycle_72/datasets/datasets_for_koen.rds")


#for(datasetIter in seq(1:5,7:10)){
for(datasetIter in 1){

  pdf(paste0("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyngen_cycle_72/dataset",datasetIter,".pdf"))

  data <- dataAll[[datasetIter]]
  counts <- t(data$counts)
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
  gamList <- fitGAM(counts, pseudotime=pseudoT, cellWeights=cWeights, verbose=TRUE)

  # Test for association of expression with the trajectory
  assocTestRes <- associationTest(gamList)
  #hist(assocTestRes$pvalue)

  ### tradeR on true pseudotime
  pst <- matrix(truePseudotime, nrow=ncol(counts), ncol=1, byrow=FALSE)
  gamListTrueTime <- fitGAM(counts, pseudotime=pst, cellWeights=cWeights, verbose=TRUE)
  assocTestTrueRes <- associationTest(gamListTrueTime)
  #hist(assocTestTrueRes$pvalue)

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
  pr_graph_test <- principalGraphTest(cds, k=3, cores=1)
  print(plot_cell_trajectory(cds,color_by="louvain_component")) # fails to discover cyclic pattern.


  # Performance plots

  ########
  # FDP-TPR
  ########
  library(iCOBRA)
  library(scales)
  truth <- as.data.frame(matrix(rep(0,nrow(counts)), dimnames=list(rownames(counts),"status")))
  truth[falseGenes,"status"] <- 1

  ### estimated pseudotime
  pval <- data.frame( tradeR_slingshot_assoc=assocTestRes$pval,
                      Monocle3=pr_graph_test$pval,
                      tradeR_slingshot_assoc_trueTime=assocTestTrueRes$pvalue,
                        row.names=rownames(counts))
  cobra <- COBRAData(pval=pval, truth=truth)
  saveRDS(cobra, file=paste0("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyngen_cycle_72/datasets/cobra",datasetIter,".rds"))
  dev.off()
}
