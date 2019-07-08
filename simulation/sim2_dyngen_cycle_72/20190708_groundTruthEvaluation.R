library(slingshot)
library(RColorBrewer)
library(mgcv)
library(tradeSeq)
library(edgeR)
library(rafalib)
library(wesanderson)
 library(BiocParallel)
library(doParallel)
NCORES <- BiocParallel::bpparam()$workers

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
for(datasetIter in 1:10){

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

  ### tradeSeq on true pseudotime
  cWeights <- rep(1,ncol(counts))
  pst <- matrix(truePseudotime, nrow=ncol(counts), ncol=1, byrow=FALSE)
  gamListTrueTime <- fitGAM(counts, pseudotime=pst, cellWeights=cWeights)
  assocTestTrueRes <- associationTest(gamListTrueTime)

  ### tradeSeq with 3 knots (Monocle default)
  gamListTrueTime_3k <- fitGAM(counts, pseudotime=pst, cellWeights=cWeights, nknots=3)
  assocTestTrueRes_3k <- associationTest(gamListTrueTime_3k)

  ### tradeSeq with 5 knots (optimal acc. to AIC)
  gamListTrueTime_5k <- fitGAM(counts, pseudotime=pst, cellWeights=cWeights, nknots=5)
  assocTestTrueRes_5k <- associationTest(gamListTrueTime_5k)

  # Performance plots

  ########
  # FDP-TPR
  ########
  library(iCOBRA)
  library(scales)
  truth <- as.data.frame(matrix(rep(0,nrow(counts)), dimnames=list(rownames(counts),"status")))
  truth[falseGenes,"status"] <- 1

  ### estimated pseudotime
  pval <- data.frame( tradeSeq_10k=assocTestTrueRes$pval,
                      tradeSeq_3k=assocTestTrueRes_3k$pval,
                      tradeSeq_5k=assocTestTrueRes_5k$pval,
                        row.names=rownames(counts))
  cobra <- COBRAData(pval=pval, truth=truth)
  saveRDS(cobra, file=paste0("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyngen_cycle_72/datasets/cobraGroundTruth",datasetIter,".rds"))
}
