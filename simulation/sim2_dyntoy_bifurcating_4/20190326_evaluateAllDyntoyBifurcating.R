library(slingshot)
library(RColorBrewer)
library(mgcv)
library(tradeR)
library(edgeR)
library(rafalib)
library(wesanderson)
palette(wes_palette("Darjeeling1", 10, type="continuous"))
datasetClusters <- read.table("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/datasetClustersSlingshot.txt", header=TRUE)

  FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
  }

#pdf("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/lineages.pdf")
for(datasetIter in 1:10){

  pdf(paste0("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/dataset",datasetIter,".pdf"))


  data <- readRDS(paste0("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/datasets/20190326_dyntoyDataset_", datasetIter, ".rds"))
  counts <- t(data$counts)
  falseGenes <- data$tde_overall$feature_id[data$tde_overall$differentially_expressed]
  nullGenes <- data$tde_overall$feature_id[!data$tde_overall$differentially_expressed]

  # get milestones
  gid <- data$prior_information$groups_id
  gid <- gid[match(colnames(counts),gid$cell_id),]

  pal <- wes_palette("Zissou1", 12, type = "continuous")
  truePseudotime <- data$prior_information$timecourse_continuous
  g <- Hmisc::cut2(truePseudotime,g=12)

  # quantile normalization
  normCounts <- FQnorm(counts)

  ## dim red
  pca <- prcomp(log1p(t(normCounts)), scale. = FALSE)
  rd <- pca$x[,1:3]
  ## cluster
  nClusters <- datasetClusters$nClusters[datasetIter]
  set.seed(5)
  cl <- kmeans(rd, centers = nClusters)$cluster
   rafalib::mypar(mfrow=c(1,2))
   plot(rd, col = wes_palette("Darjeeling1", 10, type="continuous")[cl], pch=16, asp = 1)
   legend("topleft",legend=as.character(1:nClusters),col=wes_palette("Darjeeling1", 10, type="continuous"),pch=16,cex=2/3,bty='n')
   plot(rd, col = brewer.pal(8,"Dark2")[as.numeric(as.factor(gid$group_id))+1], pch=16, asp = 1)
  #
  # rafalib::mypar(mfrow=c(1,3))
  # plot(rd, col = pal[cl], pch=16, asp = 1)
  # legend("topleft",legend=as.character(1:nClusters),col=pal,pch=16,cex=2/3,bty='n')
  # plot(rd, col = pal[g], pch=16, asp = 1)
  # plot(rd, col = as.numeric(as.factor(gid$group_id))+1, pch=16, asp = 1)

  #lineages
  lin <- getLineages(rd, cl, start.clus=datasetClusters$start[datasetIter], end.clus=c(datasetClusters$end1[datasetIter], datasetClusters$end2[datasetIter]))
  plot(rd, col = pal[g], pch=16, asp = 1)
  lines(lin,lwd=2)
  #curves
  crv <- getCurves(lin)
  plot(rd, col = pal[g], pch=16, asp = 1)
  lines(crv, lwd=2, col="black")
# } #for lineages.pdf
# dev.off()

  ### fit smoothers on raw data
  cWeights <- slingCurveWeights(crv)
  pseudoT <- slingPseudotime(crv, na=FALSE)
  gamList <- fitGAM(counts, pseudotime=pseudoT, cellWeights=cWeights, verbose=FALSE)

  # Test for differences at end point
  endRes <- diffEndTest(gamList)
  #hist(endRes$pvalue)
  deGenesEndGam <- rownames(counts)[which(p.adjust(endRes$pvalue,"fdr")<=0.05)]
  mean(falseGenes%in%deGenesEndGam) #TPR
  mean(!deGenesEndGam%in%falseGenes) #FDR

  ## tradeR pattern test
  resPattern <- patternTest(gamList)
  #hist(resPattern$pval)
  deGenesPattern <- rownames(counts)[which(p.adjust(resPattern$pval,"fdr")<=0.05)]
  length(deGenesPattern)
  mean(falseGenes%in%deGenesPattern) #TPR
  mean(!deGenesPattern%in%falseGenes) #FDR

  # edgeR analysis on final clusters
  clF <- as.factor(cl)
  design <- model.matrix(~clF)

  library(edgeR)
  d <- DGEList(counts)
  d <- calcNormFactors(d)
  d <- estimateDisp(d, design)
  fit <- glmFit(d, design)
  L <- matrix(0, nrow=ncol(fit$coefficients), ncol=1)
  rownames(L) <- colnames(fit$coefficients)
  endClusters <- unlist(c(datasetClusters[datasetIter,c("end1","end2")]))
  if(1 %in% endClusters){
    not1Cluster <- endClusters[!endClusters==1]
    L[not1Cluster,1] <- 1
  } else {
    L[endClusters,1] <- c(1,-1)
  }
  L
  lrt <- glmLRT(fit, contrast=L) #cluster 2 vs cluster 1 (end clusters)
  deGenesEdgeR <- rownames(lrt$table)[p.adjust(lrt$table$PValue,"fdr")<=0.05]
  mean(falseGenes%in%deGenesEdgeR) #TPR
  mean(deGenesEdgeR%in%nullGenes) #FDR


  # Monocle BEAM analysis
  ### old monocle BEAM analysis
  library(monocle,lib.loc="~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/methodsPaper/")
  featureInfo <- data.frame(gene_short_name=rownames(counts))
  rownames(featureInfo) <- rownames(counts)
  fd <- new("AnnotatedDataFrame", featureInfo)
  cds <- newCellDataSet(cellData=normCounts, featureData=fd, expressionFamily=negbinomial.size())
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  cds <- reduceDimension(cds)
  cds <- orderCells(cds)
  print(plot_cell_trajectory(cds, color_by = "State"))
  if(nlevels(phenoData(cds)$State)==1){
    BEAM_res <- data.frame(pval=rep(1,nrow(counts)))
  } else {
    BEAM_res <- BEAM(cds,  cores = 1)
  }

  if(nlevels(phenoData(cds)$State)==1){
    resEndMon <- data.frame(pvalue=rep(1,nrow(counts)))
    resPatternMon <- data.frame(pvalue=rep(1,nrow(counts)))
  } else {
    # tradeR downstream of Monocle 2
    phenoData(cds)$milestone <- gid$group_id
    print(plot_cell_trajectory(cds, color_by = "milestone"))
    plot(rd, col = as.numeric(as.factor(gid$group_id))+1, pch=16, asp = 1) ; legend("topleft", paste0("M",1:4), col=2:5,pch=16)
    #M2 and M4 are the branches we need to compare.
    # these correspond to State 1 and State 6.

    ptMon <- matrix(phenoData(cds)$Pseudotime, nrow=ncol(counts), ncol=2, byrow=FALSE)
    stateMon <- phenoData(cds)$State
    #plot(x=ptMon,y=stateMon)
    # use milestones to derive true weights.
    # use slingshot weights to identify curves.
    cellWeightsMon <- matrix(0, nrow=ncol(cds), ncol=2)
    rownames(cellWeightsMon) <- colnames(counts)
    cellWeightsMon[gid$group_id %in% c("M1","M3"),] <- c(1/2,1/2)
    cellWeightsMon[gid$group_id %in% "M2",which.max(colSums(cWeights[gid$group_id %in% "M2",]))] <- 1
    cellWeightsMon[gid$group_id %in% "M4",which.max(colSums(cWeights[gid$group_id %in% "M4",]))] <- 1

    gamListMon <- fitGAM(counts, pseudotime=ptMon, cellWeights=cellWeightsMon, verbose=FALSE)
    resPatternMon <- patternTest(gamListMon)
    resEndMon <- diffEndTest(gamListMon)
  }

  # tradeR on true time and weights
  ### tradeR on true pseudotime
  pst <- matrix(truePseudotime, nrow=ncol(counts),ncol=2, byrow=FALSE)
  gamListTrueTime <- fitGAM(counts, pseudotime=pst, cellWeights=slingCurveWeights(crv), verbose=FALSE)
  # end point test
  waldEndPointResTrueTime <- diffEndTest(gamListTrueTime)
  padjWaldTrueTime <- p.adjust(waldEndPointResTrueTime$pvalue,"fdr")
  sum(padjWaldTrueTime<=0.05)
  # pattern test
  patternResTrueTime <- patternTest(gamListTrueTime)
  padjPatternTrueTime <- p.adjust(patternResTrueTime$pvalue,"fdr")
  sum(padjPatternTrueTime<=0.05, na.rm=TRUE)

  # GPfates
  logCpm <- edgeR::cpm(normCounts, prior.count=.125, log=TRUE)
  sampleInfo <- data.frame(global_pseudotime=truePseudotime)
  rownames(sampleInfo) <- colnames(counts)
  write.table(logCpm, file="~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/simDyntoyLogCpm.txt", row.names=TRUE, col.names=TRUE, quote=FALSE)
  write.table(sampleInfo, file="~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/sampleInfoSimDyntoy.txt", row.names=TRUE, col.names=TRUE, quote=FALSE)
  write.table(datasetIter, file="~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/currIter.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
  #run with 20190206_runGPfates_simDyntoy.py
  system("python3 /Users/koenvandenberge/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/20190207_runGPfates_simDyntoy_bifurcating4.py")
  # import GPfates output
  GPfatesWeights <- read.table("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/GPfatesWeights.txt", header=FALSE)
  GPfatesBif <- read.table("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/GPfatesBifStats.txt", header=FALSE)
  colnames(GPfatesBif) <- c("bif_ll", "amb_ll", "shuff_bif_ll", "shuff_amb_ll", "phi0_corr", "D", "shuff_D")


  # GPfates weights + tradeR
  # based on true pseudotime
  gamListGPfatesTrueTime <- fitGAM(counts, pseudotime=pst, cellWeights=GPfatesWeights, verbose=FALSE)
  # end point test
  waldEndPointResTrueTimeGPfates <- diffEndTest(gamListGPfatesTrueTime)
  # pattern test
  patternResTrueTimeGPfates <- patternTest(gamListGPfatesTrueTime)

  ## performance
  library(iCOBRA)
  library(scales)
  truth <- as.data.frame(matrix(rep(0,nrow(counts)), dimnames=list(rownames(counts),"status")))
  truth[falseGenes,"status"] <- 1
  saveRDS(truth, file=paste0("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/datasets/truth",datasetIter,".rds"))


  pval <- data.frame( tradeR_slingshot_end=endRes$pval,
                      tradeR_slingshot_pattern=resPattern$pval,
                      BEAM=BEAM_res$pval,
                      edgeR=lrt$table$PValue,
                      tradeR_GPfates_end=waldEndPointResTrueTimeGPfates$pval,
                      tradeR_GPfates_pattern=patternResTrueTimeGPfates$pval,
                      tradeR_Monocle2_end=resEndMon$pvalue,
                      tradeR_Monocle2_pattern=resPatternMon$pvalue,
                        row.names=rownames(counts))
  score <- data.frame(GPfates=GPfatesBif$D,
                        row.names=rownames(counts))
  cobra <- COBRAData(pval=pval, truth=truth, score=score)
  saveRDS(cobra, file=paste0("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/datasets/cobra",datasetIter,".rds"))

  dev.off()
}
