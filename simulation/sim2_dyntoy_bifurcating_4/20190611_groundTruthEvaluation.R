library(slingshot)
library(RColorBrewer)
library(mgcv)
library(tradeSeq)
library(edgeR)
library(rafalib)
library(wesanderson)
library(BiocParallel)
library(doParallel)
NCORES <- 2
palette(wes_palette("Darjeeling1", 10, type="continuous"))
datasetClusters <- read.table("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/datasetClustersSlingshot.txt", header=TRUE)
source("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/20190611_helper.R")

  FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
  }

#pdf("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/lineages.pdf")
for(datasetIter in 1:10){

  #pdf(paste0("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/dataset",datasetIter,".pdf"))


  data <- readRDS(paste0("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/datasets/20190326_dyntoyDataset_", datasetIter, ".rds"))
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
  normCounts <- round(FQnorm(counts))

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
   plot(rd, col = brewer.pal(8,"Dark2")[as.numeric(as.factor(gid$group_id))], pch=16, asp = 1)
   legend("topright",paste0("M",1:length(unique(gid$group_id))), col=1:4, pch=16)

  rafalib::mypar(mfrow=c(1,3))
  plot(rd, col = pal[cl], pch=16, asp = 1)
  legend("topleft",legend=as.character(1:nClusters),col=pal,pch=16,cex=2/3,bty='n')
  plot(rd, col = pal[g], pch=16, asp = 1)
  plot(rd, col = as.numeric(as.factor(gid$group_id))+1, pch=16, asp = 1)

  #lineages
  lin <- getLineages(rd, cl, start.clus=datasetClusters$start[datasetIter], end.clus=c(datasetClusters$end1[datasetIter], datasetClusters$end2[datasetIter]))
  plot(rd, col = pal[g], pch=16, asp = 1)
  lines(lin,lwd=2)
  #curves
  crv <- getCurves(lin)
  plot(rd, col = pal[g], pch=16, asp = 1)
  lines(crv, lwd=2, col="black")

  ### tradeSeq: fit smoothers on truth data
  trueWeights <- getWeightsBifurcation(data, crv)
  trueT <- matrix(truePseudotime, nrow=length(truePseudotime), ncol=2, byrow=FALSE)
  gamListTruth <- fitGAM(counts, pseudotime=trueT, cellWeights=trueWeights)

  endRes <- diffEndTest(gamListTruth)
  patternRes <- patternTest(gamListTruth)

  # Monocle BEAM analysis
  ### Monocle 2 BEAM analysis
  library(monocle,lib.loc="~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/methodsPaper/")
  featureInfo <- data.frame(gene_short_name=rownames(counts))
  rownames(featureInfo) <- rownames(counts)
  fd <- new("AnnotatedDataFrame", featureInfo)
  cds <- newCellDataSet(cellData=normCounts, featureData=fd, expressionFamily=negbinomial.size())
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  cds <- reduceDimension(cds)
  cds <- orderCells(cds)
  set.seed(11)
  branch <- rep(NA, ncol(counts))
  branch[trueWeights[,1]==1] <- "A"
  branch[trueWeights[,2]==1] <- "B"
  branch[is.na(branch)] <- c("A","B")[sample(1:2, size=sum(is.na(branch)), replace=TRUE)]
  phenoData(cds)$Branch <- as.factor(branch)
  phenoData(cds)$original_cell_id <- colnames(counts)
  phenoData(cds)$Pseudotime <- truePseudotime
  BEAM_true <- BEAM_kvdb(cds,  cores = 1)
  #Note: original BEAM throws away cells (supposedly the root state), this function does not.

  ## ImpulseDE2
  library(ImpulseDE2)
  condition <- as.character(branch)
  condition <- ifelse(condition=="A","control","case")
  # ImpulseDE2 cannot handle datasets where every gene contains at least one zero.
  # since then it errors on the DESeq2 size factor estimation. Since the data is
  # already quantile normalized, we provide a size factor of 1 as input.
  sf <- rep(1,ncol(normCounts)) #data is already normalized
  names(sf) <- colnames(normCounts)
  # dfAnnotation
  dfAnn <- data.frame(Sample=colnames(counts), Condition=condition, Time=truePseudotime)
  # Even with supplies size factors, it errors because it still attempts to calculate
  # the size factors for dispersion estimation with DESeq2. Do it manually.
  library(DESeq2)
  dds <- suppressWarnings( DESeqDataSetFromMatrix(
                    countData = round(normCounts),
                    colData = dfAnn,
                    design = ~Condition+Condition:Time) )
  dds <- estimateSizeFactors(dds, type="poscounts")
  dds <- estimateDispersions(dds)
  vecDispersionsInv <- mcols(dds)$dispersion
  vecDispersions <- 1/vecDispersionsInv
  names(vecDispersions) <- rownames(dds)
  imp <- runImpulseDE2(matCountData=round(normCounts), dfAnnotation=dfAnn, boolCaseCtrl=TRUE, vecSizeFactorsExternal=sf, vecDispersionsExternal=vecDispersions, scaNProc=2)

  # # GPfates
  # ## for GPfates, we would have to adapt the OMGP fitting. However, it is not possible to provide true weights to the OMGP model.
  # logCpm <- edgeR::cpm(normCounts, prior.count=.125, log=TRUE)
  # sampleInfo <- data.frame(global_pseudotime=truePseudotime)
  # rownames(sampleInfo) <- colnames(counts)
  # write.table(logCpm, file="~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/simDyntoyLogCpm.txt", row.names=TRUE, col.names=TRUE, quote=FALSE)
  # write.table(sampleInfo, file="~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/sampleInfoSimDyntoy.txt", row.names=TRUE, col.names=TRUE, quote=FALSE)
  # write.table(datasetIter, file="~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/currIter.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
  # #run with 20190206_runGPfates_simDyntoy.py
  # system("python3 /Users/koenvandenberge/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/20190207_runGPfates_simDyntoy_bifurcating4.py")
  # # import GPfates output
  # GPfatesWeights <- read.table("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/GPfatesWeights.txt", header=FALSE)
  # GPfatesBif <- read.table("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/GPfatesBifStats.txt", header=FALSE)
  # colnames(GPfatesBif) <- c("bif_ll", "amb_ll", "shuff_bif_ll", "shuff_amb_ll", "phi0_corr", "D", "shuff_D")


  ## performance
  library(iCOBRA)
  library(scales)
  truth <- as.data.frame(matrix(rep(0,nrow(counts)), dimnames=list(rownames(counts),"status")))
  truth[falseGenes,"status"] <- 1


  pval <- data.frame( tradeSeq_diffEnd=endRes$pval,
                      tradeSeq_pattern=patternRes$pval,
                      BEAM=BEAM_true$pval,
                      ImpulseDE2=imp$dfImpulseDE2Results$p,
                        row.names=rownames(counts))
  cobra <- COBRAData(pval=pval, truth=truth)
  saveRDS(cobra, file=paste0("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/datasets/groundTruthCobra",datasetIter,".rds"))

  #dev.off()
}
