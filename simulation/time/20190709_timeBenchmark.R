dataset <- readRDS(paste0("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/time/bigDyntoyDataset.rds"))
library(RColorBrewer)
library(wesanderson)
library(slingshot)
library(tradeSeq)
library(ImpulseDE2)
library(profmem)
  FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
  }

counts <- as.matrix(t(dataset$counts))
normCounts <- FQnorm(counts)

pal <- wes_palette("Zissou1", 12, type = "continuous")
truePseudotime <- dataset$prior_information$timecourse_continuous
g <- Hmisc::cut2(truePseudotime,g=12)

### prepare slingshot
pca <- prcomp(log1p(t(normCounts)), scale. = FALSE)
rd <- pca$x[,1:3]
## cluster
nClusters <- 8
set.seed(5)
cl <- kmeans(rd, centers = nClusters)$cluster
rafalib::mypar(mfrow=c(1,2))

rafalib::mypar(mfrow=c(1,2))
plot(rd, col = cl, pch=16, asp = 1)
legend("topleft",legend=as.character(1:nClusters),col=1:nClusters,pch=16,cex=2/3,bty='n')
plot(rd, col = pal[g], pch=16, asp = 1)
lin <- getLineages(rd, cl, start.clus=5, end.clus=c(7,8))
plot(rd, col = pal[g], pch=16, asp = 1)
lines(lin,lwd=2, col="black")
crv <- getCurves(lin)
plot(rd, col = pal[g], pch=16, asp = 1)
lines(crv, lwd=2, col="black")

pt <- slingPseudotime(crv, na=FALSE)
cWeights <- slingCurveWeights(crv)

### prepare Monocle
library(monocle,lib.loc="~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/methodsPaper/")
featureInfo <- data.frame(gene_short_name=rownames(counts))
rownames(featureInfo) <- rownames(counts)
fd <- new("AnnotatedDataFrame", featureInfo)
cds <- newCellDataSet(cellData=normCounts, featureData=fd, expressionFamily=negbinomial.size())


nCells <- c(250, 500, 1000, 5000, 10000)
fitGamTime <- list()
fitGamMem <- list()
diffEndTradeSeqTime <- list()
diffEndTradeSeqMem <- list()
patternTradeSeqTime <- list()
patternTradeSeqMem <- list()
BEAMTime <- list()
BEAMMem <- list()
ImpulseTime <- list()
ImpulseMem <- list()
for(ii in 1:length(nCells)){

  set.seed(ii)

  n <- nCells[ii]
  id <- sample(1:nrow(counts), n)

  # tradeSeq fitting
  fitGamTime[[ii]] <- system.time({
    #fitGamMem[[ii]] <- total(profmem({
        gamList <- fitGAM(normCounts[,id], pseudotime=pt[id,], cellWeights=cWeights[id,], parallel=FALSE)
      #}))/1e6
    })

  # tradeSeq diffEndTest (similar for all other tests)
  diffEndTradeSeqTime[[ii]] <- system.time({
      diffEndTradeSeqMem[[ii]] <- total(profmem({
        det <- diffEndTest(gamList)
      }))/1e6
    })

  # tradeSeq patternTest
  patternTradeSeqTime[[ii]] <- system.time({
    patternTradeSeqMem[[ii]] <- total(profmem({
      pat <- patternTest(gamList)
    }))/1e6
  })

  # BEAM
  cdsSubset <- newCellDataSet(cellData=normCounts[,id], featureData=fd, expressionFamily=negbinomial.size())
  cdsSubset <- estimateSizeFactors(cdsSubset)
  cdsSubset <- estimateDispersions(cdsSubset)
  cdsSubset <- reduceDimension(cdsSubset)
  cdsSubset <- orderCells(cdsSubset)

  BEAMTime[[ii]] <- system.time({
    BEAMMem[[ii]] <- total(profmem({
      BEAM_res <- BEAM(cdsSubset,  cores = 1)
    }))/1e6
  })


  # ImpulseDE2
  if(n <= 1000){ #only up to 1000 cells, takes extremely long.
    m <- gamList[[1]]
    branch <- rep(NA, n)
    time <- rep(NA, n)
    branch[m$model$l1==1] <- "A"
    ptSub <- pt[id,]
    time[m$model$l1==1] <- ptSub[m$model$l1==1,1]
    branch[m$model$l2==1] <- "B"
    time[m$model$l2==1] <- ptSub[m$model$l2==1,2]
    condition <- ifelse(branch=="A","control","case")
    # ImpulseDE2 cannot handle datasets where every gene contains at least one zero.
    # since then it errors on the DESeq2 size factor estimation. Since the data is
    # already quantile normalized, we provide a size factor of 1 as input.
    sf <- rep(1,n) #data is already normalized
    names(sf) <- colnames(normCounts[,id])
    # dfAnnotation
    dfAnn <- data.frame(Sample=colnames(normCounts[,id]), Condition=condition, Time=time)
    # Even with supplied size factors, it errors because it still attempts to calculate
    # the size factors for dispersion estimation with DESeq2. Do it manually.
    library(DESeq2)
    dds <- suppressWarnings( DESeqDataSetFromMatrix(
                      countData = round(normCounts[,id]),
                      colData = dfAnn,
                      design = ~Condition+Condition:Time) )
    dds <- estimateSizeFactors(dds, type="poscounts")
    dds <- estimateDispersions(dds)
    vecDispersionsInv <- mcols(dds)$dispersion
    vecDispersions <- 1/vecDispersionsInv
    names(vecDispersions) <- rownames(dds)

    ImpulseTime[[ii]] <- system.time({
      ImpulseMem[[ii]] <- total(profmem({
        imp <- runImpulseDE2(matCountData=round(normCounts), dfAnnotation=dfAnn, boolCaseCtrl=TRUE, vecSizeFactorsExternal=sf, vecDispersionsExternal=vecDispersions, scaNProc=1)
      }))/1e6
    })

  }
}

saveRDS(fitGamTime, file=paste0("/home/compomics/Dropbox/compomicsVM/time/fitGamTime.rds"))
saveRDS(fitGamMem, file=paste0("/home/compomics/Dropbox/compomicsVM/time/fitGamMem.rds"))
saveRDS(diffEndTradeSeqTime, file=paste0("/home/compomics/Dropbox/compomicsVM/time/diffEndTradeSeqTime.rds"))
saveRDS(diffEndTradeSeqMem, file=paste0("/home/compomics/Dropbox/compomicsVM/time/diffEndTradeSeqMem.rds"))
saveRDS(patternTradeSeqTime, file=paste0("/home/compomics/Dropbox/compomicsVM/time/patternTradeSeqTime.rds"))
saveRDS(patternTradeSeqMem, file=paste0("/home/compomics/Dropbox/compomicsVM/time/patternTradeSeqMem.rds"))
saveRDS(BEAMTime, file=paste0("/home/compomics/Dropbox/compomicsVM/time/BEAMTime.rds"))
saveRDS(BEAMMem, file=paste0("/home/compomics/Dropbox/compomicsVM/time/BEAMMem.rds"))
saveRDS(ImpulseTime, file=paste0("/home/compomics/Dropbox/compomicsVM/time/ImpulseTime.rds"))
saveRDS(ImpulseMem, file=paste0("/home/compomics/Dropbox/compomicsVM/time/ImpulseMem.rds"))
