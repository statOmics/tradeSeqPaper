library(monocle)
library(DESeq2)
library(RColorBrewer)
library(mgcv)
library(slingshot)
library(tradeSeq)
library(microbenchmark)
library(edgeR)
library(here)
library(rafalib)
library(wesanderson)
library(BiocParallel)
library(doParallel)
library(profvis)

## Pre-process ----
NCORES <- 2
palette(wes_palette("Darjeeling1", 10, type = "continuous"))
source(here("simulation", "time", "20190611_helper.R"))

FQnorm <- function(counts) {
  rk <- apply(counts, 2, rank, ties.method = "min")
  counts.sort <- apply(counts, 2, sort)
  refdist <- apply(counts.sort, 1, median)
  norm <- apply(rk, 2, function(r) {
    refdist[r]
  })
  rownames(norm) <- rownames(counts)
  return(norm)
}

for (size in c("small", "big")) {
  dataset <- readRDS(paste0(here("simulation", "time",
                                 paste0(size, "DyntoyDataset.rds"))))
  counts <- t(dataset$counts)
  
  # get milestones
  gid <- dataset$prior_information$groups_id
  gid <- gid[match(colnames(counts), gid$cell_id), ]
  
  pal <- wes_palette("Zissou1", 12, type = "continuous")
  truePseudotime <- dataset$prior_information$timecourse_continuous
  g <- Hmisc::cut2(truePseudotime, g = 12)
  
  # quantile normalization
  normCounts <- round(FQnorm(counts))
  
  
  ## Running slingshot ----
  ## dim red
  pca <- prcomp(log1p(t(normCounts)), scale. = FALSE)
  rd <- pca$x[, 1:3]
  ## cluster
  cl <- as.factor(gid$group_id)
  rafalib::mypar(mfrow = c(1, 2))
  plot(rd, col = brewer.pal(8, "Dark2")[as.numeric(as.factor(gid$group_id))],
       pch = 16, asp = 1)
  legend("topright", paste0("M", 1:length(unique(gid$group_id))), col = 1:4, pch = 16)
  plot(rd, col = pal[g], pch = 16, asp = 1)
  
  # lineages
  lin <- getLineages(rd, as.numeric(cl), start.clus = 1,
                     end.clus = c(2, 4))
  plot(rd, col = pal[g], pch = 16, asp = 1)
  lines(lin, lwd = 2)
  # curves
  crv <- getCurves(lin)
  plot(rd, col = pal[g], pch = 16, asp = 1)
  lines(crv, lwd = 2, col = "black")
  
  ## Monocle BEAM analysis ----
  ### Monocle 2 BEAM analysis
  trueWeights <- getWeightsBifurcation(dataset, crv)
  featureInfo <- data.frame(gene_short_name = rownames(counts))
  rownames(featureInfo) <- rownames(counts)
  fd <- new("AnnotatedDataFrame", featureInfo)
  cds <- newCellDataSet(cellData = normCounts, featureData = fd, 
                        expressionFamily = negbinomial.size())
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  cds <- reduceDimension(cds)
  cds <- orderCells(cds)
  set.seed(11)
  branch <- rep(NA, ncol(counts))
  branch[trueWeights[, 1] == 1] <- "A"
  branch[trueWeights[, 2] == 1] <- "B"
  branch[is.na(branch)] <- c("A", "B")[sample(1:2, size = sum(is.na(branch)), replace = TRUE)]
  phenoData(cds)$Branch <- as.factor(branch)
  phenoData(cds)$original_cell_id <- colnames(counts)
  phenoData(cds)$Pseudotime <- truePseudotime
  
  ## tradeSeq  ----
  ### tradeSeq: fit smoothers on truth data
  trueT <- matrix(truePseudotime, nrow = length(truePseudotime), ncol = 2, byrow = FALSE)
  
  ## edgeR ---
  edgeR <- function(){
    clF <- as.factor(cl)
    design <- model.matrix(~clF)
    d <- DGEList(counts)
    d <- calcNormFactors(d)
    d <- estimateDisp(d, design)
    fit <- glmFit(d, design)
    L <- matrix(0, nrow = ncol(fit$coefficients), ncol = 1)
    rownames(L) <- colnames(fit$coefficients)
    endClusters <- c(2, 4)
    if (1 %in% endClusters) {
      not1Cluster <- endClusters[!endClusters == 1]
      L[not1Cluster,1] <- 1
    } else {
      L[endClusters,1] <- c(1,-1)
    }
    lrt <- glmLRT(fit, contrast = L)
  }
  
  ## Benchmark time ----
  time_benchmark <- microbenchmark(
    fitGAM(as.matrix(counts), pseudotime = trueT, cellWeights = trueWeights),
    BEAM_kvdb(cds, cores = 1),
    edgeR,
    times = 1L
  )
  save
  
  ## Benchmark memory ----
  ### tradeSeq
  if (!file.exists(here("simulation", "time", "fitGam-memory.Rprof"))) {
    profvis(fitGAM(as.matrix(counts), pseudotime = trueT, cellWeights = trueWeights),
            prof_output = here("simulation", "time",
                               paste0(size, "-fitGam-memory.Rprof")))
  }
  profvis(prof_input = here("simulation", "time",
                            paste0(size, "-fitGam-memory.Rprof")))
  
  ### BEAM
  if (!file.exists(here("simulation", "time", "BEAM-memory.Rprof"))) {
    profvis(BEAM_kvdb(cds, cores = 1),
            prof_output = here("simulation", "time",
                               paste0(size, "-BEAM-memory.Rprof")))
  }
  profvis(prof_input = here("simulation", "time",
                            paste0(size, "-BEAM-memory.Rprof")))
  
  ### edgeR
  if (!file.exists(here("simulation", "time", "edgeR-memory.Rprof"))) {
    profvis(BEAM_kvdb(cds, cores = 1),
            prof_output = here("simulation", "time",
                               paste0(size, "-edgeR-memory.Rprof")))
  }
  profvis(prof_input = here("simulation", "time",
                            paste0(size, "-edgeR-memory.Rprof")))
}
