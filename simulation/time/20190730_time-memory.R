library(princurve)
library(RColorBrewer)
library(mgcv)
library(BiocParallel)
library(doParallel)
library(microbenchmark)
library(here)
library(tidyverse)
library(slingshot)
library(DESeq2)
library(tradeSeq)
library(rafalib)
library(dyno)
library(dyntoy)
library(monocle)
library(edgeR)
library(wesanderson)

## Pre-process ----
palette(wes_palette("Darjeeling1", 10, type = "continuous"))
source(here::here("simulation", "time", "20190611_helper.R"))
set.seed(87657)

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

print("Checking stack")
print(Cstack_info()["size"] > 7971060)
## datasets ----

for (size in 2:4) {
  ## Generate dataset ----
  print(paste0("dataset of size", 10^size))
  dataset <- readRDS(
    here::here("simulation", "time", "data", paste0(size, "dataset.rds"))
    )
  ## pre-process ----
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
  print("...slingshot")
  ## dim red
  pca <- prcomp(log1p(t(normCounts)), scale. = FALSE, rank. = 4)
  rd <- pca$x[, 1:4]
  ## cluster
  cl <- as.factor(gid$group_id)
  
  # lineages
  lin <- getLineages(rd, as.numeric(cl), start.clus = 1,
                     end.clus = c(2, 4))
  # curves
  crv <- getCurves(lin)
  
  # ## Monocle BEAM analysis ----
  print("...BEAM")
  ### Monocle 2 BEAM analysis
  trueWeights <- getWeightsBifurcation(dataset, crv)
  featureInfo <- data.frame(gene_short_name = rownames(counts))
  rownames(featureInfo) <- rownames(counts)
  fd <- new("AnnotatedDataFrame", featureInfo)
  cds <- newCellDataSet(cellData = normCounts, featureData = fd,
                        expressionFamily = negbinomial.size())
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  cds <- reduceDimension(cds, reduction_method = "ICA")
  cds <- orderCells(cds)
  print("is it here?")
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
  print("...tradeSeq")  
  # ## edgeR ----
  print("...edgeR")
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
  
  ## GPFates ----
  print("...GPfates")
  if (size < 4) {
    logCpm <- edgeR::cpm(normCounts, prior.count = .125, log = TRUE)
    sampleInfo <- data.frame(global_pseudotime = truePseudotime)
    rownames(sampleInfo) <- colnames(counts)
    write.table(logCpm, file = "./timeBenchLogCpm.txt", row.names = TRUE,
                col.names = TRUE, quote = FALSE)
    write.table(sampleInfo, file = "./timeBenchSampleInfo.txt", row.names = TRUE,
                col.names = TRUE, quote = FALSE)
    system("python3 ./20190806_preprocessGPfatesTimeMemBenchmark.py",
           ignore.stdout = TRUE)
  }
  
  ## Benchmark time ----
  print("...benchmark")
  if (size < 4) {
    time_benchmark <- microbenchmark(
      fitGAM(as.matrix(counts), pseudotime = trueT, cellWeights = trueWeights,
             nknots = 4),
      BEAM_kvdb(cds, cores = 1),
      edgeR(),
      system("python3 ./20190806_analyzeGPfatesTimeBenchmark.py",
             ignore.stdout = TRUE),
      times = ifelse(size == 2, 10, 2)
    )  
  } else {
    time_benchmark <- microbenchmark(
      fitGAM(as.matrix(counts), pseudotime = trueT, cellWeights = trueWeights,
             nknots = 4),
      BEAM_kvdb(cds, cores = 1),
      edgeR(),
      times = ifelse(size == 2, 10, 2)
    )
  }
  
  write.table(x = time_benchmark,
              file = here("simulation", "time", "data",
                          paste0(size, "-time-benchmark.txt")))
  print(time_benchmark)
  
  # Benchmark memory ----
  if (size < 4) {
    mem <- rep(0, 4)
    names(mem) <- c("tradeSeq", "BEAM", "edgeR", "GPFates")
  } else {
    mem <- rep(0, 3)
    names(mem) <- c("tradeSeq", "BEAM", "edgeR")
  }
  
  ### tradeSeq
  print("...tradeSeq memory")
  Rprof(filename = here::here("simulation", "time","Rprof.out"),
        memory.profiling = TRUE)
  test <- fitGAM(as.matrix(counts), pseudotime = trueT, cellWeights = trueWeights,
                 verbose = FALSE)
  Rprof(filename = 'NULL')
  mem["tradeSeq"] <-
    summaryRprof(
    filename = here::here("simulation", "time", "Rprof.out"),
    memory = "both")$by.total[, "mem.total"] %>% 
    max()
    
  print(mem)
  
  ### BEAM
  print("...BEAM memory")
  Rprof(filename = here::here("simulation", "time","Rprof.out"),
        memory.profiling = TRUE)
  test <- BEAM_kvdb(cds, cores = 1)
  Rprof(filename = 'NULL')
  mem["BEAM"] <- summaryRprof(
    filename = here::here("simulation", "time", "Rprof.out"),
    memory = "both")$by.total[, "mem.total"] %>%
    max()

  ### edgeR
  print("...edgeR memory")
  Rprof(filename = here::here("simulation", "time","Rprof.out"),
        memory.profiling = TRUE)
  test <- edgeR()
  Rprof(filename = 'NULL')
  mem["edgeR"] <- summaryRprof(
    filename = here::here("simulation", "time", "Rprof.out"),
    memory = "both")$by.total[, "mem.total"] %>%
    max()

  write.table(x = mem,file = here("simulation", "time", "data", 
                                  paste0(size, "-mem-benchmark.txt")))

  ## GPfates
  print("...GPFates memory")
  if (size < 4) {
    memGPfatesAll <-
      system("python3 ./20190806_analyzeGPfatesMemoryBenchmark.py",
             intern = TRUE, ignore.stdout = FALSE)
    mem1 <- sapply(memGPfatesAll, strsplit, split = " [ ]+") %>% unlist()
    mem1 <- str_subset(mem1, "MiB")
    mem1 <- str_remove(mem1, " MiB")
    print(max(mem1))
    mem["GPFates"] <- max(as.numeric(mem1), na.rm = TRUE)
  }
  
  ## All together
  write.table(x = mem,file = here("simulation", "time", "data",
                                  paste0(size, "-mem-benchmark.txt")))
}

file.remove(here::here("simulation", "time", "Rprof.out"))
file.remove(here::here("simulation", "time", "Rprof.out"))
file.remove(here::here("simulation", "time", "timeBenchLogCpm.txt"))
file.remove(here::here("simulation", "time", "timeBenchSampleInfo.txt"))
file.remove(here::here("simulation", "time", "m.pkl"))