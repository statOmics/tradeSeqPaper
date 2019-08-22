.libPaths("/accounts/campus/hector.rouxdebezieux/R/x86_64-pc-linux-gnu-library/3.5")
library(RColorBrewer)
library(mgcv)
library(slingshot)
library(tradeSeq)
library(microbenchmark)
library(here)
library(rafalib)
library(dyno)
library(dyntoy)
library(tidyverse)
library(wesanderson)
library(BiocParallel)
library(doParallel)

## Pre-process ----
NCORES <- 2
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
  pca <- prcomp(log1p(t(normCounts)), scale. = FALSE)
  rd <- pca$x[, 1:4]
  ## cluster
  cl <- as.factor(gid$group_id)
  
  # lineages
  lin <- getLineages(rd, as.numeric(cl), start.clus = 1,
                     end.clus = c(2, 4))
  # curves
  crv <- getCurves(lin)
  
  ## tradeSeq  ----
  ### tradeSeq: fit smoothers on truth data
  trueT <- matrix(truePseudotime, nrow = length(truePseudotime), ncol = 2, byrow = FALSE)
  gamList <- fitGAM(as.matrix(counts), pseudotime = trueT, cellWeights = trueWeights,
                    nknots = 4)
  print("...tradeSeq")  
  ## Benchmark time ----
  print("...benchmark")
  time_benchmark <- microbenchmark(
    diffEndTest(gamList),
    startVsEndTest(gamList),
    patternTest(gamList),
    associationTest(gamList),
    earlyDETest(gamList, knots = c(1, 2)),
    times = 10L
  )
  write.table(x = time_benchmark,
              file = here("simulation", "time", "data",
                          paste0(size, "-time-tests.txt")))
  
  }