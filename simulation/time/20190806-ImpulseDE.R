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
library(ImpulseDE2)
library(tidyverse)
library(dyno)
library(dyntoy)
library(here)
source(here::here("simulation", "time", "20190611_helper.R"))

set.seed(87657)
dataset <- generate_dataset(
  model = model_bifurcating(),
  num_cells = 100,
  num_features = 5000,
  differentially_expressed_rate = .2
)

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
## dim red
pca <- prcomp(log1p(t(normCounts)), scale. = FALSE)
rd <- pca$x[, 1:3]
## cluster
cl <- as.factor(gid$group_id)

# lineages
lin <- getLineages(rd, as.numeric(cl), start.clus = 1,
                   end.clus = c(2, 4))
# curves
crv <- getCurves(lin)
trueWeights <- getWeightsBifurcation(dataset, crv)
trueT <- matrix(truePseudotime, nrow = length(truePseudotime), ncol = 2, byrow = FALSE)

gamModels <- fitGAM(as.matrix(counts)[1:10,], pseudotime = trueT,
                    cellWeights = trueWeights)
m <- gamList[[1]]
branch <- rep(NA, ncol(counts))
time <- rep(NA, ncol(counts))
branch[m$model$l1 == 1] <- "A"
time[m$model$l1 == 1] <- trueT[m$model$l1 == 1, 1]
branch[m$model$l2 == 1] <- "B"
time[m$model$l2 == 1] <- trueT[m$model$l2 == 1, 2]
condition <- ifelse(branch == "A", "control", "case")
# ImpulseDE2 cannot handle datasets where every gene contains at least one zero.
# since then it errors on the DESeq2 size factor estimation. Since the data is
# already quantile normalized, we provide a size factor of 1 as input.
sf <- rep(1, ncol(normCounts)) # data is already normalized
names(sf) <- colnames(normCounts)
# dfAnnotation
dfAnn <- data.frame(Sample = colnames(counts), Condition = condition, Time = time)
# Even with supplied size factors, it errors because it still attempts to calculate
# the size factors for dispersion estimation with DESeq2. Do it manually.
dds <- suppressWarnings(DESeqDataSetFromMatrix(
  countData = round(normCounts),
  colData = dfAnn,
  design = ~ Condition + Condition:Time
))
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- estimateDispersions(dds)
vecDispersionsInv <- mcols(dds)$dispersion
vecDispersions <- 1 / vecDispersionsInv
names(vecDispersions) <- rownames(dds)

# time_benchmark <- microbenchmark(
#   runImpulseDE2(matCountData = round(normCounts), dfAnnotation = dfAnn,
#                 boolCaseCtrl = TRUE, vecSizeFactorsExternal = sf,
#                 vecDispersionsExternal = vecDispersions, scaNProc = 2),
#   times = 10L
# )
# 
# write.table(x = time_benchmark,
#             file = here::here("simulation", "time", "ImpulseDE-time.txt"))

Rprof(here::here("simulation", "time", "small-ImpulseDE-memory.Rprof"),
      memory.profiling = TRUE)
print(system.time(
df <- runImpulseDE2(matCountData = round(normCounts), dfAnnotation = dfAnn,
                      boolCaseCtrl = TRUE, vecSizeFactorsExternal = sf,
                      vecDispersionsExternal = vecDispersions, scaNProc = 2)
))
Rprof(filename = NULL)
print(
  summaryRprof(
    filename = here::here(
      "simulation", "time", "data", "small-ImpulseDE-memory.Rprof"),
    memory = "both")$by.total[, "mem.total"] %>% 
    max(na.rm = TRUE)
  )
