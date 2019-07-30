library(monocle)
library(DESeq2)
library(RColorBrewer)
library(mgcv)
library(slingshot)
library(tradeSeq)
library(microbenchmark)
library(edgeR)
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

dataset <- readRDS(paste0(here("simulation", "time", "bigDyntoyDataset.rds")))
counts <- t(dataset$counts)
falseGenes <- dataset$tde_overall$feature_id[dataset$tde_overall$differentially_expressed]
nullGenes <- dataset$tde_overall$feature_id[!dataset$tde_overall$differentially_expressed]

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
nClusters <- datasetClusters$nClusters[datasetIter]
set.seed(5)
cl <- kmeans(rd, centers = nClusters)$cluster
rafalib::mypar(mfrow = c(1, 2))
plot(rd, col = wes_palette("Darjeeling1", 10, type = "continuous")[cl], pch = 16, asp = 1)
legend("topleft", legend = as.character(1:nClusters),
       col = wes_palette("Darjeeling1", 10, type = "continuous"), pch = 16, cex = 2 / 3,
       bty = "n")
plot(rd, col = brewer.pal(8, "Dark2")[as.numeric(as.factor(gid$group_id))],
     pch = 16, asp = 1)
legend("topright", paste0("M", 1:length(unique(gid$group_id))), col = 1:4, pch = 16)
rafalib::mypar(mfrow = c(1, 3))
plot(rd, col = pal[cl], pch = 16, asp = 1)
legend("topleft", legend = as.character(1:nClusters), col = pal, pch = 16, cex = 2 / 3,
       bty = "n")
plot(rd, col = pal[g], pch = 16, asp = 1)
plot(rd, col = as.numeric(as.factor(gid$group_id)) + 1, pch = 16, asp = 1)

# lineages
lin <- getLineages(rd, cl, start.clus = datasetClusters$start[datasetIter],
                   end.clus = c(datasetClusters$end1[datasetIter],
                                datasetClusters$end2[datasetIter]))
plot(rd, col = pal[g], pch = 16, asp = 1)
lines(lin, lwd = 2)
# curves
crv <- getCurves(lin)
plot(rd, col = pal[g], pch = 16, asp = 1)
lines(crv, lwd = 2, col = "black")

## Monocle BEAM analysis ----
### Monocle 2 BEAM analysis
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
trueWeights <- getWeightsBifurcation(dataset, crv)
trueT <- matrix(truePseudotime, nrow = length(truePseudotime), ncol = 2, byrow = FALSE)
gamListTruth <- fitGAM(counts, pseudotime = trueT, cellWeights = trueWeights)

endRes <- diffEndTest(gamListTruth)
patternRes <- patternTest(gamListTruth)

## Benchmark time ----
microbenchmark(
  fitGAM(counts, pseudotime = trueT, cellWeights = trueWeights),
  BEAM_kvdb(cds, cores = 1),
  times = 10L
)

## Benchmark memory ----
profvis(fitGAM(counts, pseudotime = trueT, cellWeights = trueWeights))
profvis(BEAM_kvdb(cds, cores = 1))