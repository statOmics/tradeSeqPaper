
### old monocle: BEAM analysis
# detach("package:monocle", unload=TRUE)
library(monocle, lib.loc = "~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/methodsPaper/")
DelayedArray:::set_verbose_block_processing(TRUE)
# Passing a higher value will make some computations faster but use more memory. Adjust with caution!
options(DelayedArray.block.size = 1000e6)
cds <- readRDS(gzcon(url("http://trapnell-lab.gs.washington.edu/public_share/valid_subset_GSE72857_cds2.RDS")))
pData(cds)$cell_type2 <- plyr::revalue(
  as.character(pData(cds)$cluster),
  c(
    "1" = "Erythrocyte",
    "2" = "Erythrocyte",
    "3" = "Erythrocyte",
    "4" = "Erythrocyte",
    "5" = "Erythrocyte",
    "6" = "Erythrocyte",
    "7" = "Multipotent progenitors",
    "8" = "Megakaryocytes",
    "9" = "GMP",
    "10" = "GMP",
    "11" = "Dendritic cells",
    "12" = "Basophils",
    "13" = "Basophils",
    "14" = "Monocytes",
    "15" = "Monocytes",
    "16" = "Neutrophils",
    "17" = "Neutrophils",
    "18" = "Eosinophls",
    "19" = "lymphoid"
  )
)

cell_type_color <- c(
  "Basophils" = "#E088B8",
  "Dendritic cells" = "#46C7EF",
  "Eosinophls" = "#EFAD1E",
  "Erythrocyte" = "#8CB3DF",
  "Monocytes" = "#53C0AD",
  "Multipotent progenitors" = "#4EB859",
  "GMP" = "#D097C4",
  "Megakaryocytes" = "#ACC436",
  "Neutrophils" = "#F5918A",
  "NA" = "#000080"
)
cds <- cds[, !phenoData(cds)$cell_type2 == "Dendritic cells"]
#### remove Eosinophls
cds <- cds[, !phenoData(cds)$cell_type2 == "Eosinophls"]

### Monocle BEAM analysis on ICA
counts <- exprs(cds)
pd <- phenoData(cds)
fd <- featureData(cds)
cds2 <- newCellDataSet(counts,
  phenoData = pd, featureData = fd,
  expressionFamily = VGAM::negbinomial.size()
)
cds2 <- estimateSizeFactors(cds2)
cds2 <- estimateDispersions(cds2)
cds2 <- reduceDimension(cds2, max_components = 2, method = "ICA")
cds2 <- orderCells(cds2, num_paths = 2)
plot_cell_trajectory(cds2, color_by = "State")
plot_cell_trajectory(cds2, color_by = "cell_type2")
BEAM_res <- BEAM(cds2, cores = 1)
sum(BEAM_res$qval < 0.05)
# save(BEAM_res,file="~/BEAM_resMonocleOldPaulEtal.rda")


## fit trajectory with slingshot
### get ICA coordinates
library(dplyr)
x <- 1
y <- 2
theta <- 0
S_matrix <- reducedDimS(cds2)
plot(t(S_matrix[1:2, ]), col = cell_type_color[phenoData(cds2)$cell_type2], pch = 16)


### slingshot
library(RColorBrewer)
gcolpal <- c(brewer.pal(8, "Dark2")[-c(2, 3, 5)], brewer.pal(12, "Paired")[c(1, 2, 8, 10, 9)], brewer.pal(12, "Set3")[c(7, 8, 12)], brewer.pal(8, "Pastel2")[8], brewer.pal(11, "BrBG")[11], brewer.pal(11, "PiYG")[1], "cyan", "darkblue", "darkorchid2", "brown1", "springgreen1", "deepskyblue4", "darkolivegreen", "antiquewhite2")

set.seed(97)
data_df <- t(S_matrix[1:2, ])
rd <- data_df
cl <- kmeans(rd, centers = 7)$cluster
plot(rd, col = brewer.pal(9, "Set1")[cl], pch = 16, asp = 1)
legend("bottomright", as.character(1:7), col = brewer.pal(9, "Set1")[1:7], pch = 16, bty = "n", cex = 2 / 3)
library(slingshot)
lin <- getLineages(rd, clusterLabels = cl, start.clus = 4)
plot(rd, col = gcolpal[cl], xlab = "ICA1", ylab = "ICA2")
lines(lin, lwd = 2)
crv <- getCurves(lin)
plot(rd, col = gcolpal[cl], main = "color by cluster", xlab = "ICA1", ylab = "ICA2")
lines(crv, lwd = 2)
plot(rd, col = cell_type_color[phenoData(cds)$cell_type2], main = "color by cell type", xlab = "ICA1", ylab = "ICA2", pch = 16)
lines(crv, lwd = 2)

######## tradeSeq analysis
library(mgcv)
library(tradeSeq)
counts <- exprs(cds)
gamListPaulICA <- fitGAM(counts, pseudotime = slingPseudotime(crv, na = FALSE), cellWeights = slingCurveWeights(crv), verbose = TRUE)

# end point test: 2543 (2266) genes
waldEndResPaul <- diffEndTest(gamListPaulICA, omnibus = TRUE, pairwise = FALSE)
sum(p.adjust(waldEndResPaul$pvalue, "fdr") <= 0.05)
endGenes <- rownames(counts)[which(p.adjust(waldEndResPaul$pvalue, "fdr") <= 0.05)]
# pattern test: 2290 genes
patternResPaul <- patternTest(gamListPaulICA)
sum(p.adjust(patternResPaul$pvalue, "fdr") <= 0.05, na.rm = TRUE)
patternGenes <- rownames(counts)[which(p.adjust(patternResPaul$pvalue, "fdr") <= 0.05)]
# start point test: 2442 genes
waldStartResPaul <- startVsEndTest(gamListPaulICA, omnibus = TRUE, pairwise = FALSE)
sum(p.adjust(waldStartResPaul$pvalue, "fdr") <= 0.05)
startGenes <- rownames(counts)[which(p.adjust(waldStartResPaul$pvalue, "fdr") <= 0.05)]


# edgeR with cell type labels
library(edgeR)
d <- DGEList(exprs(cds))
d <- calcNormFactors(d)
# clF <- as.factor(cl)
# clF <- relevel(clF,ref=4) #set progenitor as ref
ct <- as.factor(phenoData(cds)$cell_type2)
ct <- relevel(ct, ref = "Multipotent progenitors")
design <- model.matrix(~ct)
d <- estimateDisp(d, design)
plotBCV(d)
fit <- glmFit(d, design)

# leukocyte clusters: 4 6 3 5 7
# erythrocyte clusters: 4 2 1
Lleuk <- matrix(0, nrow = ncol(fit$coefficients), ncol = 5)
rownames(Lleuk) <- colnames(fit$coefficients)
Lleuk[c("ctErythrocyte", "ctBasophils"), 1] <- c(1, -1)
Lleuk[c("ctErythrocyte", "ctGMP"), 2] <- c(1, -1)
Lleuk[c("ctErythrocyte", "ctMegakaryocytes"), 3] <- c(1, -1)
Lleuk[c("ctErythrocyte", "ctMonocytes"), 4] <- c(1, -1)
Lleuk[c("ctErythrocyte", "ctNeutrophils"), 5] <- c(1, -1)
lrtLeuk <- glmLRT(fit, contrast = Lleuk)
sum(p.adjust(lrtLeuk$table$PValue, "fdr") <= 0.05)


### gene set enrichment
## use data from https://www.sciencedirect.com/science/article/pii/S221367111630131X?via%3Dihub#app3 to perform specific GSEA
library(openxlsx)
gs <- read.xlsx("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/1-s2.0-S221367111630131X-mmc6.xlsx")
# gs <- gs[gs$Lineage%in%c("Neutrophil Lineage", "Erythrocyte Lineage"),]
gs <- gs[gs$GeneSymbol %in% rownames(counts), ]
table(gs$Lineage)

gsList <- c()
for (i in 1:length(unique(gs$Lineage))) {
  gsList[[i]] <- gs$GeneSymbol[gs$Lineage == unique(gs$Lineage)[i]]
}
names(gsList) <- unique(gs$Lineage)

# gsList <- list()
# gsList[[1]] <- gs$GeneSymbol[gs$Lineage=="Neutrophil Lineage"]
# gsList[[2]] <- gs$GeneSymbol[gs$Lineage=="Erythrocyte Lineage"]
# names(gsList) <- c("Neutrophil Lineage", "Erythrocyte Lineage")

## use ranks of genes as input for fgsea.
library(fgsea)
rankstradeSeq <- rank(waldEndResPaul$waldStat)
names(rankstradeSeq) <- rownames(waldEndResPaul)
gseatradeSeq <- fgsea(gsList, rankstradeSeq, nperm = 1e5, minSize = 5)
gseatradeSeq

pvalBeam <- BEAM_res$pval
names(pvalBeam) <- rownames(BEAM_res)
ranksBEAM <- rank(pvalBeam)
gseaBEAM <- fgsea(gsList, ranksBEAM, nperm = 1e4, minSize = 5)
gseaBEAM

ranksEdgeR <- rank(lrtLeuk$table$LR)
names(ranksEdgeR) <- rownames(lrtLeuk$table)
gseaEdgeR <- fgsea(gsList, ranksEdgeR, nperm = 1e5, minSize = 5)
gseaEdgeR

## GSEA for erythrocytes is significant for tradeSeq, while it isn't for BEAM, and this is the biological contrast we are actually looking at. None of the other gene sets are significant.