library(monocle) # Load Monocle
#library(reticulate)
#use_python("/Applications/miniconda3/bin/python3.6")

library(rafalib)

cds <- readRDS(gzcon(url("http://trapnell-lab.gs.washington.edu/public_share/valid_subset_GSE72857_cds2.RDS")))

# Update the old CDS object to be compatible with Monocle 3
cds <- updateCDS(cds)
pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$cluster),
                                        c("1" = 'Erythrocyte',
                                        "2" = 'Erythrocyte',
                                        "3" = 'Erythrocyte',
                                        "4" = 'Erythrocyte',
                                        "5" = 'Erythrocyte',
                                        "6" = 'Erythrocyte',
                                        "7" = 'Multipotent progenitors',
                                        "8" = 'Megakaryocytes',
                                        "9" = 'GMP',
                                        "10" = 'GMP',
                                        "11" = 'Dendritic cells',
                                        "12" = 'Basophils',
                                        "13" = 'Basophils',
                                        "14" = 'Monocytes',
                                        "15" = 'Monocytes',
                                        "16" = 'Neutrophils',
                                        "17" = 'Neutrophils',
                                        "18" = 'Eosinophls',
                                        "19" = 'lymphoid'))

cell_type_color <- c("Basophils" = "#E088B8",
                    "Dendritic cells" = "#46C7EF",
                    "Eosinophls" = "#EFAD1E",
                    "Erythrocyte" = "#8CB3DF",
                    "Monocytes" = "#53C0AD",
                    "Multipotent progenitors" = "#4EB859",
                    "GMP" = "#D097C4",
                    "Megakaryocytes" = "#ACC436",
                    "Neutrophils" = "#F5918A",
                    'NA' = '#000080')

DelayedArray:::set_verbose_block_processing(TRUE)

# Passing a higher value will make some computations faster but use more memory. Adjust with caution!
options(DelayedArray.block.size=1000e6)

### added by kvdb: remove cells not part of trajectory
#### remove dendritic cells
table(phenoData(cds)$cell_type2)
cds <- cds[,!phenoData(cds)$cell_type2=="Dendritic cells"]
#### remove Eosinophls
cds <- cds[,!phenoData(cds)$cell_type2=="Eosinophls"]

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- preprocessCDS(cds, num_dim = 20)

cds <- reduceDimension(cds, reduction_method = 'UMAP')#, python_home="/Applications/miniconda3/bin")
cds <- partitionCells(cds)
cds <- learnGraph(cds,  RGE_method = 'SimplePPT')
# note that plot is different from vignette: we dont find branching in erythrocytes
plot_cell_trajectory(cds,
                     color_by = "cell_type2") +
                     scale_color_manual(values = cell_type_color)


## fit trajectory with slingshot
### get UMAP coordinates
x=1
y=2
theta=0
#reduced_dim_coords <- reducedDimK(cds)
S_matrix <- reducedDimS(cds)
 data_df <- data.frame(t(S_matrix[c(x, y), ]))
    colnames(data_df) <- c("data_dim_1", "data_dim_2")
    data_df$sample_name <- row.names(data_df)
  #  data_df <- merge(data_df, lib_info_with_pseudo, by.x = "sample_name",
  #      by.y = "row.names")
    return_rotation_mat <- function(theta) {
        theta <- theta/180 * pi
        matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)),
            nrow = 2)
    }
    rot_mat <- return_rotation_mat(theta)
    cn1 <- c("data_dim_1", "data_dim_2")
    # cn2 <- c("source_prin_graph_dim_1", "source_prin_graph_dim_2")
    # cn3 <- c("target_prin_graph_dim_1", "target_prin_graph_dim_2")
    data_df[, cn1] <- as.matrix(data_df[, cn1]) %*% t(rot_mat)
plot(data_df[,1],data_df[,2], col=cell_type_color[phenoData(cds)$cell_type2], pch=16)


### slingshot
library(RColorBrewer)
gcolpal <- c(brewer.pal(8,"Dark2")[-c(2,3,5)],brewer.pal(12,"Paired")[c(1,2,8,10,9)],brewer.pal(12,"Set3")[c(7,8,12)], brewer.pal(8, "Pastel2")[8], brewer.pal(11,"BrBG")[11], brewer.pal(11,"PiYG")[1], "cyan", "darkblue","darkorchid2", "brown1", "springgreen1", "deepskyblue4", "darkolivegreen","antiquewhite2")

set.seed(97)
rd <- data_df[,1:2]
cl <- kmeans(rd, centers = 7)$cluster
plot(rd, col = brewer.pal(9,"Set1")[cl], pch=16, asp = 1)
library(slingshot)
lin <- getLineages(rd, clusterLabels=cl, start.clus=4)
plot(rd,col=gcolpal[cl], xlab="UMAP1", ylab="UMAP2")
lines(lin, lwd=2)
crv <- getCurves(lin)
plot(rd,col=gcolpal[cl], main="color by cluster", xlab="UMAP1", ylab="UMAP2")
lines(crv, lwd=2)
plot(rd,col=cell_type_color[phenoData(cds)$cell_type2], main="color by cell type", xlab="UMAP1", ylab="UMAP2", pch=16)
lines(crv, lwd=2)

######## tradeR analysis
library(mgcv)
library(tradeR)
counts=exprs(cds)
#gamListPaul <- fitGAM(counts, pseudotime=slingPseudotime(crv,na=FALSE), cellWeights=slingCurveWeights(crv), verbose=TRUE)
load("~/gamListPaul.rda")
#end point test: 2207 (2266) genes
waldEndResPaul <- diffEndTest(gamListPaul, omnibus=TRUE, pairwise=FALSE)
sum(p.adjust(waldEndResPaul$pvalue,"fdr")<=0.05)
endGenes <- rownames(counts)[which(p.adjust(waldEndResPaul$pvalue,"fdr")<=0.05)]
# pattern test: 2425 genes
patternResPaul <- patternTest(gamListPaul)
sum(p.adjust(patternResPaul$pvalue,"fdr")<=0.05, na.rm=TRUE)
patternGenes <- rownames(counts)[which(p.adjust(patternResPaul$pvalue,"fdr")<=0.05)]
# start point test: 2015 genes
waldStartResPaul <- startVsEndTest(gamListPaul, omnibus=TRUE, pairwise=FALSE)
sum(p.adjust(waldStartResPaul$pvalue,"fdr")<=0.05)
startGenes <- rownames(counts)[which(p.adjust(waldStartResPaul$pvalue,"fdr")<=0.05)]

## plot 6 most significant progenitor genes
library(SummarizedExperiment)
png("/Users/koenvandenberge/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/plots/rdPlotStartGenesPaul.png", units="in", width=12, height=9, res=330)
oStart <- order(waldStartResPaul$waldStat, decreasing=TRUE)
mypar(mfrow=c(2,3))
k=0
while(k<6){
  k=k+1
  cols <- colorRampPalette(c("yellow","red"))(20)
  geneId <- rownames(waldStartResPaul)[oStart[k]]
  g <- Hmisc::cut2(log(assays(se)$counts[geneId,]+1), g=20)
  plot(rd,col=cols[g], main=paste0("color by expression of ",geneId), xlab="UMAP1", ylab="UMAP2", pch=16, cex=2/3)
  lines(crv, lwd=2, col="black")
}
dev.off()

#### combine end point with pattern test
library(tidyverse)
compare <- inner_join(patternResPaul %>% mutate(Gene = rownames(patternResPaul),
                                            pattern = waldStat) %>%
                                     select(Gene, pattern),
                      waldEndResPaul %>% mutate(Gene = rownames(waldEndResPaul),
                                        end = waldStat) %>%
                                 select(Gene, end),
                      by = c("Gene" = "Gene")) %>%
           mutate(transientScore = (min_rank(desc(end)))^2 +
                                   (min_rank(pattern) - 1)^2)

ggplot(compare, aes(x = log(pattern), y = log(end))) +
  geom_point(aes(col = transientScore)) +
  labs(x = "endPointTest Wald Statistic (log scale)",
       y = "patternTest Wald Statistic (log scale)") +
  scale_color_continuous(low = "yellow", high = "red") +
  theme_classic()

topTransient <- (compare %>% arrange(desc(transientScore)))
# first gene is a target of Irf8, according to Paul et al that refers to Marquis et al 2011
# second gene Apoe is very important for haematopoiesis and atherosclerosis: https://www.ncbi.nlm.nih.gov/pubmed/21968112, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4824305/
# third gene also target of Irf8 according to Paul et al that refers to Shin et al
# fourth gene, lamp1, is also target of Irf8 (Paul et al) and encodes for hematopoietic stem cell differentiation: https://www.genecards.org/cgi-bin/carddisp.pl?gene=LAMP1



### old monocle: BEAM analysis
detach("package:monocle", unload=TRUE)
library(monocle,lib.loc="~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/methodsPaper/")
DelayedArray:::set_verbose_block_processing(TRUE)
# Passing a higher value will make some computations faster but use more memory. Adjust with caution!
options(DelayedArray.block.size=1000e6)
cds <- readRDS(gzcon(url("http://trapnell-lab.gs.washington.edu/public_share/valid_subset_GSE72857_cds2.RDS")))
pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$cluster),
                                        c("1" = 'Erythrocyte',
                                        "2" = 'Erythrocyte',
                                        "3" = 'Erythrocyte',
                                        "4" = 'Erythrocyte',
                                        "5" = 'Erythrocyte',
                                        "6" = 'Erythrocyte',
                                        "7" = 'Multipotent progenitors',
                                        "8" = 'Megakaryocytes',
                                        "9" = 'GMP',
                                        "10" = 'GMP',
                                        "11" = 'Dendritic cells',
                                        "12" = 'Basophils',
                                        "13" = 'Basophils',
                                        "14" = 'Monocytes',
                                        "15" = 'Monocytes',
                                        "16" = 'Neutrophils',
                                        "17" = 'Neutrophils',
                                        "18" = 'Eosinophls',
                                        "19" = 'lymphoid'))

cell_type_color <- c("Basophils" = "#E088B8",
                    "Dendritic cells" = "#46C7EF",
                    "Eosinophls" = "#EFAD1E",
                    "Erythrocyte" = "#8CB3DF",
                    "Monocytes" = "#53C0AD",
                    "Multipotent progenitors" = "#4EB859",
                    "GMP" = "#D097C4",
                    "Megakaryocytes" = "#ACC436",
                    "Neutrophils" = "#F5918A",
                    'NA' = '#000080')
cds <- cds[,!phenoData(cds)$cell_type2=="Dendritic cells"]
#### remove Eosinophls
cds <- cds[,!phenoData(cds)$cell_type2=="Eosinophls"]

### Monocle BEAM analysis
counts <- exprs(cds)
pd <- phenoData(cds)
fd <- featureData(cds)
cds2 <- newCellDataSet(counts, phenoData=pd, featureData=fd,
                expressionFamily=VGAM::negbinomial.size())
cds2 <- estimateSizeFactors(cds2)
cds2 <- estimateDispersions(cds2)
cds2 <- reduceDimension(cds2, max_components = 2, method = 'ICA')
cds2 <- orderCells(cds2, num_paths=2)
plot_cell_trajectory(cds2, color_by = "State")
BEAM_res <- BEAM(cds2,  cores = 1)
sum(BEAM_res$qval<0.05)
#save(BEAM_res,file="~/BEAM_resMonocleOldPaulEtal.rda")


### get ICA coordinates and slingshot for plot
library(dplyr)
x=1
y=2
theta=0
S_matrixICA <- reducedDimS(cds2)
plot(t(S_matrixICA[1:2,]), col=cell_type_color[phenoData(cds2)$cell_type2], pch=16)

### slingshot
library(RColorBrewer)
gcolpal <- c(brewer.pal(8,"Dark2")[-c(2,3,5)],brewer.pal(12,"Paired")[c(1,2,8,10,9)],brewer.pal(12,"Set3")[c(7,8,12)], brewer.pal(8, "Pastel2")[8], brewer.pal(11,"BrBG")[11], brewer.pal(11,"PiYG")[1], "cyan", "darkblue","darkorchid2", "brown1", "springgreen1", "deepskyblue4", "darkolivegreen","antiquewhite2")

set.seed(97)
data_df <- t(S_matrixICA[1:2,])
rdICA <- data_df
clICA <- kmeans(rdICA, centers = 7)$cluster
plot(rdICA, col = brewer.pal(9,"Set1")[clICA], pch=16, asp = 1) ; legend("bottomright", as.character(1:7), col=brewer.pal(9,"Set1")[1:7], pch=16, bty="n", cex=2/3)
library(slingshot)
linICA <- getLineages(rdICA, clusterLabels=clICA, start.clus=4)
plot(rdICA,col=gcolpal[clICA], xlab="ICA1", ylab="ICA2")
lines(linICA, lwd=2)
crvICA <- getCurves(linICA)
plot(rdICA,ICAcol=gcolpal[clICA], main="color by cluster", xlab="ICA1", ylab="ICA2")
lines(crvICA, lwd=2)
plot(rdICA,col=cell_type_color[phenoData(cds2)$cell_type2], main="color by cell type", xlab="ICA1", ylab="ICA2", pch=16)
lines(crvICA, lwd=2)




# plot_cell_trajectory(cds, color_by = "State")
# cds <- orderCells(cds, root_state=3 , num_paths=2)
# plot_cell_trajectory(cds, color_by = "cell_type2")
# BEAM_res <- BEAM(cds,  cores = 1)
# sum(BEAM_res$qval<0.05)
# load("~/BEAM_resMonocleOldPaulEtal.rda") #BEAM_res object

#### edgeR analysis
# # edgeR with cluster labels
# library(edgeR)
# d <- DGEList(exprs(cds))
# d <- calcNormFactors(d)
# clF <- as.factor(cl)
# clF <- relevel(clF,ref=4) #set progenitor as ref
# design <- model.matrix(~clF)
# d <- estimateDisp(d, design)
# plotBCV(d)
# fit <- glmFit(d, design)
#
# #leukocyte clusters: 4 6 3 5 7
# #erythrocyte clusters: 4 2 1
# Lleuk <- matrix(0,nrow=ncol(fit$coefficients), ncol=8)
# rownames(Lleuk) <- colnames(fit$coefficients)
# Lleuk[c("clF6","clF2"),1] <- c(1,-1)
# Lleuk[c("clF6","clF1"),2] <- c(1,-1)
# Lleuk[c("clF3","clF2"),3] <- c(1,-1)
# Lleuk[c("clF3","clF1"),4] <- c(1,-1)
# Lleuk[c("clF5","clF2"),5] <- c(1,-1)
# Lleuk[c("clF5","clF1"),6] <- c(1,-1)
# Lleuk[c("clF7","clF2"),7] <- c(1,-1)
# Lleuk[c("clF7","clF1"),8] <- c(1,-1)
# lrtLeuk <- glmLRT(fit,contrast=Lleuk)
# sum(p.adjust(lrtLeuk$table$PValue,"fdr")<=0.05)

# edgeR with cell type labels
library(edgeR)
d <- DGEList(exprs(cds))
d <- calcNormFactors(d)
#clF <- as.factor(cl)
#clF <- relevel(clF,ref=4) #set progenitor as ref
ct <- as.factor(phenoData(cds)$cell_type2)
ct <- relevel(ct, ref="Multipotent progenitors")
design <- model.matrix(~ct)
d <- estimateDisp(d, design)
plotBCV(d)
fit <- glmFit(d, design)

#leukocyte clusters: 4 6 3 5 7
#erythrocyte clusters: 4 2 1
Lleuk <- matrix(0,nrow=ncol(fit$coefficients), ncol=5)
rownames(Lleuk) <- colnames(fit$coefficients)
Lleuk[c("ctErythrocyte","ctBasophils"),1] <- c(1,-1)
Lleuk[c("ctErythrocyte","ctGMP"),2] <- c(1,-1)
Lleuk[c("ctErythrocyte","ctMegakaryocytes"),3] <- c(1,-1)
Lleuk[c("ctErythrocyte","ctMonocytes"),4] <- c(1,-1)
Lleuk[c("ctErythrocyte","ctNeutrophils"),5] <- c(1,-1)
lrtLeuk <- glmLRT(fit,contrast=Lleuk)
sum(p.adjust(lrtLeuk$table$PValue,"fdr")<=0.05)


### gene set enrichment
## use data from https://www.sciencedirect.com/science/article/pii/S221367111630131X?via%3Dihub#app3 to perform specific GSEA
library(openxlsx)
gs <- read.xlsx("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/1-s2.0-S221367111630131X-mmc6.xlsx")
#gs <- gs[gs$Lineage%in%c("Neutrophil Lineage", "Erythrocyte Lineage"),]
gs <- gs[gs$GeneSymbol%in%rownames(counts),]
table(gs$Lineage)

gsList <- c()
for(i in 1:length(unique(gs$Lineage))){
  gsList[[i]] <- gs$GeneSymbol[gs$Lineage==unique(gs$Lineage)[i]]
}
names(gsList) <- unique(gs$Lineage)

# gsList <- list()
# gsList[[1]] <- gs$GeneSymbol[gs$Lineage=="Neutrophil Lineage"]
# gsList[[2]] <- gs$GeneSymbol[gs$Lineage=="Erythrocyte Lineage"]
# names(gsList) <- c("Neutrophil Lineage", "Erythrocyte Lineage")

## use ranks of genes as input for fgsea.
library(fgsea)
ranksTradeR <- rank(waldEndResPaul$waldStat)
names(ranksTradeR) <- rownames(waldEndResPaul)
gseaTrader <- fgsea(gsList, ranksTradeR, nperm=1e5, minSize=5)
gseaTrader

pvalBeam <- BEAM_res$pval
names(pvalBeam) <- rownames(BEAM_res)
ranksBEAM <- rank(pvalBeam)
gseaBEAM <- fgsea(gsList, ranksBEAM, nperm=1e4, minSize=5)
gseaBEAM

ranksEdgeR <- rank(lrtLeuk$table$LR)
names(ranksEdgeR) <- rownames(lrtLeuk$table)
gseaEdgeR <- fgsea(gsList, ranksEdgeR, nperm=1e5, minSize=5)
gseaEdgeR

## GSEA for erythrocytes is significant for tradeR, while it isn't for BEAM, and this is the biological contrast we are actually looking at. None of the other gene sets are significant.

########### Cluster significant pattern genes
library(clusterExperiment)
.getPredictRangeDf <- function(m, lineageId, nPoints=100){
  data <- m$model
  vars <- m$model[1, ]
  vars <- vars[!colnames(vars) %in% "y"]
  offsetId <- grep(x = colnames(vars), pattern = "offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
  names(vars)[offsetId] <- offsetName
  # set all times on 0
  vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
  # set all lineages on 0
  vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
  # duplicate to nPoints
  vars <- rbind(vars, vars[rep(1, nPoints - 1), ])
  # set range of pseudotime for lineage of interest
  lineageData <- data[data[, paste0("l", lineageId)] == 1,
                      paste0("t", lineageId)]
  vars[, paste0("t", lineageId)] <- seq(min(lineageData),
                                        max(lineageData),
                                        length = nPoints)
  # set lineage
  vars[, paste0("l", lineageId)] <- 1
  # set offset
  vars[, offsetName] <- mean(m$model[, grep(x = colnames(m$model),
                                            pattern = "offset")])
  return(vars)
}

# transpose to cluster genes.
#yHatOnlyPat <- do.call(rbind,lapply(gamListPaul[patternGenes], predict, type="link"))
#yhatOnlyPatScaled <- t(scale(t(yHatOnlyPat)))
# dont cluster the fitted values but plot the 100point smoothers
nPoints=100
df1 <- .getPredictRangeDf(gamListPaul[[1]], 1, nPoints = nPoints)
df2 <- .getPredictRangeDf(gamListPaul[[1]], 2, nPoints = nPoints)
y1 <- do.call(rbind,lapply(gamListPaul[patternGenes], predict, newdata=df1, type="link"))
y2 <- do.call(rbind,lapply(gamListPaul[patternGenes], predict, newdata=df2, type="link"))
yhatPat <- cbind(y1,y2)
yhatPatScaled <- t(scale(t(yhatPat)))

rsec <- RSEC(t(yhatPatScaled[1:500,]), isCount = FALSE,
    reduceMethod="PCA", nReducedDims=10,combineMinSize=6,
    ncores=2, random.seed=176201, verbose=TRUE)

#rsec <- RSEC(t(yhatPatScaled), isCount = FALSE,
#    reduceMethod="PCA", nReducedDims=50,combineMinSize=10,
#    ncores=1, random.seed=176201, verbose=TRUE)
clusterLabels <- primaryCluster(rsec)
head(clusterMatrix(rsec)[,1:5])
plotClusters(rsec)
plotCoClustering(rsec)
## note that features and cells have been switched here!
# so, I set clusterFeaturesData=FALSE to not cluster the cells and order along pseudotime.
# I set clusterSamplesData to TRUE to cluster the genes.
## first take the mean across genes within a cluster!
plotHeatmap(rsec, clusterFeaturesData=FALSE)

mypar(mfrow=c(3,3), bty='l')
cUniq <- unique(clusterLabels)
cUniq <- cUniq[!cUniq==-1]
for(xx in cUniq){
  cId <- which(clusterLabels==xx)
  plot(x=1:100,y=rep(range(yhatPatScaled[cId,]),50), type="n", main=paste0("Cluster ",xx), xlab="Pseudotime", ylab="Normalized expression")
  for(ii in 1:length(cId)){
    geneId <- rownames(yhatPatScaled)[cId[ii]]
    yhatGene <- yhatPatScaled[geneId,]
    lines(x=1:100, y=yhatGene[1:100], col="orange", lwd=2)
    lines(x=1:100, y=yhatGene[101:200], col="darkseagreen3", lwd=2)
  }
}




## custom plot functions
plotSmoothersIk <- function(m, nPoints = 100, ...){

  data <- m$model
  y <- data$y

  #construct time variable based on cell assignments.
  nCurves <- length(m$smooth)
  timeAll <- c()
  col <- c()
  for (jj in seq_len(nCurves)) {
    for (ii in 1:nrow(data)) {
      if (data[ii, paste0("l", jj)] == 1) {
        timeAll[ii] <- data[ii, paste0("t", jj)]
        col[ii] <- jj
      } else {
        next
      }
    }
  }

  # plot raw data
  plot(x = timeAll, y = log(y + 1), col = alpha(col,.5), pch = 16, cex = 1 / 3,
       ylab = " log(count + 1)", xlab = "Pseudotime", ...)

  #predict and plot smoothers across the range
  for (jj in seq_len(nCurves)) {
    df <- .getPredictRangeDf(m, jj, nPoints = nPoints)
  yhat <- predict(m, newdata = df, type = "response")
  lines(x = df[, paste0("t", jj)], y = log(yhat + 1), col = jj, lwd = 3)
  }
  #legend("topleft", paste0("lineage", seq_len(nCurves)),col = seq_len(nCurves),
  #       lty = 1, lwd = 2, bty = "n", cex = 2 / 3)
}

###### figure for paper
library(scales)
png("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/plots/figCasePaul_v2.png", width=9, height=6, units="in", res=200)
rafalib::mypar()
# layout(matrix(c(1,1, 2,3, 8,9,
#               1,1, 4,5, 10,11,
#               1,1, 6,7, 12,13), nrow=3, ncol=6, byrow=TRUE))
layout(matrix(c(1,1, 3,4, 9,10,
                1,1, 3,4, 9,10,
                1,1, 5,6, 11,12,
                2,2, 5,6, 11,12,
                2,2, 7,8, 13,14,
                2,2, 7,8, 13,14), nrow=6, ncol=6, byrow=TRUE))
par(cex.lab=1.25, cex.axis=1.1)

### panel A: ICA trajectory with cell types
plot(rdICA, col=alpha(cell_type_color[phenoData(cds)$cell_type2],.8), pch=16, xlab="ICA1", ylab="ICA2")
lines(crvICA@curves$curve1$s[order(crvICA@curves$curve1$lambda),], col="orange", lwd=4)
lines(crvICA@curves$curve2$s[order(crvICA@curves$curve2$lambda),], col="darkseagreen3", lwd=4)
legend("topright", unique(phenoData(cds)$cell_type2),
       col = cell_type_color[unique(phenoData(cds)$cell_type2)],
       pch = 16, cex = 4/5, bty = "n")
legend("bottomright", c("Leukocyte lineage", "Erythrocyte lineage"),
       col = c("orange", "darkseagreen3"),
      lwd=3, lty=1, bty = "n", cex=4/5)
mtext("a", at=-13, font=2, cex=4/3)

### panel B: UMAP trajectory with cell types
plot(rd,col=alpha(cell_type_color[phenoData(cds)$cell_type2],.8), pch=16, xlab="UMAP1", ylab="UMAP2")
lines(crv@curves$curve1$s[order(crv@curves$curve1$lambda),], col="orange", lwd=4)
lines(crv@curves$curve2$s[order(crv@curves$curve2$lambda),], col="darkseagreen3", lwd=4)
legend("topright", unique(phenoData(cds)$cell_type2),
       col = cell_type_color[unique(phenoData(cds)$cell_type2)],
       pch = 16, cex = 4/5, bty = "n")
legend("bottomleft", c("Leukocyte lineage", "Erythrocyte lineage"),
       col = c("orange", "darkseagreen3"),
      lwd=3, lty=1, bty = "n", cex=4/5)
mtext("b", at=-0.21, font=2, cex=4/3)

### panel C: four interesting genes involved in heamatopoiesis.
palette(c("orange","darkseagreen3"))
plotSmoothersIk(gamListPaul[["Prtn3"]], main="Prtn3")
mtext("c", at=-0.5, font=2, cex=4/3)
plotSmoothersIk(gamListPaul[["Mpo"]], main="Mpo", ylim=c(0,5))
plotSmoothersIk(gamListPaul[["Car2"]], main="Car2", ylim=c(0,5))
plotSmoothersIk(gamListPaul[["Ctsg"]], main="Ctsg", ylim=c(0,5))
plotSmoothersIk(gamListPaul[["Elane"]], main="Elane", ylim=c(0,5))
plotSmoothersIk(gamListPaul[["Car1"]], main="Car1", ylim=c(0,5))

### panel D: clusters of gene families
#mypar(mfrow=c(3,2), bty='l')
for(xx in 1:6){
  cId <- which(clusterLabels==xx)
  plot(x=1:100,y=rep(range(yhatPatScaled[cId,]),50), type="n", main=paste0("Cluster ",xx), xlab="Pseudotime", ylab="Normalized expression", ylim=c(-3,2.5))
  if(xx==1) mtext("d", at=-40, font=2, cex=4/3)
  for(ii in 1:length(cId)){
    geneId <- rownames(yhatPatScaled)[cId[ii]]
    yhatGene <- yhatPatScaled[geneId,]
    lines(x=1:100, y=yhatGene[1:100], col="orange", lwd=2)
    lines(x=1:100, y=yhatGene[101:200], col="darkseagreen3", lwd=2)
  }
}
dev.off()


############################################
#### cluster-based vs continuous DE #######
############################################
# clusters contain a mixture of cell types
table(phenoData(cds)$cell_type2,cl[sampleNames(phenoData(cds))])
# this leads to high within-cluster variability
# edgeR
library(edgeR)
d <- DGEList(exprs(cds))
d <- calcNormFactors(d)
cl <- as.factor(cl)
cl <- relevel(cl,ref=4) #set progenitor as ref
design <- model.matrix(~cl)
d <- estimateDisp(d, design)
plotBCV(d)
fit <- glmFit(d, design)

# plot clusters
plot(rd,col=gcolpal[cl], main="color by cluster", xlab="UMAP1", ylab="UMAP2", pch=16) ; legend("topright",as.character(1:7),col=gcolpal[1:7], pch=16)
lines(crv, lwd=2)

# check for association in leukocyte lineage by cluster-based DE
#leukocyte clusters: 4 6 3 5 7
Lleuk <- matrix(0,nrow=ncol(fit$coefficients), ncol=10)
rownames(Lleuk) <- colnames(fit$coefficients)
Lleuk["cl6",1] <- 1
Lleuk["cl3",2] <- 1
Lleuk["cl5",3] <- 1
Lleuk["cl7",4] <- 1
Lleuk[c("cl6","cl3"),5] <- c(1,-1)
Lleuk[c("cl6","cl5"),6] <- c(1,-1)
Lleuk[c("cl6","cl7"),7] <- c(1,-1)
Lleuk[c("cl3","cl5"),8] <- c(1,-1)
Lleuk[c("cl3","cl7"),9] <- c(1,-1)
Lleuk[c("cl5","cl7"),10] <- c(1,-1)

lrtLeuk <- glmLRT(fit,contrast=Lleuk)

### genes identified with tradeR but not with edgeR
deEdgeRLeuk <- rownames(lrtLeuk)[p.adjust(lrtLeuk$table$PValue,"fdr")<=0.05]
nullEdgeRLeuk <- rownames(lrtLeuk)[!rownames(lrtLeuk)%in%deEdgeRLeuk]

#pSmoothTradeR <- getSmootherPvalues(gamListPaul)
#statSmoothTradeR <- getSmootherTestStats(gamListPaul)
pSmoothTradeRIct <- getSmootherPvalues(gamListPaulIct)
statSmoothTradeRIct <- getSmootherTestStats(gamListPaulIct)
deTradeRLeuk <- rownames(pSmoothTradeRIct)[p.adjust(pSmoothTradeRIct[,1],"fdr")<=0.05]

onlyTraderLeuk <- deTradeRLeuk[!deTradeRLeuk%in%deEdgeRLeuk]

i=0

mypar(mfrow=c(1,4))
i=i+1 ; plotSmoothers(gamListPaul[[onlyTraderLeuk[i]]]) ; plot(gamListPaul[[onlyTraderLeuk[i]]]) ; plotGeneCount(counts=counts, curve=crv, gene=onlyTraderLeuk[i], rd=rd) ; onlyTraderLeuk[i]

write.table(onlyTraderLeuk,file="~/onlyTraderLeuk.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


# i=2, 64: gene with all-zero that's significantly associated with the lineage...

# "4930486L24Rik" has nearly all zero but still very high test statistics in smoother association
# things like gene Grhpr should not be significant in a smoothe association. are they significant because the mean of the smoother differs from the overall mean?

# interesting: i=98 (Clec10a),

# erythrocyte lineage
#erythrocyte clusters: 4 2 1
Ler <- matrix(0,nrow=ncol(fit$coefficients), ncol=3)
rownames(Ler) <- colnames(fit$coefficients)
Ler["cl2",1] <- 1
Ler["cl1",2] <- 1
Ler[c("cl2","cl1"),3] <- c(1,-1)
lrtEr <- glmLRT(fit,contrast=Ler)
deEdgeREr <- rownames(lrtEr)[p.adjust(lrtEr$table$PValue,"fdr")<=0.05]
deTradeREr <- rownames(pSmoothTradeRIct)[p.adjust(pSmoothTradeRIct[,2],"fdr")<=0.05]
onlyTraderEr <- deTradeREr[!deTradeREr%in%deEdgeREr]
onlyTraderEr <- onlyTraderEr[order(statSmoothTradeRIct[onlyTraderEr,2], decreasing=TRUE)]



i=0

mypar(mfrow=c(1,4))
i=i+1 ; plotSmoothers(gamListPaul[[onlyTraderEr[i]]]) ; plot(gamListPaul[[onlyTraderEr[i]]]) ; plotGeneCount(counts=counts, curve=crv, gene=onlyTraderEr[i], rd=rd) ; onlyTraderEr[i]

write.table(onlyTraderEr,file="~/onlyTraderEr.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
