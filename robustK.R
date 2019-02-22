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
gamListPaul10 <- fitGAM(counts, pseudotime=slingPseudotime(crv,na=FALSE), cellWeights=slingCurveWeights(crv), verbose=TRUE)
devExpl10 <- unlist(lapply(gamListPaul10, function(x) summary(x)$dev.expl))
rm(gamListPaul10)

gamListPaul8 <- fitGAM(counts, pseudotime=slingPseudotime(crv,na=FALSE), cellWeights=slingCurveWeights(crv), verbose=TRUE, nknots=8)
devExpl8 <- unlist(lapply(gamListPaul8, function(x) summary(x)$dev.expl))
rm(gamListPaul8)

gamListPaul6 <- fitGAM(counts, pseudotime=slingPseudotime(crv,na=FALSE), cellWeights=slingCurveWeights(crv), verbose=TRUE, nknots=6)
devExpl6 <- unlist(lapply(gamListPaul6, function(x) summary(x)$dev.expl))
rm(gamListPaul6)

gamListPaul12 <- fitGAM(counts, pseudotime=slingPseudotime(crv,na=FALSE), cellWeights=slingCurveWeights(crv), verbose=TRUE, nknots=12)
devExpl12 <- unlist(lapply(gamListPaul12, function(x) summary(x)$dev.expl))
rm(gamListPaul12)

gamListPaul14 <- fitGAM(counts, pseudotime=slingPseudotime(crv,na=FALSE), cellWeights=slingCurveWeights(crv), verbose=TRUE, nknots=14)
devExpl14 <- unlist(lapply(gamListPaul14, function(x) summary(x)$dev.expl))
rm(gamListPaul14)


library(RColorBrewer)
cols=palette(brewer.pal(8,"Dark2"))

par(bty='l')
plot(density(devExpl6), col=cols[1], xlab="% Deviance explained", main="")
lines(density(devExpl8),col=cols[2])
lines(density(devExpl10),col=cols[3])
lines(density(devExpl12),col=cols[4])
lines(density(devExpl14),col=cols[5])
legend("topright",paste0("k=",seq(6,14,by=2)), col=cols[1:5], lty=1, lwd=2)
