library(tradeSeq)
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
options(DelayedArray.block.size = 1000e6)

### added by kvdb: remove cells not part of trajectory
#### remove dendritic cells
table(phenoData(cds)$cell_type2)
cds <- cds[, !phenoData(cds)$cell_type2 == "Dendritic cells"]
#### remove Eosinophls
cds <- cds[, !phenoData(cds)$cell_type2 == "Eosinophls"]

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- preprocessCDS(cds, num_dim = 20)

cds <- reduceDimension(cds, reduction_method = "UMAP") # , python_home="/Applications/miniconda3/bin")
cds <- partitionCells(cds)
cds <- learnGraph(cds, RGE_method = "SimplePPT")
# note that plot is different from vignette: we dont find branching in erythrocytes
plot_cell_trajectory(cds,
                     color_by = "cell_type2") +
                     scale_color_manual(values = cell_type_color) + theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.text=element_text(size=12))


## fit trajectory with slingshot
### get UMAP coordinates
x <- 1
y <- 2
theta <- 0
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
plot(data_df[, 1], data_df[, 2], col = cell_type_color[phenoData(cds)$cell_type2], pch = 16)

### slingshot
library(RColorBrewer)
gcolpal <- c(brewer.pal(8, "Dark2")[-c(2, 3, 5)],
             brewer.pal(12, "Paired")[c(1, 2, 8, 10, 9)],
             brewer.pal(12, "Set3")[c(7, 8, 12)],
             brewer.pal(8, "Pastel2")[8], brewer.pal(11, "BrBG")[11],
             brewer.pal(11, "PiYG")[1], "cyan", "darkblue", "darkorchid2",
             "brown1", "springgreen1", "deepskyblue4", "darkolivegreen",
             "antiquewhite2")

set.seed(97)
rd <- data_df[, 1:2]
cl <- kmeans(rd, centers = 7)$cluster
plot(rd, col = brewer.pal(9, "Set1")[cl], pch = 16, asp = 1)
library(slingshot)
lin <- getLineages(rd, clusterLabels = cl, start.clus = 4)
plot(rd, col = gcolpal[cl], xlab = "UMAP1", ylab = "UMAP2")
lines(lin, lwd = 2)
crv <- getCurves(lin)
plot(rd, col = gcolpal[cl], main = "color by cluster", xlab = "UMAP1", ylab = "UMAP2")
lines(crv, lwd = 2)
plot(rd, col = cell_type_color[phenoData(cds)$cell_type2], main = "color by cell type", xlab = "UMAP1", ylab = "UMAP2", pch = 16)
lines(crv, lwd = 2)

######## tradeSeq analysis
library(mgcv)
library(tradeSeq)
counts <- exprs(cds)
for(seed in 1:10){
  assign(paste0("gamListPaul",seed),fitGAM(counts, pseudotime=slingPseudotime(crv,na=FALSE), cellWeights=slingCurveWeights(crv), verbose=TRUE,seed=seed))
}

resList <- sapply(1:10,function(seed) startVsEndTest(get(paste0("gamListPaul",seed))), simplify=FALSE)
nrDE = unlist(lapply(resList, function(x) sum(p.adjust(x$pvalue,"fdr")<=0.05, na.rm=TRUE)))
nrDE

### intersection of DE genes
library(UpSetR)
deId <- as.data.frame(do.call(cbind,lapply(resList, function(x) as.numeric(p.adjust(x$pvalue,"fdr")<=0.05))))
colnames(deId) <- paste0("iter.",1:10)
deId[is.na(deId)] <- 0
upset(deId, order.by = "freq", text.scale = c(1.3, 1.3, 1, 1, 2, 2), nsets=10, keep.order=TRUE, nintersects=20)

## what about the top 1000 genes?
orders <- as.data.frame(do.call(cbind,lapply(resList, function(x) order(x$waldStat, decreasing=TRUE))))
orders1k <- orders[1:1000,]
# pairwise sharing of 1k top genes
range(apply(combn(colnames(orders1k),2),2, function(cols) mean(orders1k[,cols[1]] %in% orders1k[,cols[2]])))
