setwd("~/Dropbox/research/PhD/research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/")
library(here)
library(slingshot)
library(RColorBrewer)
library(mgcv)
library(tradeSeq)
library(edgeR)
library(rafalib)
library(wesanderson)
library(ggplot2)
library(cowplot)
library(iCOBRA)
library(scales)
library(dplyr)
RNGversion("3.5.0")


### prepare performance plots ####
cols <- c(rep(c("#C6DBEF", "#08306B"), each = 3), "#4292C6", "#4daf4a",
          "#e41a1c", "#e78ac3", "#ff7f00", "darkgoldenrod1")
names(cols) <- c("tradeSeq_slingshot_end", "tradeSeq_GPfates_end", "tradeSeq_Monocle2_end",
                 "tradeSeq_slingshot_pattern", "tradeSeq_GPfates_pattern",
                 "tradeSeq_Monocle2_pattern", "tradeSeq_slingshot_assoc", "Monocle3_assoc",
                 "BEAM", "GPfates", "edgeR", "ImpulseDE2")
linetypes <- c(rep(c("dashed", "dotdash", "solid"), 2), rep("solid", 7))
names(linetypes) <- c("tradeSeq_slingshot_end", "tradeSeq_GPfates_end", "tradeSeq_Monocle2_end",
                      "tradeSeq_slingshot_pattern", "tradeSeq_GPfates_pattern",
                      "tradeSeq_Monocle2_pattern", "tradeSeq_slingshot_assoc", "Monocle3_assoc",
                      "BEAM", "GPfates", "edgeR", "ImpulseDE2")

theme_set(theme_bw())
theme_update(legend.position = "none",
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.grid.major.x = element_line(linetype = "dashed", colour = "black"),
             panel.grid.minor.x = element_line(linetype = "dashed", colour = "grey"),
             axis.title.x = element_text(size = rel(1)),
             axis.title.y = element_text(size = rel(1)),
             axis.text.x = element_text(size = rel(.8)),
             axis.text.y = element_text(size = rel(.8)))

### bifurcating dyntoy plot  ####
dir <- "~/Dropbox/research/PhD/research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/datasets"
cobraFiles <- list.files(dir, pattern="cobra*", full.names=TRUE)
cobraFiles <- cobraFiles[c(1,3:10,2)] #order from 1 to 10

plotPerformanceCurve <- function(cobraObject){
  colnames(cobraObject@pval) <- gsub(colnames(cobraObject@pval),pattern="tradeR",replacement="tradeSeq")
  cobraObject <- calculate_adjp(cobraObject)
  cobraObject <- calculate_performance(cobraObject, binary_truth = "status")

  DyntoyPlot <- data.frame(FDP = cobraObject@fdrtprcurve$FDR,
                           TPR = cobraObject@fdrtprcurve$TPR,
                           method = cobraObject@fdrtprcurve$method,
                         cutoff = cobraObject@fdrtprcurve$CUTOFF)
  pDyntoy <- ggplot(DyntoyPlot, aes(x = FDP, y = TPR, col = method)) +
    geom_path(size = 1, aes(linetype = method)) +
    xlab("FDP") +
    scale_x_continuous(limits = c(0, 1), breaks = c(0.01, 0.05, 0.1, 0.5, 1),
                       minor_breaks = c(0:5) * .1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = cols, breaks = names(cols)) +
    scale_linetype_manual(values = linetypes, breaks = names(linetypes))
    pDyntoy
}

## write source data file
for(ii in 1:length(cobraFiles)){
    cobra <- readRDS(cobraFiles[ii])
    write.table(pval(cobra), file=paste0("~/pvalDataset",ii,".txt"))
}


for(ii in 1:length(cobraFiles)){
    cobra <- readRDS(cobraFiles[ii])
    assign(paste0("bifplot",ii),plotPerformanceCurve(cobra))
}

p1 <- plot_grid(bifplot1, bifplot2, bifplot3, bifplot4, bifplot5,
           bifplot6, bifplot7, bifplot8, bifplot9, bifplot10,
        nrow=2, ncol=5)

plotsBif <- sapply(cobraFiles, function(file){
  cobra <- readRDS(file)
  plotPerformanceCurve(cobra)
})
resPlots <- plotsBif[seq(1,length(plotsBif),by=9)] # get relevant data frame


pAll <- ggplot(resPlots[[1]],
               aes(x = FDP, y = TPR, col = method)) +
  geom_path(size = 1, aes(linetype = method)) +
  scale_color_manual(values = cols, breaks = names(cols)) +
  scale_linetype_manual(values = linetypes, breaks = names(linetypes))  + guides(col=guide_legend(nrow=3))

legend_all <- get_legend(pAll + labs(col = "", linetype = "") +
                           theme(legend.position = "bottom",
                                 legend.key.width = unit(1.3, "cm")))

pLeg <- plot_grid(pAll, legend_all, rel_heights=c(1,0.2), nrow=2, ncol=1)
pLeg

#### plot with trajectories
datasetClusters <- read.table("~/Dropbox/research/PhD/research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/datasetClustersSlingshot.txt", header=TRUE)

  FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
  }
pal <- wes_palette("Zissou1", 12, type = "continuous")

for(datasetIter in c(1:10)){

  data <- readRDS(paste0("~/Dropbox/research/PhD/research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/datasets/20190326_dyntoyDataset_", datasetIter, ".rds"))
  counts <- t(data$counts)

  # get milestones
  gid <- data$prior_information$groups_id
  gid <- gid[match(colnames(counts),gid$cell_id),]

  pal <- wes_palette("Zissou1", 12, type = "continuous")
  truePseudotime <- data$prior_information$timecourse_continuous
  g <- Hmisc::cut2(truePseudotime,g=12)

  # quantile normalization
  normCounts <- FQnorm(counts)

  ## dim red
  pca <- prcomp(log1p(t(normCounts)), scale. = FALSE)
  rd <- pca$x[,1:3]
  ## cluster
  nClusters <- datasetClusters$nClusters[datasetIter]
  set.seed(5)
  cl <- kmeans(rd, centers = nClusters)$cluster

  #lineages
  lin <- getLineages(rd, cl, start.clus=datasetClusters$start[datasetIter], end.clus=c(datasetClusters$end1[datasetIter], datasetClusters$end2[datasetIter]))
  #curves
  crv <- getCurves(lin)
  plot(rd, col = pal[g], pch=16, asp = 1)
  lines(crv, lwd=2, col="black")

  ggTraj <- ggplot(as.data.frame(rd), aes(x = PC1, y = PC2)) +
    geom_point(col = pal[g]) +
    theme_classic() +
    geom_path(x = crv@curves$curve1$s[crv@curves$curve1$ord, 1],
              y = crv@curves$curve1$s[crv@curves$curve1$ord, 2]) +
    geom_path(x = crv@curves$curve2$s[crv@curves$curve2$ord, 1],
              y = crv@curves$curve2$s[crv@curves$curve2$ord, 2]) +
    ggtitle(paste0("Dataset",datasetIter))
  assign(paste0("trajplot",datasetIter),ggTraj)
}

p1 <- plot_grid(trajplot1 + coord_fixed(), trajplot2 + coord_fixed(), trajplot3 + coord_fixed(), trajplot4 + coord_fixed(), trajplot5 + coord_fixed(),
              bifplot1 + coord_fixed(), bifplot2 + coord_fixed(), bifplot3 + coord_fixed(), bifplot4 + coord_fixed(), bifplot5 + coord_fixed(),
        nrow=2, ncol=5, rel_heights=c(0.8,1))
#pLeg1 <- plot_grid(p1, legend_all, rel_heights=c(1,0.15), nrow=2, ncol=1)
#pLeg1
p1
ggsave("~/Dropbox/research/PhD/research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/individualPerformanceFirst1To5.pdf", width = unit(15, "in"), height = unit(9, "in"), scale = .7)
dev.new()
p2 <- plot_grid(trajplot6 + coord_fixed(), trajplot7 + coord_fixed(), trajplot8 + coord_fixed(), trajplot9 + coord_fixed(), trajplot10 + coord_fixed(),
                bifplot6 + coord_fixed(), bifplot7 + coord_fixed(), bifplot8 + coord_fixed(), bifplot9 + coord_fixed(), bifplot10 + coord_fixed(),
        nrow=2, ncol=5, rel_heights=c(0.8,1,0.8,1))
pLeg2 <- plot_grid(p2, legend_all, rel_heights=c(1,0.15), nrow=2, ncol=1)
pLeg2
ggsave("~/Dropbox/research/PhD/research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/sim2_dyntoy_bifurcating_4/individualPerformance6To10.pdf", width = unit(15, "in"), height = unit(10, "in"), scale = .7)

### mean plot
# Monocle (hence also BEAM) only finds one lineage in datasets 3 6 9
# Monocle fails to separate the correct two lineages in datasets 4 7 8 10
# GPfates only finds one lineage in datasets 3 4 8 9
# iCOBRA does not seem to be consistent in the cut-offs, so difficult to merge different datasets => calculate FDP and TPR yourself.
resList <- c()
for(ii in 1:length(cobraFiles)){
    cobra <- readRDS(cobraFiles[ii])
    pvals <- pval(cobra)
    colnames(pvals) <- gsub(colnames(pvals),pattern="tradeR",replacement="tradeSeq")
    truths <- as.logical(truth(cobra)[,1])
    # performance for all p-value based methods
    hlp <- apply(pvals,2,function(x){
      pOrder <- order(x,decreasing=FALSE)
      padj <- p.adjust(x,"fdr")
      fdr <- cumsum(!truths[pOrder])/(1:length(pOrder))
      tpr <- cumsum(truths[pOrder])/sum(truths)
      df <- data.frame(fdr=fdr, tpr=tpr, cutoff=(1:length(padj))/length(padj))
    })
    # performance for GPfates
    scores <- score(cobra)
    pOrder <- order(scores[,1],decreasing=TRUE)
    fdr <- cumsum(!truths[pOrder])/(1:length(pOrder))
    tpr <- cumsum(truths[pOrder])/sum(truths)
    df <- data.frame(fdr=fdr, tpr=tpr, cutoff=(1:nrow(scores))/nrow(scores))
    hlp$GPfates <- df
    # summarize
    dfIter <- do.call(rbind,hlp)
    dfIter$method=rep(c(colnames(pvals),"GPfates"),each=nrow(pvals))
    dfIter$dataset <- ii
    resList[[ii]] <- dfIter
}

#### across all datasets
library(tidyverse)
df <- as_tibble(do.call(rbind,resList))
df <- df %>% group_by(method,cutoff) %>%
        summarize(meanTPR=mean(tpr,na.rm=TRUE),
                meanFDR=mean(fdr,na.rm=TRUE))
pMeanAll <- ggplot(df, aes(x=meanFDR, y=meanTPR, col=method)) + geom_path(size = 1, aes(linetype = method)) +
scale_x_continuous(limits = c(0, 1), breaks = c(0.01, 0.05, 0.1, 0.5, 1),
                   minor_breaks = c(0:5) * .1) +
scale_y_continuous(limits = c(0, 1)) +
scale_color_manual(values = cols, breaks = names(cols)) +
scale_linetype_manual(values = linetypes, breaks = names(linetypes)) + ggtitle("All 10 datasets") + guides(col=guide_legend(nrow=3)) + xlab("FDR") + ylab("TPR")


pMeanLegAll <- plot_grid(pMeanAll, legend_all, rel_heights=c(1,0.15), nrow=2, ncol=1)
pMeanLegAll


#### only for datasets where each method found the correct trajectory
resListCleaned <- resList
# remove BEAM for datasets where Monocle2 did not find a branching
resListCleaned[c(3,6,9)] <- lapply(resListCleaned[c(3,6,9)], function(df){
  df[!df$method=="BEAM",]
})
# remove Monocle_tradeSeq for datasets where Monocle
# (a) did not find a branching;
# (b) wrongly assigned the two lineages to the same branch.
resListCleaned[c(3,6,9,4,7,8,10)] <- lapply(resListCleaned[c(3,6,9,4,7,8,10)], function(df){
  df[-grep(df$method, pattern="tradeSeq_Monocle2*"),]
})
# remove GPfates for datasets where they only find a single lineage
resListCleaned[c(3,4,8,9)] <- lapply(resListCleaned[c(3,4,8,9)], function(df){
  df[!df$method=="GPfates",]
})
# remove GPfates_tradeSeq for those datasets too.
resListCleaned[c(3,4,8,9)] <- lapply(resListCleaned[c(3,4,8,9)], function(df){
  df[-grep(df$method, pattern="tradeSeq_GPfates*"),]
})

library(tidyverse)
df <- as_tibble(do.call(rbind,resListCleaned))
df <- df %>% group_by(method,cutoff) %>%
        summarize(meanTPR=mean(tpr,na.rm=TRUE),
                meanFDR=mean(fdr,na.rm=TRUE))


pMean <- ggplot(df, aes(x=meanFDR, y=meanTPR, col=method)) + geom_path(size = 1, aes(linetype = method)) +
scale_x_continuous(limits = c(0, 1), breaks = c(0.01, 0.05, 0.1),
                   minor_breaks = c(0:5) * .1) +
scale_y_continuous(limits = c(0, 1)) +
scale_color_manual(values = cols, breaks = names(cols)) +
scale_linetype_manual(values = linetypes, breaks = names(linetypes)) + ggtitle("Only datasets with correct trajectory per method.")

pMeanLeg <- plot_grid(pMean, legend_all, rel_heights=c(1,0.15), nrow=2, ncol=1)
pMeanLeg

### pretending we don't know the truth: use every dataset where >1 lineage is discovered for every method.
resListBif <- resList
# remove BEAM for datasets where Monocle2 did not find a branching
resListBif[c(3,6,9)] <- lapply(resListBif[c(3,6,9)], function(df){
  df[!df$method=="BEAM",]
})
# remove Monocle_tradeSeq for datasets where Monocle
# (a) did not find a branching;
resListBif[c(3,6,9)] <- lapply(resListBif[c(3,6,9)], function(df){
  df[-grep(df$method, pattern="tradeSeq_Monocle2*"),]
})
# remove GPfates for datasets where they only find a single lineage
resListBif[c(3,4,8,9)] <- lapply(resListBif[c(3,4,8,9)], function(df){
  df[!df$method=="GPfates",]
})
# remove GPfates_tradeSeq for those datasets too.
resListBif[c(3,4,8,9)] <- lapply(resListBif[c(3,4,8,9)], function(df){
  df[-grep(df$method, pattern="tradeSeq_GPfates*"),]
})

library(tidyverse)
df <- as_tibble(do.call(rbind,resListBif))
df <- df %>% group_by(method,cutoff) %>%
        summarize(meanTPR=mean(tpr,na.rm=TRUE),
                meanFDR=mean(fdr,na.rm=TRUE))


pMean <- ggplot(df, aes(x=meanFDR, y=meanTPR, col=method)) + geom_path(size = 1, aes(linetype = method)) +
scale_x_continuous(limits = c(0, 1), breaks = c(0.01, 0.05, 0.1),
                   minor_breaks = c(0:5) * .1) +
scale_y_continuous(limits = c(0, 1)) +
scale_color_manual(values = cols, breaks = names(cols)) +
scale_linetype_manual(values = linetypes, breaks = names(linetypes)) + ggtitle("Every bifurcating trajectory for each method.")

pMeanBif <- plot_grid(pMean, legend_all, rel_heights=c(1,0.15), nrow=2, ncol=1)
pMeanBif

# for intersection of datasets where each method found a correct trajectory; only datasets 1,2,5
resListIst <- resList
resListIst <- resListIst[c(1,2,5)]
library(tidyverse)
df <- as_tibble(do.call(rbind,resListIst))
df <- df %>% group_by(method,cutoff) %>%
        summarize(meanTPR=mean(tpr,na.rm=TRUE),
                meanFDR=mean(fdr,na.rm=TRUE))


pMeanIst <- ggplot(df, aes(x=meanFDR, y=meanTPR, col=method)) + geom_path(size = 1, aes(linetype = method)) +
scale_x_continuous(limits = c(0, 1), breaks = c(0.01, 0.05, 0.1),
                   minor_breaks = c(0:5) * .1) +
scale_y_continuous(limits = c(0, 1)) +
scale_color_manual(values = cols, breaks = names(cols)) +
scale_linetype_manual(values = linetypes, breaks = names(linetypes)) + xlab("FDR") + ylab("TPR")

pMeanIstLeg <- plot_grid(pMeanIst, legend_all, rel_heights=c(1,0.15), nrow=2, ncol=1)
pMeanIstLeg
saveRDS(pMeanIst,file="~/Dropbox/research/PhD/research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/performancePlots/pMeanIstBifurcation.rds")
