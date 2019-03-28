library(here)
library(slingshot)
library(RColorBrewer)
library(mgcv)
library(tradeR)
library(edgeR)
library(rafalib)
library(wesanderson)
library(ggplot2)
library(cowplot)
library(iCOBRA)
library(scales)
library(dplyr)


### prepare performance plots ####
cols <- c(rep(c("#C6DBEF", "#08306B"), each = 3), "#4292C6", "#4daf4a",
          "#e41a1c", "#e78ac3", "#ff7f00")
names(cols) <- c("tradeR_slingshot_end", "tradeR_GPfates_end", "tradeR_Monocle2_end",
                 "tradeR_slingshot_pattern", "tradeR_GPfates_pattern",
                 "tradeR_Monocle2_pattern", "tradeR_slingshot_assoc", "Monocle3_assoc",
                 "BEAM", "GPfates", "edgeR")
linetypes <- c(rep(c("dashed", "dotdash", "solid"), 2), rep("solid", 6))
names(linetypes) <- c("tradeR_slingshot_end", "tradeR_GPfates_end", "tradeR_Monocle2_end",
                      "tradeR_slingshot_pattern", "tradeR_GPfates_pattern",
                      "tradeR_Monocle2_pattern", "tradeR_slingshot_assoc", "Monocle3_assoc",
                      "BEAM", "GPfates", "edgeR")

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
dir <- "~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/datasets"
cobraFiles <- list.files(dir, pattern="cobra*", full.names=TRUE)


plotPerformanceCurve <- function(cobraObject){
  cobraObject <- calculate_adjp(cobraObject)
  cobraObject <- calculate_performance(cobraObject, binary_truth = "status")

  DyntoyPlot <- data.frame(FDP = cobraObject@fdrtprcurve$FDR,
                           TPR = cobraObject@fdrtprcurve$TPR,
                           method = cobraObject@fdrtprcurve$method,
                         cutoff = cobraObject@fdrtprcurve$CUTOFF)
  pDyntoy <- ggplot(DyntoyPlot, aes(x = FDP, y = TPR, col = method)) +
    geom_path(size = 1, aes(linetype = method)) +
    xlab("FDP") +
    scale_x_continuous(limits = c(0, 1), breaks = c(0.01, 0.05, 0.1),
                       minor_breaks = c(0:5) * .1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = cols, breaks = names(cols)) +
    scale_linetype_manual(values = linetypes, breaks = names(linetypes))
    pDyntoy
}


for(ii in 1:length(cobraFiles)){
    cobra <- readRDS(cobraFiles[ii])
    assign(paste0("bifplot",ii),plotPerformanceCurve(cobra))
}

p1 <- plot_grid(bifplot1, bifplot2, bifplot3, bifplot4,
          bifplot5, bifplot6, bifplot7, bifplot8,
        nrow=2, ncol=4)

plotsBif <- sapply(cobraFiles, function(file){
  cobra <- readRDS(file)
  plotPerformanceCurve(cobra)
})
resPlots <- plotsBif[seq(1,length(plotsBif),by=9)] # get relevant data frame


pAll <- ggplot(resPlots[[1]],
               aes(x = FDP, y = TPR, col = method)) +
  geom_path(size = 1, aes(linetype = method)) +
  scale_color_manual(values = cols, breaks = names(cols)) +
  scale_linetype_manual(values = linetypes, breaks = names(linetypes))

legend_all <- get_legend(pAll + labs(col = "", linetype = "") +
                           theme(legend.position = "bottom",
                                 legend.key.width = unit(1.3, "cm")))

pLeg <- plot_grid(p1, legend_all, rel_heights=c(1,0.2), nrow=2, ncol=1)

#### plot with trajectories
datasetClusters <- read.table("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/datasetClustersSlingshot.txt", header=TRUE)

  FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
  }
pal <- wes_palette("Zissou1", 12, type = "continuous")

for(datasetIter in c(1,2,5:10)){

  data <- readRDS(paste0("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_bifurcating_4/datasets/20190326_dyntoyDataset_", datasetIter, ".rds"))
  counts <- t(data$counts)
  truePseudotime <- data$prior_information$timecourse_continuous
  g <- Hmisc::cut2(truePseudotime,g=12)
  # quantile normalization
  normCounts <- FQnorm(counts)
  ## dim red
  pca <- prcomp(log1p(t(normCounts)), scale. = FALSE)
  rd <- pca$x[,1:4]
  ## cluster
  set.seed(9)
  nClusters <- 6
  cl <- kmeans(rd, centers = nClusters)$cluster
  #lineages
  lin <- getLineages(rd, cl, start.clus=datasetClusters$start[datasetIter], end.clus=c(datasetClusters$end1[datasetIter], datasetClusters$end2[datasetIter]))
  crv <- getCurves(lin)

  ggTraj <- ggplot(as.data.frame(rd), aes(x = PC1, y = PC2)) +
    geom_point(col = pal[g]) +
    theme_classic() +
    geom_path(x = crv@curves$curve1$s[crv@curves$curve1$ord, 1],
              y = crv@curves$curve1$s[crv@curves$curve1$ord, 2]) +
    geom_path(x = crv@curves$curve2$s[crv@curves$curve2$ord, 1],
              y = crv@curves$curve2$s[crv@curves$curve2$ord, 2])
  assign(paste0("trajplot",datasetIter),ggTraj)
}

p1 <- plot_grid(trajplot1, trajplot10, trajplot2, trajplot5,
              bifplot1, bifplot2, bifplot3, bifplot4,
              trajplot6, trajplot7, trajplot8, trajplot9,
          bifplot5, bifplot6, bifplot7, bifplot8,
        nrow=4, ncol=4, rel_heights=c(0.8,1,0.8,1))
pLeg <- plot_grid(p1, legend_all, rel_heights=c(1,0.15), nrow=2, ncol=1)
pLeg
