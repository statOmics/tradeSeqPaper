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


FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method = 'min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}

### cyclic trajectory ####
data <- readRDS(here("simulation", "sim2_dyngen_cycle_72/72.rds"))
counts <- t(data$counts)
pal <- wes_palette("Zissou1", 12, type = "continuous")
truePseudotime <- data$prior_information$timecourse_continuous
g <- Hmisc::cut2(truePseudotime, g = 12)

# quantile normalization
normCounts <- FQnorm(counts)

# dim red
pca <- prcomp(log1p(t(normCounts)), scale. = FALSE)
rd <- pca$x[,1:2]
plot(rd, pch = 16, asp = 1, col = pal[g])

# pc
pcc <- principal_curve(rd, smoother = "periodic_lowess")
lines(x = pcc$s[order(pcc$lambda), 1], y = pcc$s[order(pcc$lambda), 2], col = "red", lwd = 2)
rdCycle <- rd
pccCycle <- pcc
gCycle <- g
rm(data, counts, pal, truePseudotime, g, normCounts, pca, rd, pcc)

### bifurcating dyntoy dataset ####
data <- readRDS(here("simulation", "sim2_dyntoy/20190206_dyntoyDataset.rds"))
counts <- t(data$counts)

pal <- wes_palette("Zissou1", 12, type = "continuous")
truePseudotime <- data$prior_information$timecourse_continuous
g <- Hmisc::cut2(truePseudotime,g = 12)

# quantile normalization
normCounts <- FQnorm(counts)

## dim red
pca <- prcomp(log1p(t(normCounts)), scale. = FALSE)
rd <- pca$x[,1:2]
plot(rd, pch = 16, asp = 1)

set.seed(9)
cl <- kmeans(rd, centers = 6)$cluster
plot(rd, col = brewer.pal(9,"Set1")[cl], pch = 16, asp = 1)
legend("topleft", legend = as.character(1:7), col = brewer.pal(9, "Set1")[1:7], pch = 16, cex = 2 / 3, bty = "n")

#lineages
lin <- getLineages(rd, cl, start.clus = 3, end.clus = c(1, 2))
plot(rd, col = pal[g], pch = 16, asp = 1)
lines(lin, lwd = 2)

#curves
crv <- getCurves(lin)
plot(rd, col = pal[g], pch = 16, asp = 1)
lines(crv, lwd = 2, col = "black")
rdDyntoy <- rd
crvDyntoy <- crv
gDyntoy <- g

rm(data, counts, pal, truePseudotime, g, normCounts, pca, rd, crv, lin, cl)

#### bifurcating dyngen dataset ####
data <- readRDS(here("simulation/sim2_dyntoy_bifurcating_4/bifurcating_4.rds"))
counts <- t(data$counts)

# this dataset has no null genes. We will simulate some.
set.seed(5)
null1 <- t(apply(counts,1,sample))
dimnames(null1) <- list(paste0("H",1:501),paste0("C",1:2011))
null2 <- t(apply(counts,1,sample))
dimnames(null2) <- list(paste0("H",502:1002),paste0("C",1:2011))
null3 <- t(apply(counts,1,sample))
dimnames(null3) <- list(paste0("H",1003:1503),paste0("C",1:2011))
counts <- rbind(counts,null1,null2,null3)
nullGenes <- rownames(counts)[substr(rownames(counts), 1, 1) == "H"]

pal <- wes_palette("Zissou1", 12, type = "continuous")
truePseudotime <- data$prior_information$timecourse_continuous[colnames(counts)]
g <- Hmisc::cut2(truePseudotime, g = 12)

# quantile normalization
normCounts <- FQnorm(counts)

## dim red
pca <- prcomp(log1p(t(normCounts)), scale. = FALSE)
rd <- pca$x[,1:3]
plot(rd, pch = 16, asp = 1)
set.seed(9)
cl <- kmeans(rd, centers = 8)$cluster
plot(rd, col = brewer.pal(9, "Set1")[cl], pch = 16, asp = 1)
legend("topleft", legend = as.character(1:7), col = brewer.pal(9, "Set1")[1:7],
       pch = 16, cex = 2 / 3, bty = "n")
plot(rd, col = pal[g], pch = 16, asp = 1)

#lineages
lin <- getLineages(rd, cl, start.clus = 6, end.clus = c(5, 2))
plot(rd, col = pal[g], pch = 16, asp = 1)
lines(lin, lwd = 2)

#curves
crv <- getCurves(lin)
plot(rd, col = pal[g], pch = 16, asp = 1)
lines(crv, lwd = 2, col = "black")
rdDyngen <- rd
crvDyngen <- crv
gDyngen <- g

rm(data, counts, truePseudotime, g, normCounts, pca, rd, crv, lin, cl,
   null1, null2, null3, nullGenes)

### prepare performance plots ####
cols <- c("#e41a1c", "#e78ac3", "#ff7f00", "#4daf4a",
          rep(c("#C6DBEF", "#08306B", "#4292C6"), each = 2))
names(cols) <- c("BEAM", "GPfates", "edgeR",  "Monocle3_assoc",
                 "tradeR_slingshot_end", "tradeR_GPfates_end",
                 "tradeR_slingshot_pattern", "tradeR_GPfates_pattern",
                 "tradeR_slingshot_assoc")
linetypes <- c(rep("solid", 5), "twodash", "solid", "twodash", "solid")
names(linetypes) <- c("BEAM", "GPfates", "edgeR",  "Monocle3_assoc",
                      "tradeR_slingshot_end", "tradeR_GPfates_end",
                      "tradeR_slingshot_pattern", "tradeR_GPfates_pattern",
                      "tradeR_slingshot_assoc")

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

### cyclic plot ####
cobraCycle <- readRDS(here("simulation", "sim2_dyngen_cycle_72", "cobraObject.rds"))
pval(cobraCycle) <- pval(cobraCycle)[,c("tradeR_slingshot_assoc", "Monocle3")]
colnames(pval(cobraCycle))[2] <- "Monocle3_assoc"
cobraCycle <- calculate_adjp(cobraCycle)
cobraCycle <- calculate_performance(cobraCycle, binary_truth = "status")

CyclePlot <- data.frame(FDP = cobraCycle@fdrtprcurve$FDR,
                         TPR = cobraCycle@fdrtprcurve$TPR,
                         method = cobraCycle@fdrtprcurve$method)
pCycle <- ggplot(CyclePlot, aes(x = FDP, y = TPR, col = method)) +
  geom_path(size = 1, aes(linetype = method)) +
  xlab("FDP") +
  scale_x_continuous(limits = c(0, 0.5), breaks = c(0.01, 0.05, 0.1),
                     minor_breaks = c(0:5) * .1) +
  scale_y_continuous(limits = c(.5, 1)) +
  scale_color_manual(values = cols, breaks = names(cols)) +
  scale_linetype_manual(values = linetypes, breaks = names(linetypes))

# cobraCyclePlot <- prepare_data_for_plot(cobraCycle,
#                                         colorscheme = cols[
#                                         match(sort(unique(cobraCycle@roc$method)),
#                                                                  names(cols))])
# facetted(cobraCyclePlot) <- FALSE
# pCycle <- plot_fdrtprcurve(cobraCyclePlot, pointsize = 0, xaxisrange = c(0, 0.5),
#                            yaxisrange = c(0.5, 1)) +
#                                         theme(legend.position = "none") + xlab("FDP") +
#                                         theme(axis.title.x = element_text(size = rel(1)),
#                                               axis.title.y = element_text(size = rel(1)),
#                                               axis.text.x = element_text(size = rel(.8)),
#                                               axis.text.y = element_text(size = rel(.8)))
rm(cobraCycle, CyclePlot)

### bifurcating dyntoy plot  ####
cobraDyntoy <- readRDS(here("simulation", "sim2_dyntoy", "cobra.rds"))
cobraDyntoy <- calculate_adjp(cobraDyntoy)
cobraDyntoy <- calculate_performance(cobraDyntoy, binary_truth = "status")

DyntoyPlot <- data.frame(FDP = cobraDyntoy@fdrtprcurve$FDR,
                         TPR = cobraDyntoy@fdrtprcurve$TPR,
                         method = cobraDyntoy@fdrtprcurve$method)
pDyntoy <- ggplot(DyntoyPlot, aes(x = FDP, y = TPR, col = method)) +
  geom_path(size = 1, aes(linetype = method)) +
  xlab("FDP") +
  scale_x_continuous(limits = c(0, 0.5), breaks = c(0.01, 0.05, 0.1),
                     minor_breaks = c(0:5) * .1) +
  scale_y_continuous(limits = c(.5, 1)) +
  scale_color_manual(values = cols, breaks = names(cols)) +
  scale_linetype_manual(values = linetypes, breaks = names(linetypes))

# cobraDyntoyPlot <- prepare_data_for_plot(cobraDyntoy,
#                                          colorscheme = cols[
#                                         match(sort(unique(cobraDyntoy@roc$method)),
#                                               names(cols))])
# facetted(cobraDyntoyPlot) <- FALSE
# pDyntoy <- plot_fdrtprcurve(cobraDyntoyPlot, pointsize = 0, xaxisrange = c(0, 0.5),
#                             yaxisrange = c(0.5, 1)) +
#                             theme(legend.position = "none") + xlab("FDP") +
#                             theme(axis.title.x = element_text(size = rel(1)),
#                                   axis.title.y = element_text(size = rel(1)),
#                                   axis.text.x = element_text(size = rel(.8)),
#                                   axis.text.y = element_text(size = rel(.8)))

rm(cobraDyntoy, DyntoyPlot)

### bifurcating dyngen ####
cobraDyntoy4 <- readRDS(here("simulation", "sim2_dyntoy_bifurcating_4", "cobra.rds"))
cobraDyntoy4 <- calculate_adjp(cobraDyntoy4)
cobraDyntoy4 <- calculate_performance(cobraDyntoy4, binary_truth = "status")

Dyntoy4Plot <- data.frame(FDP = cobraDyntoy4@fdrtprcurve$FDR,
                         TPR = cobraDyntoy4@fdrtprcurve$TPR,
                         method = cobraDyntoy4@fdrtprcurve$method)
pDyngen <- ggplot(Dyntoy4Plot, aes(x = FDP, y = TPR, col = method)) +
  geom_path(size = 1, aes(linetype = method)) +
  xlab("FDP") +
  scale_y_continuous(limits = c(.5, 1)) +
  scale_x_continuous(limits = c(0, 0.5), breaks = c(0.01, 0.05, 0.1),
                     minor_breaks = c(0:5) * .1) +
  scale_color_manual(values = cols, breaks = names(cols)) +
  scale_linetype_manual(values = linetypes, breaks = names(linetypes))

# cobraDyntoy4Plot <- prepare_data_for_plot(cobraDyntoy4,
#                                           colorscheme = cols[
#                                           match(sort(unique(cobraDyntoy4@roc$method)),
#                                                 names(cols))])
# facetted(cobraDyntoy4Plot) <- FALSE
# pDyngen <- plot_fdrtprcurve(cobraDyntoy4Plot, pointsize = 0,
#                             xaxisrange = c(0, 0.5),
#                             yaxisrange = c(0.5, 1)) +
#                                          theme(legend.position = "none") +
#                                          xlab("FDP") +
#                             theme(axis.title.x = element_text(size = rel(1)),
#                                   axis.title.y = element_text(size = rel(1)),
#                                   axis.text.x = element_text(size = rel(.8)),
#                                   axis.text.y = element_text(size = rel(.8)))

rm(cobraDyntoy4, Dyntoy4Plot)

## Get a common legend ####

cobraCycle <- readRDS(here("simulation", "sim2_dyngen_cycle_72", "cobraObject.rds"))
pval(cobraCycle) <- pval(cobraCycle)[,c("tradeR_slingshot_assoc", "Monocle3")]
colnames(pval(cobraCycle))[2] <- "Monocle3_assoc"
cobraCycle <- calculate_adjp(cobraCycle)
cobraCycle <- calculate_performance(cobraCycle, binary_truth = "status")
CyclePlot <- data.frame(FDP = cobraCycle@fdrtprcurve$FDR,
                             TPR = cobraCycle@fdrtprcurve$TPR,
                             method = cobraCycle@fdrtprcurve$method)

cobraDyntoy <- readRDS(here("simulation", "sim2_dyntoy", "cobra.rds"))
cobraDyntoy <- calculate_adjp(cobraDyntoy)
cobraDyntoy <- calculate_performance(cobraDyntoy, binary_truth = "status")

DyntoyPlot <- data.frame(FDP = cobraDyntoy@fdrtprcurve$FDR,
                         TPR = cobraDyntoy@fdrtprcurve$TPR,
                         method = cobraDyntoy@fdrtprcurve$method)

pAll <- ggplot(full_join(CyclePlot, DyntoyPlot),
               aes(x = FDP, y = TPR, col = method)) +
  geom_path(size = 1, aes(linetype = method)) +
  scale_color_manual(values = cols, breaks = names(cols)) +
  scale_linetype_manual(values = linetypes, breaks = names(linetypes))

legend_all <- get_legend(pAll + theme(legend.position = "bottom"))

### composite plot ####
ggCycle <- ggplot(as.data.frame(rdCycle), aes(x = PC1, y = PC2)) +
  geom_point(col = pal[gCycle]) +
  theme_classic() +
  geom_path(x = pccCycle$s[order(pccCycle$lambda), 1],
            y = pccCycle$s[order(pccCycle$lambda), 2]) +
  ggtitle("Cyclic dataset") +
  theme(plot.title = element_text(face = "bold", hjust = .5))

ggDyntoy <- ggplot(as.data.frame(rdDyntoy), aes(x = PC1, y = PC2)) +
  geom_point(col = pal[gDyntoy]) +
  theme_classic() +
  geom_path(x = crvDyntoy@curves$curve1$s[crvDyntoy@curves$curve1$ord, 1],
            y = crvDyntoy@curves$curve1$s[crvDyntoy@curves$curve1$ord, 2]) +
  geom_path(x = crvDyntoy@curves$curve2$s[crvDyntoy@curves$curve2$ord, 1],
            y = crvDyntoy@curves$curve2$s[crvDyntoy@curves$curve2$ord, 2]) +
  ggtitle("Dyntoy1 dataset") +
  theme(plot.title = element_text(face = "bold", hjust = .5))

ggDyngen <- ggplot(as.data.frame(rdDyngen), aes(x = PC1, y = PC2)) +
  geom_point(col = pal[gDyngen]) +
  theme_classic() +
  geom_path(x = crvDyngen@curves$curve1$s[crvDyngen@curves$curve1$ord, 1],
            y = crvDyngen@curves$curve1$s[crvDyngen@curves$curve1$ord, 2]) +
  geom_path(x = crvDyngen@curves$curve2$s[crvDyngen@curves$curve2$ord, 1],
            y = crvDyngen@curves$curve2$s[crvDyngen@curves$curve2$ord, 2]) +
  ggtitle("Dyntoy2 dataset") +
  theme(plot.title = element_text(face = "bold", hjust = .5))

# png("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/plots/wesandersonColorLegendHorizontal.png",
#     width = 2, height = 1, units = "in", res = 300)
# rafalib::mypar()
# plot(x = 0:12, y = seq(0, 1, length = 13),
#      type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
# for (i in 1:12) {
#   polygon(x = c(0, 1, 1, 0) + i - 1, y = c(0, 0, 1, 1), col = pal[i],
#           density = -1, border = NA)}
# mtext("True pseudotime", side = 3, at = 3, cex = 1.35)
# mtext("Min", side = 1, at = 0, cex = 1)
# mtext("Max", side = 1, at = 12, cex = 1)
# dev.off()

png("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/plots/wesandersonColorLegendVertical.png",
    width = 1.2, height = 2, units = "in", res = 300)
rafalib::mypar()
plot(y = 0:12, x = seq(0, 1, length = 13),
     type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
for (i in 1:12) {
  polygon(x = c(0, 1, 1, 0), y = c(0, 0, 1, 1) + i - 1, col = pal[i],
          density = -1, border = NA)
  }
mtext("True pseudotime", side = 3, at = .2, padj = -1.5, cex = .9)
mtext("Min", side = 1, at = 1, cex = .7, font = 2)
mtext("Max", side = 3, at = 1, cex = .7, font = 2)
dev.off()

p1 <- plot_grid(ggCycle,
                ggdraw() +
                  draw_image("wesandersonColorLegendVertical.png",
                             scale = 1),
                ggDyntoy, NULL,
                ggDyngen, pCycle, NULL, pDyntoy, NULL, pDyngen, 
                nrow = 2, ncol = 5, rel_heights = c(0.5, 1),
                rel_widths = c(1, 0.3, 1, 0.3, 1),
                labels = c("a", "", "c", "", "e", "b", "", "d", "", "f")
)
p1

p2 <- plot_grid(p1, legend_all, ncol = 1, rel_heights = c(1, .25))
p2
ggsave("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/plots/simPerformance.pdf")
