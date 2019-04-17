setwd("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/")
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


FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method = 'min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}

### cyclic trajectory: take first dataset ####
dataAll <- readRDS(here::here("simulation", "sim2_dyngen_cycle_72/datasets/datasets_for_koen.rds"))
data <- dataAll[[1]]
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
rm(data, counts, truePseudotime, g, normCounts, pca, rd, pcc)

### bifurcating dyntoy dataset ####
data <- readRDS(here("simulation", "sim2_dyntoy_bifurcating_4/datasets/20190326_dyntoyDataset_1.rds"))
counts <- t(data$counts)

pal <- wes_palette("Zissou1", 12, type = "continuous")
truePseudotime <- data$prior_information$timecourse_continuous
g <- Hmisc::cut2(truePseudotime,g = 12)

# quantile normalization
normCounts <- FQnorm(counts)

## dim red
pca <- prcomp(log1p(t(normCounts)), scale. = FALSE)
rd <- pca$x[,1:2]
plot(rd, pch = 16, asp = 1, col=pal[g])

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
rdBif <- rd
crvBif <- crv
gBif <- g

rm(data, counts, truePseudotime, g, normCounts, pca, rd, crv, lin, cl)

#### multifurcating dataset ####
data <- readRDS(here("simulation","/sim2_dyntoy_multifurcating_4/data.rds"))
counts <- t(data$counts)

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
lin <- getLineages(rd, cl, start.clus = 4, end.clus = c(1, 2, 8))
plot(rd, col = pal[g], pch = 16, asp = 1)
lines(lin, lwd = 2)

#curves
crv <- getCurves(lin)
plot(rd, col = pal[g], pch = 16, asp = 1)
lines(crv, lwd = 2, col = "black")
rdMulti <- rd
crvMulti <- crv
gMulti <- g

rm(data, counts, truePseudotime, g, normCounts, pca, rd, crv, lin, cl)

### prepare performance plots ####
cols <- c(rep(c("#C6DBEF", "#08306B"), each = 3), "#4292C6", "#4daf4a",
          "#e41a1c", "#e78ac3", "#ff7f00")
names(cols) <- c("tradeSeq_slingshot_end", "tradeSeq_GPfates_end", "tradeSeq_Monocle2_end",
                 "tradeSeq_slingshot_pattern", "tradeSeq_GPfates_pattern",
                 "tradeSeq_Monocle2_pattern", "tradeSeq_slingshot_assoc", "Monocle3_assoc",
                 "BEAM", "GPfates", "edgeR")
linetypes <- c(rep(c("dashed", "dotdash", "solid"), 2), rep("solid", 6))
names(linetypes) <- c("tradeSeq_slingshot_end", "tradeSeq_GPfates_end", "tradeSeq_Monocle2_end",
                      "tradeSeq_slingshot_pattern", "tradeSeq_GPfates_pattern",
                      "tradeSeq_Monocle2_pattern", "tradeSeq_slingshot_assoc", "Monocle3_assoc",
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

### cyclic performance ####
cyclePlot <- readRDS(here("simulation", "performancePlots", "pMeanCycle.rds"))

pCycle <- cyclePlot + scale_x_continuous(limits = c(0, 0.6), breaks = c(0.01, 0.05, 0.1,.5,1),
                   minor_breaks = c(0:5) * .1) +
scale_y_continuous(limits = c(.5, 1)) +
scale_color_manual(values = cols, breaks = names(cols)) +
scale_linetype_manual(values = linetypes, breaks = names(linetypes))

rm(cyclePlot)

### bifurcating performance  ####
bifPlot <- readRDS(here("simulation", "performancePlots", "pMeanIstBifurcation.rds"))

pDyntoyBif <- bifPlot +
  scale_x_continuous(limits = c(0, 0.5), breaks = c(0.01, 0.05, 0.1, .5, 1),
                     minor_breaks = c(0:5) * .1) +
  scale_y_continuous(limits = c(.5, 1)) +
  scale_color_manual(values = cols, breaks = names(cols)) +
  scale_linetype_manual(values = linetypes, breaks = names(linetypes))

rm(bifPlot)

### multifurcating performance ####
cobraDyntoy4 <- readRDS(here("simulation", "sim2_dyntoy_multifurcating_4", "cobra.rds"))
cobraDyntoy4 <- calculate_adjp(cobraDyntoy4)
cobraDyntoy4 <- calculate_performance(cobraDyntoy4, binary_truth = "status")

multiPlot <- data.frame(FDP = cobraDyntoy4@fdrtprcurve$FDR,
                         TPR = cobraDyntoy4@fdrtprcurve$TPR,
                         method = cobraDyntoy4@fdrtprcurve$method)
pMulti <- ggplot(multiPlot, aes(x = FDP, y = TPR, col = method)) +
  geom_path(size = 1, aes(linetype = method)) +
  xlab("FDP") +
  scale_y_continuous(limits = c(.5, 1)) +
  scale_x_continuous(limits = c(0, 0.5), breaks = c(0.01, 0.05, 0.1, .5, 1),
                     minor_breaks = c(0:5) * .1) +
  scale_color_manual(values = cols, breaks = names(cols)) +
  scale_linetype_manual(values = linetypes, breaks = names(linetypes))

rm(multiPlot)

## Get a common legend ####
cobraCycle <- readRDS(here("simulation", "sim2_dyngen_cycle_72", "datasets", "cobra1.rds"))
pval(cobraCycle) <- pval(cobraCycle)[,c("tradeSeq_slingshot_assoc", "Monocle3")]
colnames(pval(cobraCycle))[2] <- "Monocle3_assoc"
cobraCycle <- calculate_adjp(cobraCycle)
cobraCycle <- calculate_performance(cobraCycle, binary_truth = "status")
CyclePlot <- data.frame(FDP = cobraCycle@fdrtprcurve$FDR,
                             TPR = cobraCycle@fdrtprcurve$TPR,
                             method = cobraCycle@fdrtprcurve$method)

cobraBif <- readRDS(here("simulation", "sim2_dyntoy_bifurcating_4", "datasets", "cobra1.rds"))
colnames(pval(cobraBif)) = gsub(colnames(pval(cobraBif)), pattern="tradeR",replacement="tradeSeq")
cobraBif <- calculate_adjp(cobraBif)
cobraBif <- calculate_performance(cobraBif, binary_truth = "status")

bifPlot <- data.frame(FDP = cobraBif@fdrtprcurve$FDR,
                         TPR = cobraBif@fdrtprcurve$TPR,
                         method = cobraBif@fdrtprcurve$method)

pAll <- ggplot(full_join(CyclePlot, bifPlot),
               aes(x = FDP, y = TPR, col = method)) +
  geom_path(size = 1, aes(linetype = method)) +
  scale_color_manual(values = cols, breaks = names(cols)) +
  scale_linetype_manual(values = linetypes, breaks = names(linetypes))

legend_all <- get_legend(pAll + labs(col = "", linetype = "") +
                           theme(legend.position = "bottom",
                                 legend.key.width = unit(1.3, "cm")) +
                               guides(col=guide_legend(nrow=3)))

### composite plot ####
ggCycle <- ggplot(as.data.frame(rdCycle), aes(x = PC1, y = PC2)) +
  geom_point(col = pal[gCycle]) +
  theme_classic() +
  geom_path(x = pccCycle$s[order(pccCycle$lambda), 1],
            y = pccCycle$s[order(pccCycle$lambda), 2]) +
  ggtitle("Cyclic dataset") +
  theme(plot.title = element_text(face = "bold", hjust = .5))

ggBif <- ggplot(as.data.frame(rdBif), aes(x = PC1, y = PC2)) +
  geom_point(col = pal[gBif]) +
  theme_classic() +
  geom_path(x = crvBif@curves$curve1$s[crvBif@curves$curve1$ord, 1],
            y = crvBif@curves$curve1$s[crvBif@curves$curve1$ord, 2]) +
  geom_path(x = crvBif@curves$curve2$s[crvBif@curves$curve2$ord, 1],
            y = crvBif@curves$curve2$s[crvBif@curves$curve2$ord, 2]) +
  ggtitle("Bifurcating dataset") +
  theme(plot.title = element_text(face = "bold", hjust = .5))

ggMulti <- ggplot(as.data.frame(rdMulti), aes(x = PC1, y = PC2)) +
  geom_point(col = pal[gMulti]) +
  theme_classic() +
  geom_path(x = crvMulti@curves$curve1$s[crvMulti@curves$curve1$ord, 1],
            y = crvMulti@curves$curve1$s[crvMulti@curves$curve1$ord, 2]) +
  geom_path(x = crvMulti@curves$curve2$s[crvMulti@curves$curve2$ord, 1],
            y = crvMulti@curves$curve2$s[crvMulti@curves$curve2$ord, 2]) +
  geom_path(x = crvMulti@curves$curve3$s[crvMulti@curves$curve3$ord, 1],
            y = crvMulti@curves$curve3$s[crvMulti@curves$curve3$ord, 2]) +
  ggtitle("Multifurcating dataset") +
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
                  draw_image("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/plots/wesandersonColorLegendVertical.png",
                             scale = 1),
                ggBif, NULL,
                ggMulti, pCycle, NULL, pDyntoyBif, NULL, pMulti,
                nrow = 2, ncol = 5, rel_heights = c(0.5, 1),
                rel_widths = c(1, 0.3, 1, 0.3, 1),
                labels = c("a", "", "b", "", "c", "d", "", "e", "", "f")
)
p1

p2 <- plot_grid(p1, legend_all, ncol = 1, rel_heights = c(1, .25))
p2

# ggsave("~/Documents/simPerformance.pdf", width = unit(15, "in"), height = unit(10, "in"), scale = .7)

ggsave("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/performancePlots/simPerformance.pdf", width = unit(15, "in"), height = unit(10, "in"), scale = .7)
