libs <- c("cowplot", "tidyverse", "clusterExperiment", "RColorBrewer", "mgcv",
          "tradeR", "slingshot", "curl")
suppressMessages(suppressWarnings(
  sapply(libs, require, warn.conflicts = FALSE,
         character.only = TRUE, quietly = TRUE)
))
rm(libs)

load("~/PhD_Data/singleCell/fletcher/ShareWithKelly/E4c2b_slingshot_wsforkelly.RData")
library(slingshot)
library(rgl)
library(rafalib) ; mypar()
library(RColorBrewer)
library(mgcv)
library(tradeR)

origData <- t(pcax$x %*% t(pcax$rotation))
colpal=cc
rd <- X[,1:5]
lin <- getLineages(rd, clusterLabels = clus.labels, start.clus = "1", end.clus="4")
crv <- getCurves(lin)

#3D plot
rgl::plot3d(X[,1:3], t='p', col=colpal[as.numeric(clus.labels)],alpha=0.3, pch = 19, cex = 2, size=8, xlab="PC 1", ylab="PC 2", zlab="PC 3", aspect="iso", box=FALSE, axes=FALSE)
rgl::axes3d(tick=FALSE)
rgl::par3d(windowRect = c(20, 30, 800, 800))
for (i in seq_along(curves)){
  rgl::plot3d(crv@curves[[i]]$s[order(crv@curves[[i]]$lambda),1:3], type='l',add=TRUE, lwd=4,col=colpal[which.max(tail(lin@lineages[[i]],1)==levels(clus.labels))])
}




load("~/gamListOE_tradeR.rda")

path <- "https://www.cell.com/cms/10.1016/j.stem.2017.04.003/attachment/f2c195e6-972c-41e7-9e94-d4907b217a9e/mmc3"
curl_download(path, destfile = "~/OE_DE.xlsx")
OE_DE <- readxl::read_xlsx("~/OE_DE.xlsx")[,1:8]


## which psuedotime values to compare?
#### neuronal lineage
mypar(mfrow=c(1,2))
plot(slingPseudotime(crv)[,1],col=colpal[as.numeric(clus.labels)])
plot(1:13,pch=16,col=colpal[1:13])
ptHBC2 <- median(slingPseudotime(crv)[as.numeric(clus.labels)==5,1]) #delta HBC2 cells
ptGBC1 <- median(slingPseudotime(crv)[as.numeric(clus.labels)==3,1], na.rm=TRUE) #GBC cells, neuronal lineage

#### sustentacular lineage
plot(slingPseudotime(crv)[,3],col=colpal[as.numeric(clus.labels)])
plot(1:13,pch=16,col=colpal[1:13])
ptGBC3 <- median(slingPseudotime(crv)[as.numeric(clus.labels)==5,3], na.rm=TRUE) #GBC cells, sustentacular lineage
ptISUS <- median(slingPseudotime(crv)[as.numeric(clus.labels)==6,3], na.rm=TRUE) #iSUS cells

# neuronal lineage

OE_DE1 <- OE_DE %>% filter(str_detect(Contrast, "HBC2") &
                           !str_detect(Contrast, "HBC1"))
OE_DE1_GBC <- OE_DE %>% filter(Contrast=="GBC-∆HBC2")
OE_DE1_iSUS <- OE_DE %>% filter(Contrast=="∆HBC2-iSUS")

DE1 <- startVsEndTest(models = gamList, lineages=TRUE, pseudotimeValues=c(ptHBC2, ptGBC1))
DE1 <- DE1 %>% mutate(Feature = rownames(DE1)) %>%
               mutate(pvalue = lineage1) %>%
               select(Feature, pvalue) %>%
               mutate(p.adj = p.adjust(pvalue, "fdr")) %>%
               filter(p.adj < 0.05)
# 2063 genes
mean(OE_DE1_GBC$Feature %in% DE1$Feature) #25%
mean(DE1$Feature %in% OE_DE1_GBC$Feature) #40%


Common1 <- inner_join(DE1, OE_DE1_GBC)
Not_consistent <- Common1 %>% group_by(Feature) %>%
                              filter(n() > 1) %>%
                              select(logFC) %>%
                              mutate(rank = row_number()) %>%
                              spread(key = rank, value = logFC) %>%
                              filter(`1` * `2` < 0)
Common1 %>% select(Feature) %>% distinct() %>% nrow()
# 916 genes are common

Not_common <- anti_join(DE1, OE_DE1_GBC)
Other_contrasts <- inner_join(Not_common, OE_DE) %>%
                                  filter(str_detect(Contrast, "GBC") |
                                         str_detect(Contrast, "iSUS") |
                                         str_detect(Contrast, "2-.HBC1$"))
Other_contrasts %>% select(Feature) %>% distinct() %>% nrow()
# 268 genes are in other contrasts
round(100 * 24 / 238, 2)

Not_Detected <- anti_join(Not_common, OE_DE1_GBC)
Not_Detected %>% select(Feature) %>% distinct() %>% nrow()
# 1246 genes are not detected
round(100 * 59 / 238, 2)

Weird_Detected <- inner_join(Not_common, OE_DE1_GBC) %>%
                  anti_join(Other_contrasts, by = c("Feature" = "Feature"))
Weird_Detected %>% select(Feature) %>% distinct() %>% nrow()
# 0 genes are detected in contrasts further away
round(100 * 25 / 238, 2)
