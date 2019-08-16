load("~/PhD_Data/singleCell/fletcher/ShareWithKelly/E4c2b_slingshot_wsforkelly.RData")
library(slingshot)
library(rgl)
library(rafalib)
mypar()
library(RColorBrewer)
library(mgcv)
library(tradeSeq)

origData <- t(pcax$x %*% t(pcax$rotation))
colpal <- cc
rd <- X[, 1:5]
lin <- getLineages(rd, clusterLabels = clus.labels, start.clus = "1", end.clus = "4")
crv <- getCurves(lin)
plot(X[, 1:2], col = colpal[as.factor(clus.labels)], pch = 16)
lines(crv)


# knots are at pseudotime values of  0.00000  99.36628 170.88884
knots <-  c(0,  99.37, 170.89, 243.97, 375.2, 421.79)
pseudotime <- slingPseudotime(crv, na=FALSE)

# function that finds a cell closest to knot point
findCell <- function(knot, pseudotime){
  which.min(abs(pseudotime-knot))
}
# lineage 1
closestCell1 <- c()
for(ii in 1:length(knots)) closestCell1[ii] <- findCell(knots[ii], pseudotime[,1])
pseudotime[closestCell1,1] #OK
# lineage 2
closestCell2 <- c()
for(ii in 1:5) closestCell2[ii] <- findCell(knots[ii], pseudotime[,2])
pseudotime[closestCell2,2] #OK
# lineage 3
closestCell3 <- c()
for(ii in 1:4) closestCell3[ii] <- findCell(knots[ii], pseudotime[,3])
pseudotime[closestCell3,3] #OK

cols <- c("#FF7F00", "#1F78B4", "#E7298A")

colpal <- c("#1B9E77","antiquewhite2","cyan","#E7298A","#A6CEE3","#666666","#E6AB02","#FFED6F","darkorchid2","#B3DE69","#FF7F00","#A6761D","#1F78B4")

# 3D plot
rgl::plot3d(X[, 1:3], t = "p", col = colpal[as.factor(clus.labels)], alpha = 0.3, pch = 19, cex = 2, size = 8, xlab = "PC 1", ylab = "PC 2", zlab = "PC 3", aspect = "iso", box = FALSE, axes = FALSE)
# rgl::plot3d(X[,1:3], t='p', col=c("white","black")[hlp+1],alpha=0.3, pch = 19, cex = 2, size=8, xlab="PC 1", ylab="PC 2", zlab="PC 3", aspect="iso", box=FALSE, axes=FALSE)
rgl::axes3d(tick = FALSE)
rgl::par3d(windowRect = c(20, 30, 800, 800))
for (i in seq_along(curves)) {
  rgl::plot3d(crv@curves[[i]]$s[order(crv@curves[[i]]$lambda), 1:3], type = "l", add = TRUE, lwd = 4, col = colpal[which.max(tail(lin@lineages[[i]], 1) == levels(clus.labels))])
}
#rgl::pch3d(crv@curves[[1]]$s[closestCell,1:3], pch=2, add=TRUE)
rgl::pch3d(crv@curves[[1]]$s[closestCell,1:3], pch=as.character(1:6), add=TRUE, cex=3/2)
#rgl::pch3d(crv@curves[[2]]$s[closestCell2[-1],1:3], pch=2, add=TRUE)
rgl::pch3d(crv@curves[[2]]$s[closestCell2[-1],1:3], pch=as.character(2:5), add=TRUE, cex=3/2)
#rgl::pch3d(crv@curves[[3]]$s[closestCell3[-1],1:3], pch=2, add=TRUE, size=9)
rgl::pch3d(crv@curves[[3]]$s[closestCell3[-1],1:3], pch=as.character(2:4), add=TRUE, cex=3/2)
rgl.postscript("~/fletcher3d_knots.pdf", fmt="pdf")
