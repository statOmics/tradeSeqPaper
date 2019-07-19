## we will analyze the bulk RNA-seq time course study from Kiselov et al. https://academic.oup.com/nar/article/43/20/9663/1395034

download.file("https://raw.githubusercontent.com/daniel-spies/rna-seq_tcComp/master/Biological_data_Kiselov_et%20al.txt", destfile="~/kiselov.txt")
rl <- readLines("~/kiselov.txt")
headerIds <- grep(pattern="30M_", x=rl, fixed=TRUE)
rl[headerIds]

# Comparison between PIK3CA mutant stimulated with EGF (30M_3rep_6TP_GSE69822_pik3ca_treat) versus wild type stimulated with EGF (30M_3rep_6TP_GSE69822_pik3ca_contr)

# wild type data
hlpWT <- rl[(headerIds[1]+2):(headerIds[2]-3)]
listWT <- lapply(hlpWT,function(x) strsplit(x, split=" "))
# get gene names
geneNames <- lapply(listWT, function(x) x[[1]][1])
# drop gene names from list, and convert to matrix
matWT <- do.call(rbind,lapply(listWT, function(x) as.numeric(x[[1]][-1])))
colnames(matWT) <- strsplit(rl[(headerIds[1]+1)],split=" ")[[1]]
rownames(matWT) <- geneNames
head(matWT)

# mutant data: PI3KA
hlpTreat <- rl[(headerIds[4]+2):(headerIds[5]-3)]
listTreat <- lapply(hlpTreat,function(x) strsplit(x, split=" "))
# get gene names
geneNames <- lapply(listTreat, function(x) x[[1]][1])
# drop gene names from list, and convert to matrix
matTreat <- do.call(rbind,lapply(listTreat, function(x) as.numeric(x[[1]][-1])))
colnames(matTreat) <- strsplit(rl[(headerIds[4]+1)],split=" ")[[1]]
rownames(matTreat) <- geneNames
head(matTreat)

# make count matrix
counts <- cbind(matWT, matTreat)
time <- as.numeric(unlist(lapply(strsplit(colnames(counts),split="_"),"[",2)))
treat <- factor(unlist(lapply(strsplit(colnames(counts),split="_"),"[",1)), levels=c("WT", "PIK3CA"))

# filter counts: at least 5 counts in 3 replicats
keep <- rowSums(counts>5)>=3
counts <- counts[keep,]

## tradeSeq
library(tradeSeq)
time <- matrix(time, nrow=ncol(counts), ncol=2, byrow=FALSE)
rownames(time) <- colnames(counts)
weights <- matrix(0, nrow=ncol(counts), ncol=2)
weights[1:(ncol(counts)/2), 1] <- 1
weights[(ncol(counts)/2+1):ncol(counts), 2] <- 1

## evaluate optimal K
infMat <- evaluateK(counts, pseudotime=time, cellWeights=weights, nGenes=250, k=3:6)

## fit GAM
gamList <- fitGAM(counts, pseudotime=time, cellWeights=weights, nknots=6)

## look at fit for six random genes
set.seed(81)
id <- sample(1:length(gamList), size=6)
par(mfrow=c(2,3))
for(ii in 1:6) plotSmoothers(gamList[[id[ii]]], main=resPat[id[ii],"pvalue"], ylim=c(3,9))


## test for different expression pattern
resPat <- patternTest(gamList, nPoints=8)

## statistical significance
hist(resPat$pvalue)
sum(p.adjust(resPat$pvalue, "fdr") <= 0.01, na.rm=TRUE)
# in the paper, they did pairwise comparisons with DESeq2 and found 7486 DE genes. We have 7184, which is very comparable.
deTradeSeq <- rownames(counts)[p.adjust(resPat$pvalue, "fdr") <= 0.01]

### compare with list from paper
download.file("https://github.com/wikiselev/rnaseq.mcf10a/blob/master/data/diff.expr.all.rda?raw=true", destfile="~/Downloads/diff.expr.all.rda")
load("~/Downloads/diff.expr.all.rda")
sum(diff.expr.pten.wt$padj < 0.01, na.rm=TRUE) #7825: PTEN-KO
sum(diff.expr.ki.wt$padj < 0.01, na.rm=TRUE)  #7486
deOriginal <- rownames(diff.expr.ki.wt)[which(diff.expr.ki.wt$padj <= 0.01)]


## overlap
mean(deTradeSeq %in% deOriginal)


## plot six most significant genes.
oo <- order(resPat$pvalue, decreasing=FALSE)
par(mfrow=c(2,3))
for(ii in seq(1,6)) plotSmoothers(gamList[[oo[ii]]], main=rownames(resPat)[oo[ii]], ylim=c(1,9))
dev.new()
par(mfrow=c(2,3))
for(ii in seq(7,12)) plotSmoothers(gamList[[oo[ii]]], main=rownames(resPat)[oo[ii]], ylim=c(1,9))

## cluster
deGenes <- rownames(counts)[p.adjust(resPat$pvalue, "fdr") <= 0.01]
deGenes <- deGenes[deGenes %in% names(gamList)]
clusRes <- clusterExpressionPatterns(gamList, nPoints=8, genes=deGenes, nReducedDims=4)
