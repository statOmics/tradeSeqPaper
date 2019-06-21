library(plyr)
# download data from http://www.gs.washington.edu/~xqiu/proj2/RGE_analysis_data.tar.gz
data <- read.csv("~/Downloads/RGE_analysis_data/Nature_hta_paper/input/Olsson_RSEM_SingleCellRNASeq.csv", row.names=1)
library(stringr)
sample_sheet <- data.frame(groups = str_split_fixed(colnames(data), "\\.+", 3), row.names = colnames(data))
dataCleaned <- data[, sample_sheet$groups.1 %in% c('Lsk', 'Cmp', 'Gmp', 'LK')]
ssCleaned <- sample_sheet[sample_sheet$groups.1 %in% c('Lsk', 'Cmp', 'Gmp', 'LK'),]

## data from Monocle 2 paper
fig1b <- read.csv("~/Downloads/RGE_analysis_data/Nature_hta_paper/input/fig1b.txt",row.names=1, sep = '\t')
colClusters <- fig1b[1,-1]
rowClusters <- fig1b[-1,1]

# subset to Gmp
colClusters <- unlist(c(colClusters[substr(names(colClusters),1,3)=="Gmp"]))
sum(sample_sheet$groups.1=="Gmp")
names(colClusters) <- unlist(lapply(lapply(strsplit(names(colClusters),split=".", fixed=TRUE), function(x) x[-1]),paste,collapse="."))

# subset data to Gmp
dataGmp <- dataCleaned[,ssCleaned$groups.1=="Gmp"]
ssGmp <- ssCleaned[ssCleaned$groups.1=="Gmp",]








#match up the column name in fig1b to the colnames in data
#note that you should not run this mutliple times
fig1b_names <- colnames(fig1b)
match_id <- which(str_split_fixed(colnames(fig1b), "\\.+", 2)[, 2] %in% colnames(dataCleaned) == T)
fig1b_names[match_id] <- str_split_fixed(colnames(fig1b), "\\.+", 2)[match_id, 2]
no_match_id <- which(str_split_fixed(colnames(fig1b), "\\.+", 2)[, 2] %in% colnames(dataCleaned) == F)
fig1b_names[no_match_id] <- str_split_fixed(colnames(fig1b), "\\.\\.", 2)[no_match_id, 2]
colnames(fig1b)[2:383] <- fig1b_names[2:383]
clusters <- as.numeric(fig1b[1, 2:383])
cluster_assignments <- revalue(as.factor(clusters), c("1" = "HSCP-1", "2" = "HSCP-2", "3" = "Meg", "4" = "Eryth",
                                                                 "5" = "Multi-Lin", "6" = "MDP", "7" = "Mono", "8" = "Gran", "9" = "Myelocyte"))



## subset data for granulocyte and monocyte progenitor cell FACS sorted experiment
ss <- sample_sheet[sample_sheet$groups.1 %in% c('Lsk', 'Cmp', 'Gmp', 'LK'),]
data <- dataCleaned[,ss$groups.1=="Gmp"]
ssGmp <- ss[ss$groups.1=="Gmp",]
clustGmp <- droplevels(cluster_assignments[ss$groups.1=="Gmp"])

### PCA
FQnorm <- function(counts, type="median"){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  if(type=="mean"){
    refdist <- apply(counts.sort,1,mean)
  } else if(type=="median"){
    refdist <- apply(counts.sort,1,median)
  }
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
dataNorm <- FQnorm(dataGmp)
pc <- prcomp(log1p(t(dataNorm)))
plot(pc$x[,1:2], col=colClusters, pch=16)
