library(tidyverse)
library(dyno)
library(dyntoy)
library(patchwork)


  dataset <- generate_dataset(
    model = model_bifurcating(),
    num_cells = 10000,
    num_features = 5000,
    differentially_expressed_rate = .2
  )
  saveRDS(dataset, file=paste0("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeSeqPaper/simulation/time/bigDyntoyDataset.rds"))
