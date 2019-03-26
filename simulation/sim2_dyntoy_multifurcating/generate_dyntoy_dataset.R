library(tidyverse)
library(dyno)
library(dyntoy)
library(patchwork)

set.seed(1)
for (i in 1:10) {
  dataset <- generate_dataset(
    model = model_multifurcating(),
    num_cells = 500,
    num_features = 5000,
    differentially_expressed_rate = .2
  )
}

saveRDS(dataset, file=paste0("~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy_multifurcating/20190326_dyntoyDataset_", i, ".rds"))
