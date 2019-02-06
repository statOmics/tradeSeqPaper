library(tidyverse)
library(dyno)
library(dyntoy)
library(patchwork)

set.seed(1)
dataset <- generate_dataset(
  model = model_bifurcating(),
  num_cells = 500,
  num_features = 5000,
  differentially_expressed_rate = .2
)

g1 <- plot_graph(dataset)
g2 <- plot_heatmap(dataset %>% add_root(), features_oi = 100)

dimred <- dyndimred::dimred_landmark_mds(get_expression(dataset), num_landmarks = 1000)
g3 <- plot_dimred(dataset, dimred = dimred)

g4 <- plot_dimred(dataset, dimred = dyndimred::dimred_ica)

# you might have to press 'Zoom' to be able to view thos plot
patchwork::wrap_elements(g1 + g2 + g3 + g4) + ggtitle(glue::glue("ID = {dataset$id}; trajectory type = {dataset$trajectory_type}"))

dataset$tde_overall

saveRDS(dataset,file="~/Dropbox/PhD/Research/singleCell/trajectoryInference/trajectoryDE/tradeRPaper/simulation/sim2_dyntoy/20190206_dyntoyDataset.rds")
