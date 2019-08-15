library(here)
library(dyno)
library(dyntoy)
library(tidyverse)
library(wesanderson)

## Pre-process ----
NCORES <- 2
palette(wes_palette("Darjeeling1", 10, type = "continuous"))
source(here::here("simulation", "time", "20190611_helper.R"))
set.seed(87657)

## datasets ----

for (size in 2:4) {
  ## Generate dataset ----
  print("dataset")
  dataset <- generate_dataset(
    model = model_bifurcating(),
    num_cells = 10^size,
    num_features = 5000,
    differentially_expressed_rate = .2
  )
  saveRDS(dataset, here::here("simulation", "time", "data",
                              paste0(size, "dataset.rds")))
  
}