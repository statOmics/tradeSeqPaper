# you won't be able to run this unless you have access to the PRISM server

library(dynutils)
library(tidyverse)
library(qsub)

# load in design from dynbenchmark scripts/01-datasets/02a-simulate_dyngen_datasets.R

design2 <- bind_rows(design[rep(72, 10, ])
design2$seed <- sample(1:10000, nrow(design2))

qsub_config <- override_qsub_config(memory = "10G", max_wall_time = "24:00:00", num_cores = 1, name = "dyngen", wait = F, stop_on_error = FALSE, execute_before = c("source ~/activate_dynverse.sh no"), remove_tmp_folder = FALSE)

handle <- qsub_pmap(
  design2,
  simulate_dyngen,
  use_cache = TRUE,
  qsub_config = qsub_config
)

results <- qsub_retrieve(handle)

results %>% write_rds("~/datasets_for_koen.rds")
