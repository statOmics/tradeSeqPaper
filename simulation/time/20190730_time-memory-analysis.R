library(here)
library(rafalib)
library(wesanderson)
library(profvis)
library(tidyverse)

palette(wes_palette("Darjeeling1", 10, type = "continuous"))

for (size in c("small", "big")) {
  ## Benchmark time ----
  time_benchmark <- read.table(
    file = here::here("simulation", "time", paste0(size, "-time-benchmark.txt")))
  time_benchmark <- time_benchmark %>%
    mutate(expr = case_when(str_detect(expr, "edgeR") ~ "edgeR",
                            str_detect(expr, "fitGAM") ~ "tradeSeq",
                            str_detect(expr, "BEAM") ~ "BEAM",
                            str_detect(expr, "ImpulseDE2") ~ "ImpulseDE2"),
           time = duration(round(time / 60^5,2), units = "seconds")) %>%
    group_by(expr) %>%
    summarise(min = duration(min(time), units = "seconds"),
              firstQuartile = quantile(time, 0.25),
              median = duration(median(time), units = "seconds"),
              mean = duration(mean(time), units = "seconds"),
              thirdQuartile = quantile(time, 0.75),
              max = duration(max(time), units = "seconds"))
  
  ## Benchmark memory ----
  ### tradeSeq
  print(profvis(prof_input = here::here("simulation", "time",
                            paste0(size, "-fitGam-memory.Rprof"))))
  ### For small dataset: -28559.2 / 28902.8 (MB)
  
  ### BEAM
  profvis(prof_input = here::here("simulation", "time",
                            paste0(size, "-BEAM-memory.Rprof")))
  ### For small dataset: -3689.8 / 3952.0 (MB)
  
  ### edgeR
  profvis(prof_input = here::here("simulation", "time",
                            paste0(size, "-edgeR-memory.Rprof")))
  ### For small dataset: -635.4	/ 785.7 (MB)
  
  ### ImpusleDE
  profvis(prof_input = here::here("simulation", "time",
                            paste0(size, "-ImpusleDE-memory.Rprof")))
}
