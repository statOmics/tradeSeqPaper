library(here)
library(rafalib)
library(wesanderson)
library(profvis)
library(tidyverse)
library(lubridate)

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
  summaryRprof(filename = here::here("simulation", "time",
                                     paste0(size, "-fitGam-memory.Rprof")),
               memory = "both", lines = "show")$by.line[, "mem.total"] %>%
    print()
  
  ### BEAM
  summaryRprof(filename = here::here("simulation", "time",
                                     paste0(size, "-BEAM-memory.Rprof")),
               memory = "both", lines = "show")$by.line[, "mem.total"] %>%
    print()
  
  ### edgeR
  summaryRprof(filename = here::here("simulation", "time",
                                     paste0(size, "-edgeR-memory.Rprof")),
               memory = "both", lines = "show")$by.line[, "mem.total"] %>%
    print()
  
  ### ImpusleDE
  summaryRprof(filename = here::here("simulation", "time",
                                     paste0(size, "-ImpulseDE-memory.Rprof")),
               memory = "both", lines = "show", chunksize = 10000) %>%
    print()
}
