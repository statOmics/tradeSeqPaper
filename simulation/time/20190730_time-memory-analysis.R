library(here)
library(rafalib)
library(wesanderson)
library(tidyverse)
library(lubridate)

palette(wes_palette("Darjeeling1", 10, type = "continuous"))

time_benchmarks <- list()
mem_benchmarks <- list()
for (size in 2:5) {
  ## Benchmark time ----
  time_benchmark <- read.table(
    file = here::here("simulation", "time", paste0(size, "-time-benchmark.txt")))
  time_benchmark <- time_benchmark %>%
    mutate(expr = case_when(str_detect(expr, "edgeR") ~ "edgeR",
                            str_detect(expr, "fitGAM") ~ "tradeSeq",
                            str_detect(expr, "BEAM") ~ "BEAM",
                            str_detect(expr, "GPfates") ~ "GPfates"),
           time = duration(round(time / 60^5,2), units = "seconds")) %>%
    mutate(n = 10^size)
  time_benchmarks[[size - 1]] <- time_benchmark
  
  ## Benchmark memory ----
  mem_benchmarks[size - 1] <- read.table(here::here("simulation", "time",
                                         paste0(size, "-mem-benchmark.txt")))
}

time_benchmarks <- do.call("rbind", time_benchmarks) %>% as.data.frame()
mem_benchmarks <- do.call("rbind", mem_benchmarks) %>% as.data.frame()


