library(here)
library(rafalib)
library(wesanderson)
library(tidyverse)
library(lubridate)

palette(wes_palette("Darjeeling1", 10, type = "continuous"))

time_benchmarks <- list()
mem_benchmarks <- list()
for (size in 2:3) {
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
  mem_benchmark <- read.table(
    here::here("simulation", "time", paste0(size, "-mem-benchmark.txt")))
  mem_benchmark <- mem_benchmark %>%
    mutate(method = rownames(mem_benchmark),
           n = 10^size) %>%
    rename("mem" = x)
  mem_benchmarks[[size - 1]] <- mem_benchmark
}

time_benchmarks <- do.call("rbind", time_benchmarks) %>% as.data.frame()
# time_benchmarks <- time_benchmarks %>% add_row("expr" = "ImpulseDE",
#                                                "time" = duration(3.5, units = "hours"),
#                                                n = 100)
mem_benchmarks <- do.call("rbind", mem_benchmarks) %>% as.data.frame() %>%
#   add_row(mem = 148.2,
#           "method" = "ImpulseDE",
#           n = 100) %>%
  identity()

ggplot(time_benchmarks, aes(x = n, y = as.numeric(time) / 60,
                            group = interaction(expr, n), col = expr)) +
  geom_point() +
  theme_classic() +
  labs(x = "number of cells", y = "time (minutes)") +
  scale_x_log10()

ggplot(mem_benchmarks, aes(x = n, y = mem / 10 ^ 3,
                           group = interaction(method, n),
                           col = method)) +
  geom_point() +
  theme_classic() +
  labs(x = "number of cells", y = "memory (kB)") +
  scale_x_log10()
