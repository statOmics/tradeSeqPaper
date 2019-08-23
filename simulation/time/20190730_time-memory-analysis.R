library(here)
library(rafalib)
library(wesanderson)
library(tidyverse)
library(lubridate)
library(cowplot)

palette(wes_palette("Darjeeling1", 10, type = "continuous"))

time_benchmarks <- list()
mem_benchmarks <- list()
for (size in 2:4) {
  ## Benchmark time ----
  time_benchmark <- read.table(
    file = here::here("simulation", "time", "data",
                      paste0(size, "-time-benchmark.txt")))
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
    here::here("simulation", "time", "data",
               paste0(size, "-mem-benchmark.txt")))
  mem_benchmark <- mem_benchmark %>%
    mutate(method = rownames(mem_benchmark),
           n = 10^size) %>%
    rename("mem" = x)
  mem_benchmarks[[size - 1]] <- mem_benchmark
}

# Plots ----
time_benchmarks <- do.call("rbind", time_benchmarks) %>% as.data.frame()
time_benchmarks <- time_benchmarks %>% add_row("expr" = "ImpulseDE2",
                                               "time" = duration(3.5, units = "hours"),
                                               n = 100) %>%
  group_by(expr, n) %>%
  summarise(time = mean(time))
mem_benchmarks <- do.call("rbind", mem_benchmarks) %>% as.data.frame() %>%
  add_row(mem = 148.2,
          "method" = "ImpulseDE2",
          n = 100) %>%
  group_by(method, n) %>%
  summarise(mem = mean(mem))

cols <- c("#4292C6", "#e41a1c", "#e78ac3", "#ff7f00", "darkgoldenrod1")
names(cols) <- c("tradeSeq", "BEAM", "GPfates", "edgeR", "ImpulseDE2")
p1 <- ggplot(time_benchmarks, aes(x = n, y = as.numeric(time) / 60,
                            group = expr, col = expr)) +
  geom_point(size = 5) +
  geom_line() +
  theme_classic() +
  labs(x = "number of cells", y = "time (minutes)", col = 'method') +
  scale_x_log10() +
  guides(col = FALSE) +
  scale_color_manual(values = cols, breaks = names(cols))
p2 <- ggplot(mem_benchmarks, aes(x = n, y = mem / 10 ^ 3,
                           group = method, n,
                           col = method)) +
  geom_point(size = 5) +
  theme_classic() +
  geom_line() +
  labs(x = "number of cells", y = "memory (MB)") +
  scale_x_log10() +
  scale_color_manual(values = cols, breaks = names(cols))
p <- plot_grid(p1, p2, rel_widths = c(0.45, 0.6), rel_heights = c(0.8, 0.8))
ggsave(filename = here::here('simulation', 'time', 'figures', 'benchmark.pdf'),
       plot = p, width = 8, height = 4)

# Test benchmarks ----
test_benchmarks <- list()
for (size in 2:4) {
  ## Benchmark time ----
  test_benchmark <- read.table(
    file = here::here("simulation", "time", "data",
                      paste0(size, "-time-tests.txt")))
  test_benchmark <- test_benchmark %>%
    mutate(expr = case_when(str_detect(expr, "diffEndTest") ~ "diffEndTest",
                            str_detect(expr, "patternTest") ~ "patternTest",
                            str_detect(expr, "earlyDETest") ~ "earlyDETest",
                            str_detect(expr, "associationTest") ~ "associationTest",
                            str_detect(expr, "startVsEndTest") ~ "startVsEndTest"),
           time = duration(round(time / 60^5,2), units = "seconds")) %>%
    mutate(n = 10^size)
  test_benchmarks[[size - 1]] <- test_benchmark
}
test_benchmarks <- do.call("rbind", test_benchmarks) %>% as.data.frame()
ggplot(test_benchmarks, aes(x = n, y = as.numeric(time) / 60,
                            group = interaction(expr, n), col = expr)) +
  geom_boxplot() +
  # geom_line() +
  theme_classic() +
  labs(x = "number of cells", y = "time (minutes)", col = 'Tests') +
  scale_x_log10()
ggsave(here::here('simulation', 'time', 'figures', 'test.pdf'))
