libs <- c("cowplot", "tidyverse", "clusterExperiment", "RColorBrewer", "mgcv", 
          "tradeR", "slingshot", "curl")
suppressMessages(suppressWarnings(
  sapply(libs, require, warn.conflicts = FALSE,
         character.only = TRUE, quietly = TRUE)
))
rm(libs)

setwd("../../OE")
load("data/gamListOE_tradeR.rda")

path <- "https://www.cell.com/cms/10.1016/j.stem.2017.04.003/attachment/f2c195e6-972c-41e7-9e94-d4907b217a9e/mmc3"
curl_download(path, destfile = "clustering/data/OE_DE.xlsx")
OE_DE <- readxl::read_xlsx("data/OE_DE.xlsx")[,1:8]

# 1st branching ----

OE_DE1 <- OE_DE %>% filter(str_detect(Contrast, "HBC2") &
                           !str_detect(Contrast, "HBC1"))

DE1 <- earlyDETest(models = gamList, knots = c(2, 4),
                  nPoints = 50, omnibus = T, pairwise = T)
DE1 <- DE1 %>% mutate(Feature = rownames(DE1)) %>%
               mutate(pvalue = 1 - pchisq(waldStat_1vs3 + waldStat_2vs3,
                                          df = df_1vs3 + df_2vs3)) %>%
               select(Feature, pvalue) %>%
               mutate(p.adj = p.adjust(pvalue)) %>%
               filter(p.adj < 0.05)
# 238 genes

Common1 <- inner_join(DE1, OE_DE1)
Not_consistent <- Common1 %>% group_by(Feature) %>%
                              filter(n() > 1) %>%
                              select(logFC) %>%
                              mutate(rank = row_number()) %>%
                              spread(key = rank, value = logFC) %>%
                              filter(`1` * `2` < 0)
Common1 %>% select(Feature) %>% distinct() %>% nrow()
round(100 * 130 / 238, 2)
# 130 genes are common

Not_common <- anti_join(DE1, OE_DE1)
Other_contrasts <- inner_join(Not_common, OE_DE) %>%
                                  filter(str_detect(Contrast, "GBC") |
                                         str_detect(Contrast, "iSUS") |
                                         str_detect(Contrast, "2-.HBC1$")) 
Other_contrasts %>% select(Feature) %>% distinct() %>% nrow()
# 24 genes are in other contrasts
round(100 * 24 / 238, 2)

Not_Detected <- anti_join(Not_common, OE_DE)
Not_Detected %>% select(Feature) %>% distinct() %>% nrow()
# 59 genes are not detected
round(100 * 59 / 238, 2)

Weird_Detected <- inner_join(Not_common, OE_DE) %>%
                  anti_join(Other_contrasts, by = c("Feature" = "Feature"))
Weird_Detected %>% select(Feature) %>% distinct() %>% nrow()
# 25 genes are detected in contrasts further away
round(100 * 25 / 238, 2)

# 2nd branching ----
DE2 <- earlyDETest(models = gamList, knots = c(5, 7),
                   nPoints = 50, pairwise = T)

DE2 <- DE2 %>% mutate(Feature = rownames(DE2)) %>%
               select(Feature, pvalue_1vs2) %>%
               mutate(p.adj = p.adjust(pvalue_1vs2)) %>%
               filter(p.adj < 0.05)
# 126 genes
OE_DE2 <- OE_DE %>% filter(str_detect(Contrast, "GBC") &
                             !str_detect(Contrast, "HBC2"))
Common <- inner_join(DE2, OE_DE2)
Not_consistent <- Common %>% group_by(Feature) %>%
                             filter(n() > 1) %>%
                             select(logFC) %>%
                             mutate(rank = row_number()) %>%
                             spread(key = rank, value = logFC) %>%
                             filter(`1` * `2` > 0)
Common %>% select(Feature) %>% distinct() %>% nrow()
# 58 are common
round(100 * 58 / 126, 2)

Not_common <- anti_join(DE2, OE_DE2)
Other_contrasts <- inner_join(Not_common, OE_DE) %>%
                   filter(str_detect(Contrast, "INP2-INP1") |
                          str_detect(Contrast, "MV1-MV2") |
                          str_detect(Contrast, "^GBC-.HBC2$")) 
Other_contrasts %>% select(Feature) %>% distinct() %>% nrow()
# 42 are in close contrasts
round(100 * 42 / 126, 2)

Not_Detected <- anti_join(Not_common, OE_DE)
Not_Detected %>% select(Feature) %>% distinct() %>% nrow()
# 13 are not detected
round(100 * 13 / 126, 2)

Weird_Detected <- inner_join(Not_common, OE_DE) %>%
                  anti_join(Other_contrasts, by = c("Feature" = "Feature"))
Weird_Detected %>% select(Feature) %>% distinct() %>% nrow()
# 13 are detected in contrasts further away.
round(100 * 13 / 126, 2)