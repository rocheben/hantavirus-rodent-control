## ---------------------------------------------------------------------------
## Figure S7 -- Same as Fig. 5 (non-targeted control) with the control campaign
## applied from day 1, at the very beginning of the breeding season, for the
## two communities.
## ---------------------------------------------------------------------------

source("code/R/paths.R"); source("code/R/parameters.R"); source("code/R/model.R")
source("code/R/simulate.R"); source("code/R/theme.R"); source("code/R/figures.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })

fig <- seov_community_pair("all", SPRING_CAMPAIGN)
ggsave("figures/FigS7_community_spring.png", fig, width = 11, height = 9, dpi = 300)
