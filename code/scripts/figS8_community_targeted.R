## ---------------------------------------------------------------------------
## Figure S8 -- Same as Fig. 5 under targeted control: the culling is
## restricted to the two competent hosts of the virus (Rattus spp. in the
## Lyon community, Peromyscus spp. in the Mexican one), with the campaign
## applied from day 181.
## ---------------------------------------------------------------------------

source("code/R/paths.R"); source("code/R/parameters.R"); source("code/R/model.R")
source("code/R/simulate.R"); source("code/R/theme.R"); source("code/R/figures.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })

fig <- seov_community_pair("competent", AUTUMN_CAMPAIGN)
ggsave("figures/FigS8_community_targeted.png", fig, width = 11, height = 9, dpi = 300)
