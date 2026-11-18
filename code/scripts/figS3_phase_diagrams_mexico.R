## ---------------------------------------------------------------------------
## Figure S3 -- Same as Fig. S2, for the Mexican community.
## ---------------------------------------------------------------------------

source("code/R/paths.R"); source("code/R/theme.R"); source("code/R/figures.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

fig <- seov_phase_figure("Mexico")
ggsave("figures/FigS3_phase_diagrams_mexico.png", fig, width = 10, height = 6.5, dpi = 300)
