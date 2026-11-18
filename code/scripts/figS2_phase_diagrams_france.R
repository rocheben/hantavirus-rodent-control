## ---------------------------------------------------------------------------
## Figure S2 -- Phase diagrams (infected individuals against total community
## abundance) during the 30th year in the Lyon community, for three timings
## of targeted control and four control intensities, in a balanced community.
## ---------------------------------------------------------------------------

source("code/R/paths.R"); source("code/R/theme.R"); source("code/R/figures.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

fig <- seov_phase_figure("France")
ggsave("figures/FigS2_phase_diagrams_france.png", fig, width = 10, height = 6.5, dpi = 300)
