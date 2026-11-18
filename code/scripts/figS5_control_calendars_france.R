## ---------------------------------------------------------------------------
## Figure S5 -- The control grid of Fig. 3 for the Lyon community, for the
## four combinations of strategy (targeted, non-targeted) and calendar (autumn
## campaign from day 181, spring campaign from day 1). Panels A and B are the
## panels A and B of Fig. 3.
## ---------------------------------------------------------------------------

source("code/R/paths.R"); source("code/R/theme.R"); source("code/R/figures.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork) })

fig <- seov_calendar_figure("France")
ggsave("figures/FigS5_control_calendars_france.png", fig, width = 11, height = 10, dpi = 300)
