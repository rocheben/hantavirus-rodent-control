## ---------------------------------------------------------------------------
## Figure S6 -- Same as Fig. S5, for the Mexican community. Panels A and B are
## the panels C and D of Fig. 3.
## ---------------------------------------------------------------------------

source("code/R/paths.R"); source("code/R/theme.R"); source("code/R/figures.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork) })

fig <- seov_calendar_figure("Mexico")
ggsave("figures/FigS6_control_calendars_mexico.png", fig, width = 11, height = 10, dpi = 300)
