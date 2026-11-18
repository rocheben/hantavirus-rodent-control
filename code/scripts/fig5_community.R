## ---------------------------------------------------------------------------
## Figure 5 -- Abundance of every species and number of infected individuals
## during the 30th year under three intensities of non-targeted control, in a
## balanced community (ratio = 1) with the autumn campaign (from day 181), for
## the Lyon (A) and the Mexican (B) community.
##
## Fig. S7 is the same figure with the spring campaign, Fig. S8 the same figure
## under targeted control.
## ---------------------------------------------------------------------------

source("code/R/paths.R"); source("code/R/parameters.R"); source("code/R/model.R")
source("code/R/simulate.R"); source("code/R/theme.R"); source("code/R/figures.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(patchwork) })

fig5 <- seov_community_pair("all", AUTUMN_CAMPAIGN)
ggsave("figures/Fig5_community.png", fig5, width = 11, height = 9, dpi = 300)

## summary quoted in the Results, for every panel of Figs. 4, S7 and S8 -------
summary_tbl <- bind_rows(lapply(SEOV_COMMUNITIES, function(cm)
  seov_community_summary(seov_read_sim("community", cm), seov_parameters(cm))))
write.csv(summary_tbl, "results/community_summary.csv", row.names = FALSE)

print(as.data.frame(filter(summary_tbl, control_target == "all",
                           control_start == AUTUMN_CAMPAIGN) %>%
                      select(community, control_power, species, abundance_min,
                             end_vs_no_control, infected_year, infected_max,
                             shannon_at_peak)))
