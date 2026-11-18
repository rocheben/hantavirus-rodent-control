## ---------------------------------------------------------------------------
## Figure 3 -- Hantavirus prevalence as a function of small mammal diversity:
## SEOV in the Lyon community (A), SNV in the Mexican community (B).
##
## The thin paths are the stationary trajectories, one per simulated year after
## day 8000: every year the community travels the same loop in the
## prevalence-diversity plane. The points are the local extrema of the number
## of infected individuals, the quantity summarised in the Results.
##
## "Infected individuals" means every infected individual of the community,
## exactly as in Figs. 4 and 5. The prevalence among the competent hosts alone,
## which is what field seroprevalences refer to, is written to
## results/prevalence_in_hosts.csv.
## ---------------------------------------------------------------------------

source("code/R/paths.R"); source("code/R/parameters.R"); source("code/R/model.R")
source("code/R/simulate.R"); source("code/R/theme.R"); source("code/R/figures.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork) })

france <- seov_read_sim("dilution", "France")
mexico <- seov_read_sim("dilution", "Mexico")

fig3 <- seov_dilution_panel(france, "A. Lyon (SEOV)") +
        seov_dilution_panel(mexico, "B. Mexico (SNV)") +
        plot_layout(guides = "collect") &
        theme(legend.position = "right")

ggsave("figures/Fig3_dilution.png", fig3, width = 10, height = 4.4, dpi = 300)

## values quoted in the Results
summary_tbl <- bind_rows(seov_dilution_summary(seov_dilution_extrema(france)),
                         seov_dilution_summary(seov_dilution_extrema(mexico)))
write.csv(summary_tbl, "results/fig3_dilution.csv", row.names = FALSE)

## prevalence among the competent hosts, over the whole stationary window
hosts <- bind_rows(france, mexico) %>%
  group_by(community, ratio) %>%
  summarise(prev_hosts_min  = min(Ires / Nres),
            prev_hosts_mean = mean(Ires / Nres),
            prev_hosts_max  = max(Ires / Nres),
            prev_com_max    = max(Itot / Ntot),
            .groups = "drop")
write.csv(hosts, "results/prevalence_in_hosts.csv", row.names = FALSE)

print(as.data.frame(summary_tbl))
print(as.data.frame(hosts))
