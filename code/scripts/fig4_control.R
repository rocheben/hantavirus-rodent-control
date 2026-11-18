## ---------------------------------------------------------------------------
## Figure 4 -- Number of infected individuals during the 30th year, for
## different initial relative abundances of commensal species and different
## intensities of rodent control: control targeted at the competent hosts and
## non-targeted control, in the Lyon community, where SEOV circulates in
## the two Rattus species (A, B), and in the Mexican community, where SNV
## circulates in the two Peromyscus species (C, D). 
##
## The four panels use the same control calendar, the autumn campaign from day
## 181 described in the Methods, so that the two strategies and the two
## communities are directly comparable, and the same colour scale. The four
## combinations of strategy and calendar are in Figs. S5 and S6.
##
## The red segment marks the lowest intensity at which the infection is
## eradicated. It is not a monotone threshold: the burden can rise again above
## it (non-targeted control at the highest intensities).
## ---------------------------------------------------------------------------

source("code/R/paths.R"); source("code/R/theme.R"); source("code/R/figures.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork) })

sim <- bind_rows(seov_control_grid("France"), seov_control_grid("Mexico"))
lims <- range(log10(sim$infected_year + 1))

panel <- function(cm, target, title)
  seov_control_panel(sim, target, AUTUMN_CAMPAIGN, title,
                     fill_limits = lims, keep_community = cm)

fig4 <- (panel("France", "competent", "A. Lyon (SEOV), targeted control (Rattus spp.)") |
         panel("France", "all",       "B. Lyon (SEOV), non-targeted control")) /
        (panel("Mexico", "competent", "C. Mexico (SNV), targeted control (Peromyscus spp.)") |
         panel("Mexico", "all",       "D. Mexico (SNV), non-targeted control")) +
        plot_layout(guides = "collect") &
        theme(legend.position = "right")

ggsave("figures/Fig4_control.png", fig4, width = 11, height = 10, dpi = 300)

## tables quoted in the Results ---------------------------------------------
burden <- sim %>%
  transmute(community, campaign, control_target, control_power, ratio,
            infected_year   = round(infected_year, 1),
            infected_max    = round(infected_max, 2),
            shannon_mean    = round(shannon_mean, 4),
            shannon_at_peak = round(shannon_at_peak, 4)) %>%
  arrange(community, campaign, control_target, control_power, ratio)
write.csv(burden, "results/control_burden.csv", row.names = FALSE)

thr <- seov_eradication_thresholds(sim) %>%
  select(community, campaign, control_target, ratio, eradication_threshold) %>%
  arrange(community, campaign, control_target, ratio)
write.csv(thr, "results/eradication_thresholds.csv", row.names = FALSE)
print(as.data.frame(thr))
