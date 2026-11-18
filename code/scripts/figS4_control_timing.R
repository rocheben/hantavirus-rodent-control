## ---------------------------------------------------------------------------
## Figure S4 -- Infection burden of the 30th year for each combination of
## timing and intensity of targeted control, in a balanced community, for the
## two communities.
## ---------------------------------------------------------------------------

source("code/R/paths.R"); source("code/R/parameters.R"); source("code/R/theme.R"); source("code/R/figures.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

burden <- bind_rows(lapply(SEOV_COMMUNITIES, function(cm)
  seov_read_sim("timing", cm))) %>%
  mutate(start_label = factor(control_start, levels = c(1, 91, 181),
                              labels = SEOV_TIMING_LEVELS),
         panel = ifelse(community == "France", "A. Lyon (SEOV)",
                        "B. Mexico (SNV)")) %>%
  group_by(panel, start_label, control_power) %>%
  summarise(infected_year = sum(Itot), .groups = "drop")

figS4 <- ggplot(burden, aes(x = factor(control_power), y = start_label,
                            fill = log10(infected_year + 1))) +
  geom_tile(colour = "grey60", linewidth = 0.6) +
  geom_text(aes(label = round(infected_year)), size = 3.1) +
  facet_wrap(~ panel, ncol = 1) +
  scale_y_discrete(limits = rev(SEOV_TIMING_LEVELS)) +
  scale_fill_gradient(low = "#FFF4D8", high = "darkred",
                      name = "Infected individuals\nduring the 30th year\n(log10 scale)") +
  labs(x = expression(paste("Control intensity (", day^-1, ")")), y = NULL) +
  seov_theme() +
  theme(panel.grid = element_blank())

ggsave("figures/FigS4_control_timing.png", figS4, width = 8, height = 6, dpi = 300)
write.csv(burden %>% mutate(infected_year = round(infected_year, 2)),
          "results/timing_burden.csv", row.names = FALSE)
print(as.data.frame(burden))
