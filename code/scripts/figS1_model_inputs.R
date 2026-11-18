## ---------------------------------------------------------------------------
## Figure S1 -- Model inputs: seasonal birth rate of the two rat species and
## the two control calendars, over two consecutive years taken after the
## burn-in. The figure shows that the reproductive calendar repeats identically
## from one year to the next.
##
## Also writes results/birth_death_balance.csv, the birth / death balance of
## every species of both communities.
## ---------------------------------------------------------------------------

source("code/R/paths.R"); source("code/R/parameters.R"); source("code/R/model.R")
source("code/R/simulate.R"); source("code/R/theme.R")
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(tidyr) })

p    <- seov_parameters("France")
days <- 3650 + seq_len(2 * 365)

birth_curve <- function(latin) {
  i <- which(p$species$latin == latin)
  data.frame(day   = days - 3650,
             value = vapply(days, function(t) seov_birth_rate(t, p)[i], numeric(1)),
             panel = paste0("Birth rate of ", latin))
}

control_curve <- function(start) {
  q <- p
  q$control_power <- seov_control_power(0.07, "competent", q)
  q$control_start <- start
  data.frame(day   = days - 3650,
             value = vapply(days, function(t) seov_control_rate(t, q)[1], numeric(1)),
             panel = ifelse(start == 181,
                            "Autumn campaign (day 181), control at 0.07",
                            "Spring campaign (day 1), control at 0.07"))
}

levels_panel <- c("Birth rate of Rattus norvegicus", "Birth rate of Rattus rattus",
                  "Spring campaign (day 1), control at 0.07",
                  "Autumn campaign (day 181), control at 0.07")

figS1 <- bind_rows(birth_curve("Rattus norvegicus"), birth_curve("Rattus rattus"),
                   control_curve(1), control_curve(181)) %>%
  mutate(panel = factor(panel, levels = levels_panel)) %>%
  ggplot(aes(x = day, y = value)) +
  geom_line(linewidth = 0.4) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  labs(x = "Day of simulation (two years, after the 10-year burn-in)",
       y = expression(paste("Rate (", day^-1, ")"))) +
  seov_theme()

ggsave("figures/FigS1_model_inputs.png", figS1, width = 8, height = 7, dpi = 300)

## birth / death balance of both communities; the calendar anchor makes it
## year-independent, which is checked here
balance <- lapply(SEOV_COMMUNITIES, function(cm) {
  q <- seov_parameters(cm)
  stopifnot(all.equal(seov_mean_birth_rate(q, 1), seov_mean_birth_rate(q, 30)))
  data.frame(community        = cm,
             species          = q$species$code,
             latin            = q$species$latin,
             pulses_per_year  = seov_pulses_per_year(q),
             litters_per_year = q$LPY,
             mean_birth_rate  = round(seov_mean_birth_rate(q, 1), 5),
             death_rate       = round(q$d, 5))
}) %>% bind_rows() %>%
  mutate(birth_over_death = round(mean_birth_rate / death_rate, 3))

write.csv(balance, "results/birth_death_balance.csv", row.names = FALSE)
print(balance)
