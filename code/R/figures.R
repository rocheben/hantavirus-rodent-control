## ---------------------------------------------------------------------------
## Building blocks shared by several figures
##
## The figures of the main text put the Lyon and the Mexican community side
## by side as panels of the same figure, and the supplementary figures repeat
## the same plots for other strategies and calendars. Each plot is therefore
## defined once here and every figure script only assembles the panels it
## needs.
## ---------------------------------------------------------------------------

AUTUMN_CAMPAIGN <- 181
SPRING_CAMPAIGN <- 1

seov_campaign_label <- function(start) {
  ifelse(start == AUTUMN_CAMPAIGN, "autumn campaign (day 181)",
         "spring campaign (day 1)")
}

## Cached simulations --------------------------------------------------------
seov_read_sim <- function(experiment, community) {
  path <- sprintf("results/sim_%s_%s.rds", experiment, tolower(community))
  if (!file.exists(path)) {
    stop("missing ", path, " -- run code/scripts/01_simulate.R first")
  }
  readRDS(path)
}

seov_control_grid <- function(community) {
  sim <- seov_read_sim("control", community)
  sim$ratio_f  <- factor(sim$ratio)
  sim$power_f  <- factor(sim$control_power)
  sim$campaign <- seov_campaign_label(sim$control_start)
  sim
}

## Dilution effect (Fig. 3) -------------------------------------------------
## Local extrema of the number of infected individuals, after the community has
## reached its stationary regime, expressed as a prevalence in the community.
seov_dilution_extrema <- function(sim) {
  suppressPackageStartupMessages(library(dplyr))
  sim %>%
    group_by(ratio) %>%
    filter(seov_local_extrema(Itot, span = 400)) %>%
    ungroup() %>%
    mutate(prevalence         = Itot / Ntot,
           prevalence_in_hosts = Ires / Nres,
           ratio_f            = factor(ratio, levels = c(5, 1, 0.5)))
}

## One panel of the dilution figure. The whole stationary trajectory is drawn
## as a thin translucent path, which shows the seasonal loop travelled every
## year in the prevalence-diversity plane, and the local extrema used in the
## Results are drawn on top of it.
seov_dilution_panel <- function(sim, title = NULL, span = 400) {
  suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
  traj <- sim %>%
    mutate(prevalence = Itot / Ntot,
           ratio_f    = factor(ratio, levels = c(5, 1, 0.5)))
  extrema <- seov_dilution_extrema(sim)

  ggplot(traj, aes(x = Shannon, y = prevalence, colour = ratio_f)) +
    geom_path(aes(group = interaction(ratio_f, year)), linewidth = 0.3, alpha = 0.25) +
    geom_point(data = extrema, size = 1.9, alpha = 0.9) +
    scale_colour_manual(values = seov_ratio_colours,
                        name = "Relative abundance\nof commensal species",
                        drop = FALSE) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = title, x = "Shannon index",
         y = "Prevalence in the community") +
    seov_theme() +
    theme(plot.title = element_text(size = 10))
}

seov_dilution_summary <- function(extrema) {
  suppressPackageStartupMessages(library(dplyr))
  extrema %>%
    group_by(community, ratio) %>%
    summarise(n_points       = n(),
              shannon_min    = min(Shannon),
              shannon_max    = max(Shannon),
              prevalence_min = min(prevalence),
              prevalence_max = max(prevalence),
              .groups = "drop")
}

## Control grids (Figs. 4, S5, S6) ------------------------------------------
## Lowest control intensity at which the infection is eradicated during the
## focal year (fewer than one infected individual over the whole year), for
## each combination of strategy, calendar and abundance ratio. Absent when the
## infection is never eradicated on the intensity grid, and absent as well
## where the infection does not establish at all without control, since there
## is then nothing to eradicate.
seov_eradication_thresholds <- function(sim) {
  suppressPackageStartupMessages(library(dplyr))
  established <- sim %>%
    filter(control_power == 0, infected_year >= 1) %>%
    distinct(community, ratio)
  sim %>%
    semi_join(established, by = c("community", "ratio")) %>%
    filter(infected_year < 1) %>%
    group_by(community, campaign, control_target, control_start, ratio, ratio_f) %>%
    summarise(eradication_threshold = min(control_power), .groups = "drop")
}

## One tile panel: infection burden of the focal year for one strategy and one
## calendar. The red segment marks the lowest intensity at which the infection
## is eradicated; it is not a monotone threshold, the burden can rise again
## above it (non-targeted control at high intensity), which is why the panel
## also prints every value.
seov_control_panel <- function(sim, target, start, title, fill_limits = NULL,
                               keep_community = NULL) {
  suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
  if (is.null(fill_limits)) fill_limits <- range(log10(sim$infected_year + 1))
  if (!is.null(keep_community)) sim <- sim[sim$community == keep_community, ]
  d  <- filter(sim, control_target == target, control_start == start)
  th <- seov_eradication_thresholds(sim) %>%
    filter(control_target == target, control_start == start) %>%
    mutate(power_f = factor(eradication_threshold, levels = levels(sim$power_f)))

  p <- ggplot(d, aes(x = ratio_f, y = power_f, fill = log10(infected_year + 1))) +
    geom_tile(colour = "grey60", linewidth = 0.6) +
    geom_text(aes(label = round(infected_year)), size = 3.1, colour = "black") +
    scale_fill_gradient(low = "#FFF4D8", high = "darkred", limits = fill_limits,
                        name = "Infected individuals\nduring the 30th year\n(log10 scale)") +
    labs(title = title,
         x = "Relative initial abundance of commensal species",
         y = expression(paste("Control intensity (", day^-1, ")"))) +
    seov_theme() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 10))

  if (nrow(th) > 0) {
    p <- p + geom_segment(
      data = th,
      aes(x = as.numeric(ratio_f) - 0.5, xend = as.numeric(ratio_f) + 0.5,
          y = as.numeric(power_f) - 0.5, yend = as.numeric(power_f) - 0.5),
      inherit.aes = FALSE, colour = "red", linewidth = 1.2)
  }
  p
}

## Community dynamics (Figs. 5, S7, S8) --------------------------------------
## Abundance of every species during the focal year, with the number of
## infected individuals, for one strategy and one calendar.
seov_community_figure <- function(sim, target, start, p, title = NULL) {
  suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
  d <- filter(sim, control_target == target, control_start == start)
  d$species <- factor(d$species, levels = p$species$code)
  ## the campaign markers only make sense on the controlled panels
  window <- d %>%
    filter(control_power > 0) %>%
    distinct(control_power) %>%
    tidyr::crossing(x = c(start, start + p$control_duration))

  ggplot(d, aes(x = day)) +
    geom_vline(data = window, aes(xintercept = x, linetype = "Control campaign"),
               colour = "red", linewidth = 0.4) +
    geom_line(aes(y = Abundance, colour = species), linewidth = 0.5) +
    geom_point(aes(y = Itot, shape = "Infected individuals"),
               colour = "darkred", size = 0.7) +
    facet_wrap(~ control_power, nrow = 1,
               labeller = labeller(control_power = function(x)
                 paste("Control intensity:", x))) +
    scale_colour_brewer(palette = "Dark2", name = "Species",
                        labels = p$species$latin) +
    scale_shape_manual(values = 16, name = NULL) +
    scale_linetype_manual(values = "dashed", name = NULL) +
    labs(title = title, x = "Day of the 30th year", y = "Abundance") +
    seov_theme() +
    theme(legend.key.height = unit(0.9, "lines"),
          legend.text = element_text(face = "italic"),
          plot.title  = element_text(size = 10))
}

## Summary of the community experiment, for the tables of the Results.
seov_community_summary <- function(sim, p) {
  suppressPackageStartupMessages(library(dplyr))
  community <- sim %>%
    filter(species == p$species$code[1]) %>%
    group_by(community, control_target, control_start, control_power) %>%
    summarise(shannon_mean    = mean(Shannon),
              shannon_min     = min(Shannon),
              shannon_max     = max(Shannon),
              shannon_at_peak = Shannon[which.max(Itot)],
              infected_year   = sum(Itot),
              infected_max    = max(Itot),
              .groups = "drop")

  sim %>%
    group_by(community, control_target, control_start, control_power, species) %>%
    summarise(abundance_min   = min(Abundance),
              abundance_start = first(Abundance),
              abundance_end   = last(Abundance),
              .groups = "drop") %>%
    left_join(community,
              by = c("community", "control_target", "control_start", "control_power")) %>%
    group_by(community, control_target, control_start, species) %>%
    mutate(end_vs_no_control = abundance_end / abundance_end[control_power == 0]) %>%
    ungroup() %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))
}

## Timing of control (Figs. S1 to S3) ----------------------------------------
SEOV_TIMING_LEVELS <- c("Day 1 (before breeding)",
                        "Day 91 (during breeding)",
                        "Day 181 (after breeding)")

seov_phase_figure <- function(community) {
  suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
  sim <- seov_read_sim("timing", community) %>%
    mutate(start_label = factor(control_start, levels = c(1, 91, 181),
                                labels = SEOV_TIMING_LEVELS))
  ggplot(sim, aes(x = Ntot, y = Itot, colour = phase)) +
    geom_point(size = 0.6) +
    facet_grid(start_label ~ control_power,
               labeller = labeller(control_power = function(x)
                 paste("intensity", x))) +
    scale_colour_brewer(palette = "Set1", name = "Phase", drop = FALSE) +
    labs(x = "Total abundance of the community", y = "Infected individuals") +
    seov_theme() +
    guides(colour = guide_legend(override.aes = list(size = 2)))
}

## The four strategy x calendar combinations of one community (Figs. S5, S6).
seov_calendar_figure <- function(community) {
  suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
  sim  <- seov_control_grid(community)
  lims <- range(log10(sim$infected_year + 1))
  panel <- function(target, start, title)
    seov_control_panel(sim, target, start, title, fill_limits = lims)
  (panel("competent", AUTUMN_CAMPAIGN, "A. Targeted, autumn campaign (day 181)") |
   panel("all",       AUTUMN_CAMPAIGN, "B. Non-targeted, autumn campaign (day 181)")) /
  (panel("competent", SPRING_CAMPAIGN, "C. Targeted, spring campaign (day 1)") |
   panel("all",       SPRING_CAMPAIGN, "D. Non-targeted, spring campaign (day 1)")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")
}

## The community dynamics of both communities, stacked (Figs. 4, S7, S8).
seov_community_pair <- function(target, start) {
  suppressPackageStartupMessages(library(patchwork))
  panels <- lapply(SEOV_COMMUNITIES, function(cm) {
    ttl <- if (cm == "France") "A. Lyon (SEOV)" else "B. Mexico (SNV)"
    seov_community_figure(seov_read_sim("community", cm), target, start,
                          seov_parameters(cm), ttl)
  })
  panels[[1]] / panels[[2]]
}
