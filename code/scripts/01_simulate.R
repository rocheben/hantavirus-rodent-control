## ---------------------------------------------------------------------------
## 01 -- Run every simulation of the paper and cache the results.
##
## Four experiments, run identically for both communities:
##   A. dilution       -- prevalence vs diversity, no control         (Fig. 3)
##   B. control grid   -- targeted vs non-targeted control (Figs. 4, S2, S3, S5, S6)
##   C. community      -- community dynamics under control (Figs. 5, S7, S8)
##   D. control timing -- effect of the timing of control             (Fig. S4)
##
## Output: results/sim_<experiment>_<community>.rds. Every figure script reads
## one of these files and produces exactly one figure, so nothing has to be
## re-simulated to redraw a panel.
##
## Conventions used throughout the project:
##   * time starts at t = 1 and the n-th year is made of the 365 time steps
##     (n - 1) * 365 + 1 ... n * 365, so year 30 is days 10 586 to 10 950;
##   * the day of the year runs from 1 to 365 (seov_day_of_year());
##   * a control campaign covers [control_start, control_start + 30), i.e.
##     exactly 30 days, and starts only after the 10-year burn-in;
##   * "infected individuals" always means every infected individual of the
##     community, i.e. both competent hosts (column Itot).
##
## Run from the root of the project: Rscript code/scripts/01_simulate.R (15-30 min)
## ---------------------------------------------------------------------------

source("code/R/paths.R")
source("code/R/parameters.R")
source("code/R/model.R")
source("code/R/simulate.R")
suppressPackageStartupMessages({
  library(deSolve); library(vegan); library(dplyr); library(tidyr)
})

## Analysis windows ----------------------------------------------------------
EQUILIBRIUM_DAY <- 8000   # first day used for the prevalence/diversity summary
FOCAL_YEAR      <- 30     # year reported in Figs. 4, 5 and S1 to S8

## First day of the 30-day control campaign. Days 1-182 of the model year are
## the high-reproduction season and days 183-365 the low-reproduction season,
## so a campaign starting on day 181 is applied at the end of the breeding
## season (autumn campaign) and one starting on day 1 at its very beginning
## (spring campaign). Both calendars are simulated for both strategies.
AUTUMN_CAMPAIGN <- 181
SPRING_CAMPAIGN <- 1

## Control intensities. The grid is regular (step 0.021) up to 0.147 and is
## extended to 0.16, the highest intensity of the original analysis, so that
## no eradication threshold is hidden by a gap in the grid.
CONTROL_POWERS <- c(seq(0, 0.147, by = 0.021), 0.16)

## Control strategies simulated in both communities: "competent" culls the two
## competent hosts of the virus, "all" every species of the community. In the
## Lyon community the competent hosts are the two Rattus species, in the
## Mexican one the two Peromyscus species.
TARGETS <- c("competent", "all")

## ---------------------------------------------------------------------------
## The four experiments, for one community
## ---------------------------------------------------------------------------
simulate_community <- function(community) {

  p <- seov_parameters(community)
  tag <- tolower(community)
  message("\n#### community: ", community, " ####")

  ## A. dilution effect -------------------------------------------------------
  keep_equilibrium <- function(sim) {
    sim %>%
      filter(time > EQUILIBRIUM_DAY) %>%
      transmute(time, year, community, ratio, Itot, Ires, Nres, Ntot, Shannon)
  }
  message("== A. dilution effect ==")
  sim_dilution <- seov_run_grid(expand.grid(ratio = c(0.5, 1, 5)), p, keep_equilibrium)
  saveRDS(sim_dilution, sprintf("results/sim_dilution_%s.rds", tag))

  ## B. control grid ----------------------------------------------------------
  grid_control <- expand.grid(
    ratio          = c(0.8, 1, 10),
    control_power  = CONTROL_POWERS,
    control_target = TARGETS,
    control_start  = c(AUTUMN_CAMPAIGN, SPRING_CAMPAIGN),
    stringsAsFactors = FALSE
  )
  ## with no control the scenario depends neither on the target nor on the
  ## calendar: keep a single run per abundance ratio
  grid_control <- grid_control[!(grid_control$control_power == 0 &
                                 !(grid_control$control_target == "competent" &
                                   grid_control$control_start == AUTUMN_CAMPAIGN)), ]

  burden_focal_year <- function(sim) {
    seov_year(sim, FOCAL_YEAR) %>%
      summarise(community = first(community),
                ratio = first(ratio),
                control_power   = first(control_power),
                control_target  = first(control_target),
                control_start   = first(control_start),
                infected_year   = sum(Itot),
                infected_max    = max(Itot),
                shannon_mean    = mean(Shannon),
                shannon_at_peak = Shannon[which.max(Itot)])
  }
  message("== B. rodent control grid ==")
  sim_control <- seov_run_grid(grid_control, p, burden_focal_year)

  ## replicate the no-control row across targets and calendars, so that every
  ## panel shows its own baseline
  baseline <- sim_control %>% filter(control_power == 0)
  sim_control <- bind_rows(
    sim_control %>% filter(control_power > 0),
    expand.grid(control_target = TARGETS,
                control_start  = c(AUTUMN_CAMPAIGN, SPRING_CAMPAIGN),
                stringsAsFactors = FALSE) %>%
      tidyr::crossing(baseline %>% select(-control_target, -control_start))
  )
  saveRDS(sim_control, sprintf("results/sim_control_%s.rds", tag))

  ## C. community dynamics ----------------------------------------------------
  grid_community <- expand.grid(
    ratio          = 1,
    control_power  = c(0, 0.1, 0.2),
    control_target = TARGETS,
    control_start  = c(AUTUMN_CAMPAIGN, SPRING_CAMPAIGN),
    stringsAsFactors = FALSE
  )
  grid_community <- grid_community[!(grid_community$control_power == 0 &
                                     grid_community$control_target != "all"), ]

  focal_year_long <- function(sim) {
    seov_long(seov_year(sim, FOCAL_YEAR), p) %>%
      transmute(day = seov_day_of_year(time), community, ratio, control_power,
                control_target, control_start, species, guild, Abundance, I,
                Itot, Ires, Nres, Ntot, Shannon, phase)
  }
  message("== C. community dynamics ==")
  sim_community <- seov_run_grid(grid_community, p, focal_year_long)
  baseline_com <- sim_community %>%
    filter(control_power == 0, control_target == "all")
  others <- setdiff(TARGETS, "all")
  sim_community <- bind_rows(sim_community,
    do.call(rbind, lapply(others, function(tg)
      baseline_com %>% mutate(control_target = tg))))
  saveRDS(sim_community, sprintf("results/sim_community_%s.rds", tag))

  ## D. timing of targeted control --------------------------------------------
  {
    grid_timing <- expand.grid(
      ratio          = 1,
      control_power  = seq(0, 0.105, by = 0.035),
      control_start  = c(1, 91, 181),
      control_target = "competent",
      stringsAsFactors = FALSE
    )
    focal_year_rat <- function(sim) {
      seov_year(sim, FOCAL_YEAR) %>%
        transmute(day = seov_day_of_year(time), community, ratio, control_power,
                  control_start, Itot, Ires, Nres, Ntot, Shannon, phase)
    }
    message("== D. timing of control ==")
    sim_timing <- seov_run_grid(grid_timing, p, focal_year_rat)
    saveRDS(sim_timing, sprintf("results/sim_timing_%s.rds", tag))
  }
}

for (cm in SEOV_COMMUNITIES) simulate_community(cm)

message("\ndone -- simulations cached in results/")
