## ---------------------------------------------------------------------------
## Reproduce every figure and table of the paper, for both communities.
##
##   cd <root of this project>
##   Rscript run_all.R
##
## Every path used by the project is relative to this directory, so the code
## runs unchanged on any machine. Simulations are cached in results/; delete
## that folder to force a re-run (15-30 minutes for the two communities).
##
## One script per figure: once code/scripts/01_simulate.R has been run, any figure
## can be redrawn on its own, for instance
##   Rscript code/scripts/fig4_control.R
##
## Required R packages: deSolve, vegan, dplyr, tidyr, ggplot2, patchwork.
## ---------------------------------------------------------------------------

source("code/R/paths.R")

scripts <- c(
  "code/scripts/00_species_traits.R",
  "code/scripts/01_simulate.R",
  "code/scripts/fig1_model_diagram.R",
  "code/scripts/fig2_communities.R",
  "code/scripts/fig3_dilution.R",
  "code/scripts/fig4_control.R",
  "code/scripts/fig5_community.R",
  "code/scripts/figS1_model_inputs.R",
  "code/scripts/figS2_phase_diagrams_france.R",
  "code/scripts/figS3_phase_diagrams_mexico.R",
  "code/scripts/figS4_control_timing.R",
  "code/scripts/figS5_control_calendars_france.R",
  "code/scripts/figS6_control_calendars_mexico.R",
  "code/scripts/figS7_community_spring.R",
  "code/scripts/figS8_community_targeted.R",
  "code/scripts/table1_parameters.R",
  "code/scripts/supp_trapping_success.R"
)

for (s in scripts) {
  message("\n>>> ", s)
  source(s, local = new.env(), echo = FALSE)
}

message("\nFigures in figures/, tables in results/")
print(sessionInfo())
