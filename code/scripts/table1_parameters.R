## ---------------------------------------------------------------------------
## Table 1 -- every parameter of the model, with its value and its source.
## Table S1 -- life-history traits of the two communities and their provenance.
## Table S2 -- the two contact matrices, in long format.
## ---------------------------------------------------------------------------

source("code/R/paths.R"); source("code/R/parameters.R"); source("code/R/model.R"); source("code/R/simulate.R")
suppressPackageStartupMessages({ library(dplyr) })

pf <- seov_parameters("France")
pm <- seov_parameters("Mexico")
bd <- c(seov_mean_birth_rate(pf, 1) / pf$d, seov_mean_birth_rate(pm, 1) / pm$d)

row <- function(symbol, description, value, unit, source)
  data.frame(symbol = symbol, description = description, value = value,
             unit = unit, source = source, stringsAsFactors = FALSE)

table1 <- rbind(
  row("$m$", "Number of host species per community", as.character(pf$m), "-",
      "This study"),
  row("$N_{tot}(0)$", "Initial community abundance", as.character(pf$Ntot0),
      "individuals",
      "Empirical data (Hussein et al. 2023; Vigueras-Galvan et al. 2019)"),
  row("$K_i$", "Carrying capacity (reference, scaled by sqrt(ratio))",
      as.character(pf$K_ref), "individuals", "This study"),
  row("$LS_i$", "Litter size", sprintf("%.2f-%.2f", min(c(pf$LS, pm$LS)), max(c(pf$LS, pm$LS))),
      "pups", "PanTHERIA (Jones et al. 2009), see Table S1"),
  row("$LPY_i$", "Litters per year",
      sprintf("%.2f-%.2f", min(c(pf$LPY, pm$LPY)), max(c(pf$LPY, pm$LPY))),
      "year^-1", "PanTHERIA (Jones et al. 2009), see Table S1"),
  row("$Gesta_i$", "Gestation length",
      sprintf("%d-%d", as.integer(min(c(pf$Gesta, pm$Gesta))),
              as.integer(max(c(pf$Gesta, pm$Gesta)))),
      "days", "PanTHERIA (Jones et al. 2009), see Table S1"),
  row("$d_i$", "Intrinsic death rate, from litter size and litters per year",
      paste(sprintf("%.4f", range(c(pf$d, pm$d))), collapse = "-"), "day^-1",
      "Calibrated on the mean intrinsic birth rate"),
  row("$\\bar{b}_i / d_i$", "Ratio of the mean intrinsic birth rate to the death rate",
      paste(sprintf("%.2f", range(bd)), collapse = "-"), "-",
      "Computed from the life-history traits"),
  row("$C_{i,j}$", "Contact rate: intra-specific / high / moderate / low",
      "1 / 0.1 / 0.01 / 0.001", "-",
      "Habitat sharing and agonistic behaviour, see Table S2"),
  row("$\\tau_i$", "Susceptibility (primary reservoir / secondary / others)",
      paste(c(sprintf("%.3f", c(0.13, 0.13 * 0.8)), "0"), collapse = " / "), "-",
      "Chosen so that prevalence stays in the observed range"),
  row("$\\sigma, \\alpha$", "Recovery rate and virulence (SI model)", "0", "day^-1",
      "This study"),
  row("$\\rho_i$", "Duration of the control campaign", as.character(pf$control_duration),
      "days", "This study"),
  row("", "First day of the control campaign", "1 or 181", "day of the year",
      "Spring and autumn calendars, both simulated"),
  row("", "Burn-in before the first control campaign",
      as.character(pf$control_burn_in), "days",
      "Demographic equilibrium reached after 10 years"),
  row("", "Simulated horizon", as.character(365 * pf$n_years), "days",
      "40 years, daily time step")
)
write.csv(table1, "results/table1_parameters.csv", row.names = FALSE)

## Table S1: traits and their provenance, in a compact form. A value that is
## not the species value of PanTHERIA carries a mark: (g) for the mean of the
## genus, (a) for the value declared in data/trait_overrides.csv.
traits <- read.csv("data/species_traits.csv", stringsAsFactors = FALSE)

mark <- function(value, source) {
  tag <- ifelse(grepl("genus mean", source), " (g)",
                ifelse(source == "PanTHERIA", "", " (a)"))
  paste0(format(value, trim = TRUE), tag)
}

tableS1 <- traits %>%
  transmute(community, virus, species = latin, guild,
            competence = ifelse(competence == "none", "-", competence),
            `litter size`      = mark(litter_size, litter_size_source),
            `litters per year` = mark(litters_per_year, litters_per_year_source),
            `gestation (days)` = mark(gestation_days, gestation_days_source))
write.csv(tableS1, "results/tableS1_species_traits.csv", row.names = FALSE)

## Table S2: the two contact matrices in long format
long_matrix <- function(cm) {
  C <- seov_contact_matrix(cm)
  data.frame(community = cm,
             species_i = rep(rownames(C), times = ncol(C)),
             species_j = rep(colnames(C), each = nrow(C)),
             contact   = as.vector(C), stringsAsFactors = FALSE)
}
tableS2 <- bind_rows(lapply(SEOV_COMMUNITIES, long_matrix)) %>%
  filter(species_i <= species_j)
write.csv(tableS2, "results/tableS2_contact_matrices.csv", row.names = FALSE)

print(table1[, c("symbol", "value", "unit")])
print(as.data.frame(tableS1))
