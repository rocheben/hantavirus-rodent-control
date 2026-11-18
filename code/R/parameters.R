## ---------------------------------------------------------------------------
## Parameters of the multi-host SEOV transmission model
##
## Bouilloud, Canet, Vigueras-Galvan, Arnal, Charbonnel & Roche --
## "Contributions of epidemiological modelling to the study of rodent
## management with reciprocal benefits for public health and biodiversity"
##
## Nothing is hard-coded here that belongs to a host community: the species,
## their guild, their competence and their life-history traits are read from
## data/species_traits.csv (built by code/scripts/00_species_traits.R from
## PanTHERIA), and the contact rates from data/contact_matrix_<community>.csv.
## Two communities are defined, each with the hantavirus that circulates in it:
##   "France" -- the small mammal community of an urban park in Lyon, where
##               Seoul orthohantavirus (SEOV) circulates in the two Rattus
##               species;
##   "Mexico" -- a peri-urban community of northeastern Mexico, where Sin
##               Nombre orthohantavirus (SNV) circulates in two Peromyscus
##               species, so that the competent hosts are wild species and not
##               commensal ones.
## ---------------------------------------------------------------------------

SEOV_COMMUNITIES <- c("France", "Mexico")

## Species of one community, in the order of the contact matrix --------------
seov_species <- function(community = "France") {
  community <- match.arg(community, SEOV_COMMUNITIES)
  sp <- read.csv("data/species_traits.csv", stringsAsFactors = FALSE)
  sp <- sp[sp$community == community, ]
  if (nrow(sp) == 0) stop("no species for community ", community)
  ## the two competent hosts of the virus: the Rattus species in the Lyon
  ## community, the Peromyscus species in the Mexican one
  sp$is_competent <- sp$competence %in% c("primary", "secondary")
  rownames(sp) <- NULL
  sp
}

## Contact matrix C ----------------------------------------------------------
## C[i, j] is the contact rate between species i and j. Intra-specific contacts
## are set to 1; inter-specific contacts are low (0.001), moderate (0.01) or
## relatively high (0.1) depending on habitat sharing and agonistic behaviour.
## The matrix is symmetric by construction and its rows are in the order of
## data/species_traits.csv; data/contact_matrix_rationale.csv documents each
## value.
seov_contact_matrix <- function(community = "France") {
  community <- match.arg(community, SEOV_COMMUNITIES)
  path <- sprintf("data/contact_matrix_%s.csv", tolower(community))
  raw  <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  C    <- as.matrix(raw[, -1, drop = FALSE])
  rownames(C) <- raw$species
  codes <- seov_species(community)$code
  if (!setequal(codes, rownames(C)) || !setequal(codes, colnames(C))) {
    stop("the contact matrix of ", community, " does not match its species list")
  }
  C <- C[codes, codes, drop = FALSE]
  stopifnot(isSymmetric(C), all(diag(C) == 1))
  C
}

## Master parameter list -----------------------------------------------------
## community      : "France" or "Mexico"; selects the species and the contacts.
## tau_primary    : susceptibility of the primary reservoir (R. norvegicus for
##                  SEOV, P. maniculatus for SNV).
## tau_secondary_f: the secondary reservoir is assumed slightly less
##                  susceptible, 80% of the primary one. The same two values
##                  are used in both applications, so that the two communities
##                  differ by their composition and not by the virus
##                  parameters.
## K_ref          : reference carrying capacity, identical for all species.
## N_per_species  : total community size / number of species (Ntot = 120).
## sigma, alpha   : recovery rate and disease-induced mortality; both 0 in this
##                  study, so the SIR skeleton implemented in R/model.R reduces
##                  to the SI model described in the manuscript.
## control_*      : rodent control is applied over `control_duration` days from
##                  day `control_start` of every simulated year, but only after
##                  `control_burn_in` days (10 years) so that the community has
##                  reached its demographic equilibrium first.
## breeding_days  : length of the high-reproduction half-year.
## seasonal_anchor: how the reproductive pulses are placed in time.
##   "calendar" (default) anchors the pulses on the day of the year, so the
##      reproductive calendar is identical in every simulated year. This is
##      what the Methods describe and what every figure of the paper uses.
##   "absolute" indexes them on absolute simulation time instead, which is what
##      the original notebooks did and which makes the calendar drift from year
seov_parameters <- function(community      = "France",
                            n_years        = 40,
                            tau_primary    = 0.13,
                            tau_secondary_f = 0.80,
                            K_ref          = 100,
                            N_per_species  = 20,
                            sigma          = 0,
                            alpha          = 0,
                            control_duration = 30,
                            control_burn_in  = 3650,
                            breeding_days    = 182,
                            seasonal_anchor  = c("calendar", "absolute")) {

  community       <- match.arg(community, SEOV_COMMUNITIES)
  seasonal_anchor <- match.arg(seasonal_anchor)

  sp <- seov_species(community)
  m  <- nrow(sp)

  ## two competent species per community, a primary and a secondary reservoir
  tau <- rep(0, m)
  tau[sp$competence == "primary"]   <- tau_primary
  tau[sp$competence == "secondary"] <- tau_primary * tau_secondary_f
  if (sum(tau > 0) != 2) stop("community ", community,
                              " must have one primary and one secondary host")

  ## the species carrying the `ratio` weight at t = 0 are the two commensal
  ## species that dominate the trapping records of the study sites
  focal <- c("Rattus norvegicus", "Mus musculus")

  list(
    community = community,
    virus     = unique(sp$virus),
    m       = m,
    species = sp,
    C       = seov_contact_matrix(community),
    tau     = tau,
    focal_species = focal,
    LS      = sp$litter_size,
    LPY     = sp$litters_per_year,
    Gesta   = sp$gestation_days,
    ## intrinsic death rate, calibrated on the mean intrinsic birth rate
    d       = (sp$litter_size * sp$litters_per_year) / (365 * 2),
    K_ref   = K_ref,
    Ntot0   = N_per_species * m,
    sigma   = sigma,
    alpha   = alpha,
    control_duration = control_duration,
    control_burn_in  = control_burn_in,
    breeding_days    = breeding_days,
    seasonal_anchor  = seasonal_anchor,
    n_years = n_years,
    times   = seq_len(365 * n_years)
  )
}
