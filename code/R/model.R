## ---------------------------------------------------------------------------
## Multi-host SEOV model: rates, right-hand side and scenario construction
## ---------------------------------------------------------------------------

## Seasonal birth rate -------------------------------------------------------
## Reproduction is concentrated in the first half of the year (spring-summer,
## 80% of the births) and reduced during the second half (autumn-winter, 20%).
## Within the high-reproduction half-year, births occur as successive pulses of
## length Gesta_i repeated every cycle_i = floor(breeding_days / LPY_i) days,
## each of height 0.8 * LS_i / (2 * Gesta_i); the factor 2 accounts for a
## balanced sex ratio. Outside those pulses the rate is the mean
## low-reproduction rate 0.2 * LS * LPY * (365 - LPY * Gesta) / 365^2.
##
## The pulses are anchored on the day of the year (1 to 365), so the
## reproductive calendar repeats identically every year. Setting
## `seasonal_anchor = "absolute"` in seov_parameters() indexes them on absolute
## simulation time instead, which is what the original notebooks did and which
## makes the calendar drift by 365 %% cycle days per year.
seov_pulse_cycle <- function(p) as.integer(p$breeding_days %/% p$LPY)

seov_birth_rate <- function(t, p) {
  LS <- p$LS; LPY <- p$LPY; G <- p$Gesta
  b <- (0.2 * LS * LPY) * (365 - LPY * G) / 365^2
  cycle <- seov_pulse_cycle(p)
  if (p$seasonal_anchor == "calendar") {
    dy <- seov_day_of_year(t)              # 1 ... 365
    in_pulse <- ((dy - 1) %% cycle > 0) & ((dy - 1) %% cycle < G) &
                (dy <= p$breeding_days)
  } else {
    ## legacy behaviour of the original notebooks: pulses indexed on absolute
    ## simulation time, hence drifting from year to year
    in_pulse <- (t %% cycle > 0) & (t %% cycle < G) & (t %% 365 < p$breeding_days)
  }
  b[in_pulse] <- 0.8 * LS[in_pulse] / (2 * G[in_pulse])
  b
}

## Number of reproductive pulses realised in one calendar year, per species.
## Only meaningful for the "calendar" anchor, where the year repeats itself.
seov_pulses_per_year <- function(p) {
  cycle <- seov_pulse_cycle(p)
  dy    <- 1:365
  vapply(seq_len(p$m), function(i) {
    on <- ((dy - 1) %% cycle[i] > 0) & ((dy - 1) %% cycle[i] < p$Gesta[i]) &
          (dy <= p$breeding_days)
    sum(on & !c(FALSE, on[-length(on)]))
  }, numeric(1))
}

## Mean intrinsic birth rate over one calendar year, per species.
seov_mean_birth_rate <- function(p, year = 1) {
  days <- ((year - 1) * 365 + 1):(year * 365)
  rowMeans(vapply(days, function(d) seov_birth_rate(d, p), numeric(p$m)))
}

## Rodent control (additional mortality) -------------------------------------
## Returns the per-species culling rate applied at time t. The campaign covers
## the half-open interval [control_start, control_start + control_duration),
## i.e. exactly `control_duration` days.
seov_control_rate <- function(t, p) {
  if (t <= p$control_burn_in) return(rep(0, p$m))
  if (p$seasonal_anchor == "calendar") {
    day <- seov_day_of_year(t)
    if (day < p$control_start || day >= p$control_start + p$control_duration) {
      return(rep(0, p$m))
    }
  } else {
    ## legacy window of the original notebooks: closed on both ends
    day <- t %% 365
    if (day < p$control_start || day > p$control_start + p$control_duration) {
      return(rep(0, p$m))
    }
  }
  p$control_power
}

## Right-hand side of the ODE system -----------------------------------------
## State vector: S_1..S_m, I_1..I_m, R_1..R_m, Ntot (integrated total
## abundance, kept as a state variable so that prevalence is computed on the
## same quantity as in the original analysis).
## Transmission is frequency-dependent on the *whole community*:
##   lambda_i = tau_i * sum_j C[i, j] * I_j / Ntot
## Demography is logistic (Verhulst): births are damped by (1 - N_i / K_i) and
## deaths increase with N_i / K_i. The recovered class follows the same
## density-dependent mortality as S and I; it stays empty here since sigma = 0.
seov_rhs <- function(t, state, p) {
  m <- p$m
  S <- state[1:m]
  I <- state[m + (1:m)]
  R <- state[2 * m + (1:m)]

  N    <- S + I + R
  Ntot <- sum(N)

  b    <- seov_birth_rate(t, p)
  cull <- seov_control_rate(t, p)

  lambda <- p$tau * as.vector(p$C %*% (I / Ntot))

  dS <- b * (1 - N / p$K) * N - (lambda + p$d * N / p$K) * S - cull * S
  dI <- lambda * S - (p$sigma + p$d * N / p$K + p$alpha) * I - cull * I
  dR <- p$sigma * I - (p$d * N / p$K) * R - cull * R

  list(c(dS, dI, dR, sum(dS + dI + dR)))
}

## Scenario helpers ----------------------------------------------------------

## Carrying capacities for a given commensal:wild abundance ratio.
## Commensal species get K_ref * sqrt(ratio), wild species K_ref / sqrt(ratio),
## so that the total carrying capacity of the community stays comparable.
seov_carrying_capacity <- function(ratio, p) {
  K <- rep(p$K_ref, p$m)
  commensal <- p$species$guild == "commensal"
  K[commensal]  <- p$K_ref * sqrt(ratio)
  K[!commensal] <- p$K_ref / sqrt(ratio)
  K
}

## Initial abundances for a given ratio.
## Ntot0 individuals are split between species; Rattus norvegicus and
## Mus musculus (the two species observed in the park) carry the `ratio`
## weight (p$focal_species), every other species receives the baseline weight
## 1. The realised
## commensal:wild abundance ratio at t = 0 is therefore (2 * ratio + 1) / 3,
## the nominal ratio being reached only once the community has converged
## towards its carrying capacities.
seov_initial_abundance <- function(ratio, p) {
  w <- rep(1, p$m)
  w[p$species$latin %in% p$focal_species] <- ratio
  p$Ntot0 * w / sum(w)
}

## Per-species culling rate for a control strategy.
##   "none"      : no control
##   "competent" : targeted control (trapping) of the two competent hosts of
##                 the virus -- the two Rattus species for SEOV in the Lyon
##                 community, the two Peromyscus species for SNV in the Mexican
##                 one
##   "all"       : non-targeted control (chemical) -- every small mammal species
seov_control_power <- function(power, target, p) {
  switch(target,
    none      = rep(0, p$m),
    competent = ifelse(p$species$is_competent, power, 0),
    all       = rep(power, p$m),
    stop("unknown control target: ", target)
  )
}

## Season / management phase of each time step, used to colour phase diagrams.
seov_phase <- function(t, p, controlled) {
  dy    <- seov_day_of_year(t)
  phase <- ifelse(dy <= p$breeding_days, "Peaks of reproduction", "Winter")
  if (controlled) {
    in_control <- dy >= p$control_start & dy < p$control_start + p$control_duration
    phase[in_control] <- "Deratization"
  }
  factor(phase, levels = c("Deratization", "Peaks of reproduction", "Winter"))
}
