## ---------------------------------------------------------------------------
## Running the model and shaping its output
## ---------------------------------------------------------------------------

## Integrate one scenario ----------------------------------------------------
##
## ratio          : initial and asymptotic abundance of commensal species
##                  relative to wild small mammals
## control_power  : culling rate added to the mortality of targeted species
## control_target : "none", "competent" (targeted at the competent hosts) or
##                  "all" (non-targeted)
## control_start  : first day of the year at which control is applied
##
## Returns one row per time step with, for every species, S / I / R / N, plus
## the community-level quantities Itot (all infected individuals), Nres and
## Ires (abundance and infected individuals of the two competent hosts), Ntot
## and the Shannon diversity index of the community.
seov_simulate <- function(ratio,
                          control_power  = 0,
                          control_target = "none",
                          control_start  = 1,
                          p) {

  p$K             <- seov_carrying_capacity(ratio, p)
  p$control_power <- seov_control_power(control_power, control_target, p)
  p$control_start <- control_start

  N0 <- seov_initial_abundance(ratio, p)

  ## one infected individual of the primary reservoir is introduced at t = 0
  seed <- which(p$species$competence == "primary")
  S0 <- N0; S0[seed] <- S0[seed] - 1
  I0 <- rep(0, p$m); I0[seed] <- 1
  R0 <- rep(0, p$m)

  state0 <- c(S0, I0, R0, sum(N0))
  names(state0) <- c(paste0("S_", 1:p$m), paste0("I_", 1:p$m),
                     paste0("R_", 1:p$m), "Ntot_state")

  out <- deSolve::ode(y = state0, times = p$times, func = seov_rhs, parms = p)
  out <- as.data.frame(out)

  S <- as.matrix(out[, paste0("S_", 1:p$m)])
  I <- as.matrix(out[, paste0("I_", 1:p$m)])
  R <- as.matrix(out[, paste0("R_", 1:p$m)])
  N <- S + I + R
  colnames(N) <- paste0("N_", 1:p$m)

  res_sp <- which(p$species$is_competent)

  res <- data.frame(
    time  = out$time,
    year  = seov_year_index(out$time),
    community = p$community,
    ratio = ratio,
    control_power  = control_power,
    control_target = control_target,
    control_start  = control_start,
    Itot    = rowSums(I),
    Ires    = rowSums(I[, res_sp, drop = FALSE]),
    Nres    = rowSums(N[, res_sp, drop = FALSE]),
    Ntot    = out$Ntot_state,
    Shannon = vegan::diversity(N),
    phase   = seov_phase(out$time, p, controlled = control_power > 0)
  )
  cbind(res, S, I, R, N)
}

## Long format --------------------------------------------------------------
## One row per (time, species); keeps the community-level columns.
seov_long <- function(sim, p) {
  idx <- seq_len(p$m)
  scen <- c("time", "year", "community", "ratio", "control_power", "control_target",
            "control_start", "Itot", "Ires", "Nres", "Ntot", "Shannon", "phase")
  long <- lapply(idx, function(i) {
    d <- sim[, scen, drop = FALSE]
    d$species    <- p$species$code[i]
    d$guild      <- p$species$guild[i]
    d$S          <- sim[[paste0("S_", i)]]
    d$I          <- sim[[paste0("I_", i)]]
    d$Abundance  <- sim[[paste0("N_", i)]]
    d
  })
  long <- do.call(rbind, long)
  long$species <- factor(long$species, levels = p$species$code)
  long
}

## Run a grid of scenarios --------------------------------------------------
## `grid` must be a data.frame whose columns are arguments of seov_simulate().
## `summarise_fun` is applied to each simulation before the results are bound
## together, which keeps memory use low for large grids.
seov_run_grid <- function(grid, p, summarise_fun = identity,
                          verbose = TRUE) {
  args <- setdiff(names(grid), "")
  res  <- vector("list", nrow(grid))
  for (k in seq_len(nrow(grid))) {
    call_args <- as.list(grid[k, args, drop = FALSE])
    if (verbose) {
      message(sprintf("[%d/%d] %s", k, nrow(grid),
                      paste(names(call_args), unlist(call_args),
                            sep = "=", collapse = ", ")))
    }
    sim <- do.call(seov_simulate, c(call_args, list(p = p)))
    res[[k]] <- summarise_fun(sim)
  }
  do.call(rbind, res)
}

## Local extrema ------------------------------------------------------------
## TRUE where x is the strict maximum of the window of `span` points centred
## on it (same convention as splus2R::peaks, which the original analysis used;
## points closer than span/2 to either end are never flagged).
seov_peaks <- function(x, span = 400) {
  n <- length(x)
  h <- span %/% 2
  out <- rep(FALSE, n)
  if (n <= 2 * h) return(out)
  for (i in (h + 1):(n - h)) {
    w <- x[(i - h):(i + h)]
    out[i] <- all(x[i] > w[-(h + 1)])
  }
  out
}

## TRUE at both local maxima and local minima of x.
seov_local_extrema <- function(x, span = 400) {
  seov_peaks(x, span) | seov_peaks(-x, span)
}

## Window helpers ----------------------------------------------------------
## Simulation time starts at t = 1, so the n-th simulated year is made of the
## 365 time steps (n - 1) * 365 + 1 ... n * 365. Control starts after a
## burn-in of 3650 days, i.e. from the 11th year onwards.
seov_year_index <- function(t) as.integer(ceiling(t / 365))

seov_year <- function(sim, year) {
  sim[sim$year == year, , drop = FALSE]
}

## Day of the year, 1 to 365, of each time step.
seov_day_of_year <- function(t) as.integer((t - 1) %% 365 + 1)

## Legacy analysis window of the original notebooks: days [year * 365,
## (year + 1) * 365], i.e. 366 time steps, the last one being the first day of
seov_year_legacy <- function(sim, year) {
  sim[sim$time >= year * 365 & sim$time <= (year + 1) * 365, , drop = FALSE]
}
