# Rodent management, biodiversity and hantavirus circulation

Data, code, results and figures for

> Bouilloud M., Canet A., Vigueras-Galván A.L., Arnal A., Charbonnel N. &
> Roche B. *Who is the reservoir? Modelling win-win rodent control for
> hantavirus hazard and small mammal biodiversity in France and Mexico.*

The model follows the circulation of a hantavirus in a community of six urban
small mammals and is used to (i) relate prevalence to small mammal diversity
and (ii) compare targeted and non-targeted rodent control. It is applied to two
systems: the Lyon community of an urban park of Lyon, where Seoul
orthohantavirus (SEOV) circulates in the two *Rattus* species, and a peri-urban
community of northeastern Mexico, where Sin Nombre orthohantavirus (SNV)
circulates in two *Peromyscus* species. In the second system the competent
hosts are wild species rather than commensal ones, which is what the comparison
is about.

## Contents

```
data/                        every input of the model (see data/README.md)
code/R/                      the model: parameters, equations, integration, shared plots
code/scripts/                one script per figure, plus the traits and the simulations
figures/                     Figs. 1-5 and Figs. S1-S8, PNG at 300 dpi
results/                     cached simulations (.rds) and every table quoted in the text (.csv)
run_all.R                    regenerates everything, in order
```

## Reproducing the results

```bash
cd <this directory>
Rscript run_all.R                 # 15-30 min, regenerates figures/ and results/
```

Every path used by the code is relative to this directory, so nothing has to be
edited to run it on another machine. There is no `setwd()` in any script and no
absolute path anywhere; the scripts stop with an explicit message if they are
started from elsewhere. In RStudio, open `seov-rodent-control.Rproj` (or
Session > Set Working Directory > To Project Directory) and run
`source("run_all.R")`.

Required R packages: `deSolve`, `vegan`, `dplyr`, `tidyr`, `ggplot2`,
`patchwork`.

```r
install.packages(c("deSolve", "vegan", "dplyr", "tidyr", "ggplot2", "patchwork"))
```

The figures and tables in this repository were produced with R 4.3.3. A full
re-run with R 4.6.1 (deSolve 1.42, ggplot2 4) gives the same results: every
number quoted in the paper is unchanged, the tables differ by less than 0.5 %
(ODE solver round-off), and the figures differ only in
rendering details such as legend order and colour shades.

## One script per figure

`code/scripts/01_simulate.R` runs every simulation once and caches it in `results/`;
each figure script then reads a cached simulation and writes exactly one
figure, so a single panel can be redrawn in a few seconds:

| script | output |
|---|---|
| `code/scripts/00_species_traits.R` | `data/species_traits.csv` (traits of the two communities) |
| `code/scripts/01_simulate.R` | `results/sim_*.rds` (every simulation of both communities, 15-30 min) |
| `code/scripts/fig1_model_diagram.R` | Fig. 1, model diagram |
| `code/scripts/fig2_communities.R` | Fig. 2, the two communities: species, guilds, competence, contacts |
| `code/scripts/fig3_dilution.R` | Fig. 3, prevalence vs diversity, both communities |
| `code/scripts/fig4_control.R` | Fig. 4, control strategies, both communities |
| `code/scripts/fig5_community.R` | Fig. 5, community dynamics, both communities |
| `code/scripts/figS1_model_inputs.R` | Fig. S1, birth rates and control calendars |
| `code/scripts/figS2_phase_diagrams_france.R`, `code/scripts/figS3_phase_diagrams_mexico.R` | Figs. S2, S3, phase diagrams |
| `code/scripts/figS4_control_timing.R` | Fig. S4, timing of the campaign |
| `code/scripts/figS5_control_calendars_france.R`, `code/scripts/figS6_control_calendars_mexico.R` | Figs. S5, S6, the four strategy x calendar cases |
| `code/scripts/figS7_community_spring.R`, `code/scripts/figS8_community_targeted.R` | Figs. S7, S8 |
| `code/scripts/table1_parameters.R` | Table 1 and Tables S1, S2 |
| `code/scripts/supp_trapping_success.R` | the summary of the trapping tables |

Every figure of the main text carries both communities as panels: the Lyon
community first, the Mexican one last.

## The model in one page

Each of the *m* = 6 species of a community is described by a susceptible and an
infected compartment. Demography is logistic and seasonal, transmission is
frequency-dependent on the whole community, and rodent control is an extra
mortality applied over a 30-day campaign each year.

```
dS_i/dt = b_i(t) (1 - N_i/K_i) N_i - (lambda_i + d_i N_i/K_i) S_i - rho_i(t) S_i
dI_i/dt = lambda_i S_i - (d_i N_i/K_i) I_i - rho_i(t) I_i
lambda_i = tau_i * sum_j C_ij I_j / N_tot
```

* `b_i(t)` — seasonal birth rate built from the life-history traits of
  `data/species_traits.csv`: 80 % of the births fall in the first half of the
  year, as pulses of length `Gesta_i` repeated every `floor(182 / LPY_i)` days,
  the remaining 20 % are spread over the second half. The pulses are anchored
  on the day of the year, so the reproductive calendar is the same in every
  simulated year (Fig. S1).
* `d_i = LS_i × LPY_i / (2 × 365)` — intrinsic death rate, calibrated on the
  mean intrinsic birth rate; the realised ratio ranges from 0.97 to 1.12
  (`results/birth_death_balance.csv`).
* `C` — symmetric contact matrix read from
  `data/contact_matrix_<community>.csv`: 1 within species, 0.1 / 0.01 / 0.001
  between species depending on habitat sharing and agonistic interactions.
* `tau_i` — susceptibility to the virus: 0.13 for the primary reservoir
  (*R. norvegicus* for SEOV, *P. maniculatus* for SNV), 0.104 for the secondary
  one (*R. rattus*, *P. leucopus*), 0 for every other species.
* `K_i` — carrying capacity, `100 × sqrt(ratio)` for commensal species and
  `100 / sqrt(ratio)` for wild species, where `ratio` is the relative abundance
  of commensal species that defines the epidemiological context. Initial
  abundances give the weight `ratio` to *R. norvegicus* and *M. musculus* only,
  so the realised ratio at t = 0 is `(2 × ratio + 1) / 3`.
* `rho_i(t)` — control intensity, applied from day `control_start` for 30 days
  every year, and only after a 10-year burn-in. Targeted control
  ("competent") affects the two competent hosts, non-targeted control ("all")
  every species of the community.

The recovery rate and the disease-induced mortality are zero, so the SIR
skeleton coded in `R/model.R` reduces to the SI model of the manuscript.

## Conventions

* Time starts at `t = 1`; the n-th year is made of the 365 time steps
  `(n-1) * 365 + 1 ... n * 365`, so year 30 is days 10 586 to 10 950.
* The day of the year runs from 1 to 365 (`seov_day_of_year()`).
* A campaign covers `[control_start, control_start + 30)`, exactly 30 days.
* "Infected individuals" always means every infected individual of the
  community (column `Itot`); prevalence among the competent hosts alone is
  `Ires / Nres`.

## Experiments

| experiment | scenarios | communities |
|---|---|---|
| dilution effect | ratio ∈ {0.5, 1, 5}, no control | both |
| control grid | ratio ∈ {0.8, 1, 10} × intensity ∈ {0, 0.021, …, 0.147, 0.16} × {targeted, non-targeted} × campaign start ∈ {day 1, day 181} | both |
| community dynamics | ratio = 1, intensity ∈ {0, 0.1, 0.2} × {targeted, non-targeted} × campaign start ∈ {day 1, day 181} | both |
| timing of control | ratio = 1, targeted at the competent hosts, start day ∈ {1, 91, 181} × intensity ∈ {0, 0.035, 0.07, 0.105} | both |

Simulations run for 40 years at a daily time step. Control starts after year
10; Figures 4, 5 and S2 to S8 describe year 30; Figure 3 uses the trajectory
after day 8000 (year 22 onwards) and the local extrema of the infected time
series over that window.
