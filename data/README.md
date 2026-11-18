# Data

Every input of the model. Nothing else is needed to regenerate the figures and
the tables of the paper; everything in `results/` and `figures/` is derived
from these files by the scripts.

## `communities.csv`

The two host communities. One row per species and per community.

| column | content |
|---|---|
| `community` | `France` (urban park of Lyon) or `Mexico` (peri-urban northeastern Mexico) |
| `virus` | the hantavirus of that community: `SEOV` or `SNV` |
| `code` | short name used in the code and in the figures |
| `latin` | species name |
| `guild` | `commensal` or `wild` |
| `competence` | `primary`, `secondary` or `none`: the two competent hosts of the virus are the two *Rattus* species for SEOV and the two *Peromyscus* species for SNV |
| `notes` | why the species is in the community, with its source |

Sources: Hussein et al. (2023) and Bouilloud et al. (in preparation) for the
European community; Vigueras-Galván et al. (2019) for the hantaviruses of
Mexico, their host species and the states where SNV has been detected.

## `pantheria_subset.csv`

Life-history traits of the genera used here, extracted from PanTHERIA
(Jones, K. E. et al. 2009, *Ecology* 90:2648,
https://doi.org/10.1890/08-1494.1), downloaded from
https://esapubs.org/archive/ecol/E090/184/PanTHERIA_1-0_WR05_Aug2008.txt.
Empty cells are the `-999` of the original file, i.e. missing values. Columns:
`binomial`, `genus`, `family`, `gestation_days`, `litter_size`,
`litters_per_year`, `body_mass_g`.

## `trait_overrides.csv`

Trait values taken from another source than PanTHERIA, one row each, with the
justification in the `source` column. It currently holds a single value, the
number of litters per year of *Crocidura russula*, which is documented for no
species of the genus in PanTHERIA.

## `species_traits.csv` (generated)

Written by `code/scripts/00_species_traits.R` from the three files above: the traits
actually used by the model, with, for each of them, a `*_source` column saying
whether the value is the species value of PanTHERIA, the mean of its genus in
PanTHERIA, or an override. This is the file `R/parameters.R` reads. It is
reproduced as Table S1.

## `contact_matrix_france.csv`, `contact_matrix_mexico.csv`

The contact matrices, one row and one column per species, in the order of the
`species` column. They are symmetric, their diagonal is 1, and their
off-diagonal values are 0.1, 0.01 or 0.001. These coefficients are assumptions
derived from the habitats and the behaviour of the species, not measurements.

## `contact_matrix_rationale.csv`

The justification of every level of the two matrices: which pairs receive the
high, moderate and low levels, and why.

## `TrappingSuccess_building.csv`, `TrappingSuccess_forest.csv`, `TrappingSuccess_Rats.csv`

Trapping success recorded in buildings, in forest and for rats, as collected.
They document the order of magnitude of the initial community
(`N_tot(0) = 120`) and are described in Supplementary Material S1. No parameter
of the model is estimated from them.
