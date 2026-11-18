## ---------------------------------------------------------------------------
## 00 -- Life-history traits of the two host communities
##
## Reads data/communities.csv (which species belong to which community) and
## data/pantheria_subset.csv (PanTHERIA, Jones et al. 2009), and writes
## data/species_traits.csv, the file every other script reads.
##
## Rule for a missing trait, applied in this order and recorded in the
## `*_source` columns of the output:
##   1. value of the species in PanTHERIA                 -> "PanTHERIA"
##   2. value declared in data/trait_overrides.csv        -> the source given there
##   3. mean of the genus in PanTHERIA                    -> "PanTHERIA genus mean (n = k)"
## Nothing is filled in by hand outside data/trait_overrides.csv, and
## gestation lengths are rounded to the day, as in the original analysis.
##
## Run from the root of the project: Rscript code/scripts/00_species_traits.R
## ---------------------------------------------------------------------------

source("code/R/paths.R")

communities <- read.csv("data/communities.csv", stringsAsFactors = FALSE)
pantheria   <- read.csv("data/pantheria_subset.csv", stringsAsFactors = FALSE)
overrides   <- read.csv("data/trait_overrides.csv", stringsAsFactors = FALSE)

TRAITS <- c("litter_size", "litters_per_year", "gestation_days")

lookup <- function(binomial, trait) {
  sp <- pantheria[pantheria$binomial == binomial, ]
  if (nrow(sp) == 1 && !is.na(sp[[trait]])) {
    return(list(value = sp[[trait]], source = "PanTHERIA"))
  }
  ov <- overrides[overrides$binomial == binomial & overrides$trait == trait, ]
  if (nrow(ov) == 1) return(list(value = ov$value, source = ov$source))
  genus <- if (nrow(sp) == 1) sp$genus else strsplit(binomial, " ")[[1]][1]
  pool  <- pantheria[[trait]][pantheria$genus == genus]
  pool  <- pool[!is.na(pool)]
  if (length(pool) > 0) {
    return(list(value = mean(pool),
                source = sprintf("PanTHERIA genus mean (n = %d)", length(pool))))
  }
  stop("no value and no override for ", binomial, " / ", trait)
}

traits <- communities
for (trait in TRAITS) {
  got <- lapply(communities$latin, lookup, trait = trait)
  value <- vapply(got, function(x) as.numeric(x$value), numeric(1))
  if (trait == "gestation_days") value <- round(value)
  traits[[trait]] <- round(value, 2)
  traits[[paste0(trait, "_source")]] <- vapply(got, function(x) x$source, character(1))
}

## the order of the species inside a community is the order of the columns of
## the contact matrices, and is fixed here once and for all
traits <- traits[order(traits$community != "France", seq_len(nrow(traits))), ]

write.csv(traits, "data/species_traits.csv", row.names = FALSE)
print(traits[, c("community", "virus", "code", "latin", "guild", "competence", TRAITS)])
message("\n-> data/species_traits.csv")
