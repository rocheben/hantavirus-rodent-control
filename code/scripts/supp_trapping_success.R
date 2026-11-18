## ---------------------------------------------------------------------------
## Supplementary Material -- trapping tables behind the initial abundances.
##
## Reads the trapping tables delivered in data/ and writes the summary table
## quoted in the Supplementary Material. Nothing is recomputed from them: they
## document the order of magnitude of the initial community.
## ---------------------------------------------------------------------------

source("code/R/paths.R")
suppressPackageStartupMessages({ library(dplyr) })

files <- c(building = "data/TrappingSuccess_building.csv",
           forest   = "data/TrappingSuccess_forest.csv",
           rats     = "data/TrappingSuccess_Rats.csv")

shape <- lapply(names(files), function(habitat) {
  path <- files[[habitat]]
  if (!file.exists(path)) return(NULL)
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                fileEncoding = "UTF-8-BOM")
  if (ncol(d) == 1) d <- read.csv(path, sep = ";", stringsAsFactors = FALSE,
                                  check.names = FALSE, fileEncoding = "UTF-8-BOM")
  data.frame(habitat = habitat, file = basename(path), rows = nrow(d),
             columns = paste(names(d), collapse = "; "), stringsAsFactors = FALSE)
})
shape <- bind_rows(shape)
write.csv(shape, "results/S1_trapping_tables.csv", row.names = FALSE)

print(shape)
