## ---------------------------------------------------------------------------
## Working directory and paths
##
## Every path in this project is relative to the root of the project, the
## directory that contains data/, code/ and figures/. Nothing else is needed for
## the code to run on any machine: no absolute path, no setwd() inside a
## script, no environment variable.
##
##   cd <project root>
##   Rscript run_all.R                 # everything
##   Rscript code/scripts/fig3_dilution.R   # one figure
##
## In RStudio, open the .Rproj file at the root, or use Session > Set Working
## Directory > To Project Directory.
## ---------------------------------------------------------------------------

if (!all(dir.exists(c("data", "code/R", "code/scripts")))) {
  stop("This script must be run from the root of the project, the directory ",
       "that contains data/, code/ and figures/. Current directory: ", getwd())
}

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
