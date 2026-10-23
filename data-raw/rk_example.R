# Builds data/rk_example.rda from data-raw/data_example_rktool.csv.
# Run from the package root: source("data-raw/rk_example.R")

rk_example <- utils::read.csv("data-raw/data_example_rktool.csv",
                              stringsAsFactors = FALSE)

dir.create("data", showWarnings = FALSE)
save(rk_example, file = "data/rk_example.rda", compress = "xz")
