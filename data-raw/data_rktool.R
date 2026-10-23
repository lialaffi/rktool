# Builds data/data_rktool.rda from data-raw/data_example_rktool.csv.
# Run from the package root: source("data-raw/data_rktool.R")

data_rktool <- utils::read.csv("data-raw/data_example_rktool.csv",
                              stringsAsFactors = FALSE)

dir.create("data", showWarnings = FALSE)
save(data_rktool, file = "data/data_rktool.rda", compress = "xz")
