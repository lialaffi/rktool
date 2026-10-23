# Builds data/data_rktool.rda from data-raw/rktool_dataset_example.csv.
# Run from the package root: source("data-raw/data_rktool.R")

data_rktool <- utils::read.csv("data-raw/rktool_dataset_example.csv",
                              stringsAsFactors = FALSE)

dir.create("data", showWarnings = FALSE)
save(data_rktool, file = "data/data_rktool.rda", compress = "xz")
