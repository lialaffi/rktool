# rktool

Analysis of **rhythmic ratios (rk)** and **rhythmic categories** from any
sequence of events (onsets or inter-onset intervals).

Terminology: a **category** (or peak) is a ratio such as `"1_2"` (1:2),
`"1_1"` or `"2_1"`; each category has two **bins**, `"on"` (central part)
and `"off"` (outer parts), e.g. `"1_2_on"` / `"1_2_off"`.

**Package:** `rktool`  
**Author / maintainer:** Lia Laffi

## Installation

```r
# from a local folder
install.packages("path/to/rktool", repos = NULL, type = "source")
# or, during development
devtools::load_all()

# with the vignette (a step-by-step guide to the whole workflow)
devtools::install("path/to/rktool", build_vignettes = TRUE)
vignette("rktool")
```

## Workflow

The package includes an example dataset, `rk_example` (event onsets of
17 individuals in 202 recordings, see `?rk_example`), so the whole workflow can
be run as is:

```r
library(rktool)

# 1. intervals tk from onsets (skip if tk are already computed), then rk
y <- calculate_tk(rk_example, onset_col = onset, by = c(file_id, ind_id))
y <- calculate_rk(y, tk_col = tk, by = c(file_id, ind_id))

# 2. classification into bins (1_5_off ... 1_1_on ... 5_1_off)
y$bins <- classify_rk(y$rk)

# 3. bin counts per group (wide or long)
b <- count_rk(y, by = c(file_id, ind_id, sex), bin_col = bins)

# 4. normalization: null distribution or bin width
norm_df <- compute_null_distributions(y, tk_col = tk, type_col = sex,
                                      distribution = "lognormal")
b_norm  <- normalize_rk_counts(b, norm_df, type_col = sex)
b_width <- normalize_rk_counts(b, method = "width")

# 5. on vs off test repeated over multiple null simulations
runs <- simulate_null_runs(y, tk_col = tk, type_col = sex,
                           distribution = "lognormal", runs = 100)
res  <- test_peaks(b, runs, type_col = sex)
res$summary

# density peaks of rk
find_rk_peaks(y$rk, peaks = c("1_1", "1_2", "2_1"))
```

## Functions

| Function | Purpose |
|---|---|
| `calculate_tk()` | intervals tk = onset - previous onset, by group |
| `calculate_rk()` | rk = tk / (tk + tk+1) from tk, by group |
| `classify_rk()` | classifies rk into the on/off bins of each category |
| `count_rk()` | counts bins per group (wide/long) |
| `compute_null_distributions()` | expected bin proportions under a null distribution |
| `simulate_null_runs()` | repeats the null simulation multiple times |
| `normalize_rk_counts()` | normalizes bin counts (null distribution or bin width) |
| `test_peaks()` | paired Wilcoxon on > off for each category, across multiple normalizations |
| `find_rk_peaks()` | mode of the rk density for each category |
| `plot_rk_density()` | density of tk and/or rk, with facets and category reference lines |
| `plot_rk_ternary()` | ternary plot of triplets of consecutive tk (requires ggtern) |
| `rk_example` | example dataset of event onsets |
