library(data.table)
library(mlr3resampling)

qsip.dt <- fread("/projects/genomic-ml/necromass/data-qsip-jeff/qsip_pc2_all_new.csv")
setwd("/projects/genomic-ml/da2343/soak/")

gene.names <- grep("^K", names(qsip.dt), value = TRUE)

compare <- function(comparison, group.var, DT) {
  data.table(comparison, group.var, DT = list(DT))
}

# Edit this block to change/add subsets that should be compared in one combined run.
# Each entry can filter on any qsip.dt column, and values can be one string or a vector.
subset.config <- list(
  list(
    subset_id = "GL",
    filters = list(experiment = "dim", treatment = "CN", site = "GL")
  ),
  list(
    subset_id = "PJ",
    filters = list(experiment = "dim", treatment = "CN", site = "PJ")
  ),
  list(
    subset_id = "PP",
    filters = list(experiment = "dim", treatment = "CN", site = "PP")
  ),
  list(
    subset_id = "MC",
    filters = list(experiment = "dim", treatment = "CN", site = "MC")
  )
)
comparison.name <- "dim.CN.between.sites.downsample.sizes0.compare.subsets"
comparison.group.var <- "subset_id"

apply_filters <- function(DT, filters) {
  out <- DT
  for (col.name in names(filters)) {
    vals <- as.character(filters[[col.name]])
    out <- if (length(vals) == 1L) {
      out[get(col.name) == vals]
    } else {
      out[get(col.name) %chin% vals]
    }
  }
  out
}

build_subset <- function(cfg, DT) {
  filter.cols <- names(cfg$filters)
  missing.cols <- setdiff(filter.cols, names(DT))
  if (length(missing.cols) > 0L) {
    stop("Unknown filter columns: ", paste(missing.cols, collapse = ", "))
  }

  subset.dt <- apply_filters(DT, cfg$filters)
  if (nrow(subset.dt) == 0L) {
    warning("Subset has 0 rows for subset_id: ", cfg$subset_id)
  }
  subset.dt[, (comparison.group.var) := cfg$subset_id]
  subset.dt
}

combined.dt <- rbindlist(
  lapply(subset.config, build_subset, DT = qsip.dt),
  use.names = TRUE,
  fill = TRUE
)
if (uniqueN(combined.dt[[comparison.group.var]]) < 2L) {
  stop("Need at least 2 non-empty subsets in subset.config for a comparison.")
}

comparison.dt <- compare(comparison.name, comparison.group.var, combined.dt)

make_downsample_cv <- function(folds = 10L, sizes = 0L) {
  downsample_cv <- mlr3resampling::ResamplingSameOtherSizesCV$new()
  downsample_cv$param_set$values$folds <- folds
  downsample_cv$param_set$values$sizes <- sizes
  downsample_cv
}

# Optional: minimal debug run (no batchtools submit).
if (F) {
  compare.row <- comparison.dt[1]
  task.names <- c(compare.row$group.var, "growth_per_day", gene.names)
  debug_n_per_group <- 200L
  group.var <- compare.row$group.var
  task.dt <- data.table(compare.row$DT[[1]][, task.names, with = FALSE])[
    , .SD[sample(.N, min(.N, debug_n_per_group))], by = group.var
  ]
  task.id <- paste0("debug_", compare.row$comparison)
  task <- mlr3::TaskRegr$new(
    task.id, task.dt, target = "growth_per_day"
  )$set_col_roles(compare.row$group.var, c("subset", "stratum"))

  downsample_cv <- make_downsample_cv(folds = 2L, sizes = 0L)
  cvg <- mlr3learners::LearnerRegrCVGlmnet$new()
  cvg$param_set$values$nfolds <- 3
  debug.learners <- list(cvg, mlr3::LearnerRegrFeatureless$new())

  debug.grid <- mlr3::benchmark_grid(task, debug.learners, downsample_cv)
  debug.result <- mlr3::benchmark(debug.grid)
  print(debug.result$aggregate())

  debug.score <- as.data.table(mlr3resampling::score(debug.result, mlr3::msr("regr.mse")))
  if ("train" %in% names(debug.score)) {
    debug.score[, n.train := vapply(train, length, integer(1))]
  }
  if ("n.train.groups" %in% names(debug.score)) {
    print(debug.score[, .(
      min_n_train = min(n.train.groups),
      max_n_train = max(n.train.groups),
      n_sizes = uniqueN(n.train.groups)
    ), by = train.subsets])
  }
}

# Optional: full batchtools workflow for downsampling analysis.
if (T) {
  for (comparison.i in 1:nrow(comparison.dt)) {
    compare.row <- comparison.dt[comparison.i]
    task.names <- c(compare.row$group.var, "growth_per_day", gene.names)
    task.dt <- data.table(compare.row$DT[[1]][, task.names, with = FALSE])
    task.id <- paste0("qsip_pc2_all_new_", compare.row$comparison)
    task <- mlr3::TaskRegr$new(
      task.id, task.dt, target = "growth_per_day"
    )$set_col_roles(compare.row$group.var, c("subset", "stratum"))

    downsample_cv <- make_downsample_cv(folds = 10L, sizes = 0L)
    cvg <- mlr3learners::LearnerRegrCVGlmnet$new()
    cvg$param_set$values$nfolds <- 5
    reg.learner.list <- list(
      cvg,
      mlr3::LearnerRegrFeatureless$new()
    )

    reg.bench.grid <- mlr3::benchmark_grid(task, reg.learner.list, downsample_cv)
    reg.dir <- paste0("qsip_pc2_all_new-", compare.row$comparison)
    unlink(reg.dir, recursive = TRUE)
    reg <- batchtools::makeExperimentRegistry(
      file.dir = reg.dir,
      seed = 1,
      packages = "mlr3verse"
    )
    mlr3batchmark::batchmark(reg.bench.grid, store_models = FALSE, reg = reg)
    job.table <- batchtools::getJobTable(reg = reg)
    chunks <- data.frame(job.table, chunk = 1)
    batchtools::submitJobs(chunks, resources = list(
      walltime = 60 * 60,
      memory = 32000,
      ncpus = 1,
      ntasks = 1,
      chunks.as.arrayjobs = TRUE
    ), reg = reg)
  }
}
