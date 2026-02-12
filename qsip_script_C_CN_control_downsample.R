library(data.table)
library(mlr3resampling)

qsip.dt <- fread("/projects/genomic-ml/necromass/data-qsip-jeff/qsip_pc2_all_new.csv")
setwd("/projects/genomic-ml/da2343/soak/")

gene.names <- grep("^K", names(qsip.dt), value = TRUE)

compare <- function(comparison, group.var, DT) {
  data.table(comparison, group.var, DT = list(DT))
}

keep.treatments <- c("C", "CN", "control")
dim.ccc <- qsip.dt[
  experiment == "dim" & treatment %chin% keep.treatments
]
comparison.dt <- rbind(
  compare("dim.C_CN_control.downsample.sizes0.compare.treatments", "treatment", dim.ccc)
)

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
  debug_n_per_treatment <- 200L
  task.dt <- data.table(compare.row$DT[[1]][, task.names, with = FALSE])[
    , .SD[sample(.N, min(.N, debug_n_per_treatment))], by = treatment
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
