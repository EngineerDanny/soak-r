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
  compare("dim.C_CN_control.compare.treatments", "treatment", dim.ccc)
)

# Optional: minimal debug run (no batchtools submit).
if (T) {
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
  same_other_cv <- mlr3resampling::ResamplingSameOtherCV$new()
  same_other_cv$param_set$values$folds <- 2
  cvg <- mlr3learners::LearnerRegrCVGlmnet$new()
  cvg$param_set$values$nfolds <- 3
  debug.learners <- list(cvg, mlr3::LearnerRegrFeatureless$new())
  debug.grid <- mlr3::benchmark_grid(task, debug.learners, same_other_cv)
  debug.result <- mlr3::benchmark(debug.grid)
  print(debug.result$aggregate())
}

# Optional: run the same benchmark workflow as qsip_script.R for this one comparison.
if (F) {
  for (comparison.i in 1:nrow(comparison.dt)) {
    compare.row <- comparison.dt[comparison.i]
    task.names <- c(compare.row$group.var, "growth_per_day", gene.names)
    same_other_cv <- mlr3resampling::ResamplingSameOtherCV$new()
    same_other_cv$param_set$values$folds <- 10
    task.dt <- data.table(compare.row$DT[[1]][, task.names, with = FALSE])
    task.id <- paste0("qsip_pc2_all_new_", compare.row$comparison)
    task <- mlr3::TaskRegr$new(
      task.id, task.dt, target = "growth_per_day"
    )$set_col_roles(compare.row$group.var, c("subset", "stratum"))
    cvg <- mlr3learners::LearnerRegrCVGlmnet$new()
    cvg$param_set$values$nfolds <- 5
    reg.learner.list <- list(
      cvg,
      mlr3::LearnerRegrFeatureless$new()
    )
    reg.bench.grid <- mlr3::benchmark_grid(task, reg.learner.list, same_other_cv)
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
