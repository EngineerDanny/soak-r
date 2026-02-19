local.mlr3resampling <- "/projects/genomic-ml/da2343/mlr3resampling"
if (dir.exists(local.mlr3resampling) && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(local.mlr3resampling, quiet = TRUE)
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mlr3)
  library(mlr3resampling)
})

rdata.file <- "qsip_pc2_all_new-dim.control_vs_qme.15.downsample.sizes0.compare.subsets.RData"
rdata.stem <- sub("\\.RData$", "", basename(rdata.file))
pvalue.plot.dir <- paste0(rdata.stem, ".downsample-pvalue-plots")

save_plot <- function(path, plot_obj, width, height) {
  ggsave(path, plot_obj, width = width, height = height, dpi = 400, bg = "white")
}

sanitize_file_part <- function(x) {
  gsub("[^A-Za-z0-9._-]", "_", x)
}

if (!exists("pvalue_downsample", where = asNamespace("mlr3resampling"), inherits = FALSE)) {
  stop("mlr3resampling::pvalue_downsample() is not available.")
}

if (!file.exists(rdata.file)) {
  stop(sprintf("Missing RData file: %s", rdata.file))
}

load(rdata.file)

if (!exists("bmr")) {
  stop("Expected object `bmr` in RData, but it was not found.")
}

soak_sizes_score <- as.data.table(mlr3resampling::score(bmr, mlr3::msr("regr.rmse")))

required.cols <- c("test.subset", "algorithm")
missing.cols <- setdiff(required.cols, names(soak_sizes_score))
if (length(missing.cols) > 0L) {
  stop(sprintf("Missing required columns in score data: %s", paste(missing.cols, collapse = ", ")))
}

combo.dt <- unique(soak_sizes_score[, .(
  test.subset = as.character(test.subset),
  algorithm = as.character(algorithm)
)])

if (nrow(combo.dt) == 0L) {
  stop("No subset/algorithm combinations found in score data.")
}

dir.create(pvalue.plot.dir, showWarnings = FALSE, recursive = TRUE)
saved.pvalue.plots <- character()

for (combo.i in seq_len(nrow(combo.dt))) {
  subset_name <- combo.dt$test.subset[[combo.i]]
  model_name <- combo.dt$algorithm[[combo.i]]

  down.obj <- tryCatch(
    mlr3resampling::pvalue_downsample(soak_sizes_score, subset_name, model_name),
    error = function(e) {
      message(sprintf(
        "Skipping pvalue_downsample for subset=%s, model=%s: %s",
        subset_name, model_name, conditionMessage(e)
      ))
      NULL
    }
  )
  if (is.null(down.obj)) {
    next
  }

  down.plot <- plot(down.obj)
  out.pvalue <- file.path(
    pvalue.plot.dir,
    sprintf(
      "pvalue-%s-%s.png",
      sanitize_file_part(subset_name),
      sanitize_file_part(model_name)
    )
  )
  save_plot(out.pvalue, down.plot, width = 10.4, height = 3.8)
  saved.pvalue.plots <- c(saved.pvalue.plots, out.pvalue)
}

if (length(saved.pvalue.plots) == 0L) {
  message("No pvalue_downsample plots were created.")
} else {
  for (plot.path in saved.pvalue.plots) {
    cat(sprintf("Saved plot: %s\\n", plot.path))
  }
}
