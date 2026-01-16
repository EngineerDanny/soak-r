library(data.table)
library(mlr3)
library(mlr3resampling)
library(ggplot2)
library(grid)

compute_pooling_gain_summary <- function(bmr, task, grouping, label_col,
                                         algorithms = c("cv_glmnet", "ranger"),
                                         train_groups = c("other", "all"),
                                         use_log10 = FALSE) {
  score.dt <- as.data.table(mlr3resampling::score(bmr, mlr3::msr("regr.mse")))
  
  if ("test.subset" %in% names(score.dt))   setnames(score.dt, "test.subset",  "test.group")
  if ("train.subsets" %in% names(score.dt)) setnames(score.dt, "train.subsets","train.groups")
  if (!"test.fold" %in% names(score.dt) && "iteration" %in% names(score.dt)) {
    setnames(score.dt, "iteration", "test.fold")
  }
  
  required_cols <- c("regr.mse","algorithm","test.group","train.groups","test.fold")
  miss <- setdiff(required_cols, names(score.dt))
  if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse=", "))
  
  score.dt[algorithm == "regr.rpart", algorithm := "rpart"]
  score.dt[algorithm == "regr.ranger", algorithm := "ranger"]
  score.dt[algorithm %in% c("regr.kknn", "kknn"), algorithm := "knn"]
  
  score.dt <- score.dt[
    algorithm %in% algorithms & train.groups %in% c("same", train_groups)
  ]
  
  if (use_log10) {
    score.dt[, mse_value := log10(regr.mse)]
  } else {
    score.dt[, mse_value := regr.mse]
  }
  
  wide.dt <- dcast(
    score.dt,
    algorithm + test.group + test.fold ~ train.groups,
    value.var = "mse_value"
  )
  if (!("same" %in% names(wide.dt))) stop("Need 'same' results for paired deltas.")
  
  delta.list <- lapply(train_groups, function(group_name) {
    if (!(group_name %in% names(wide.dt))) return(NULL)
    wide.dt[!is.na(same) & !is.na(get(group_name)), .(
      algorithm,
      test.group,
      train.groups = group_name,
      test.fold,
      delta = get(group_name) - same
    )]
  })
  delta.dt <- rbindlist(delta.list, use.names = TRUE)
  if (nrow(delta.dt) == 0) stop("No paired deltas found for train groups.")
  
  summary.dt <- delta.dt[, .(
    delta_mean = mean(delta),
    delta_sd = sd(delta),
    nfold = .N,
    p_value = if (.N > 1) t.test(delta, mu = 0)$p.value else NA_real_
  ), by = .(algorithm, test.group, train.groups)]
  summary.dt[, delta_se := delta_sd / sqrt(nfold)]
  summary.dt[, delta_ci_low := ifelse(nfold > 1,
                                      delta_mean + qt(0.025, nfold - 1) * delta_se,
                                      NA_real_)]
  summary.dt[, delta_ci_high := ifelse(nfold > 1,
                                       delta_mean + qt(0.975, nfold - 1) * delta_se,
                                       NA_real_)]
  
  summary.dt[, `:=`(
    task = task,
    grouping = grouping
  )]
  
  if (!is.null(label_col) && label_col %in% names(score.dt)) {
    label.map <- score.dt[, .(label = get(label_col)[1L]), by = test.group]
    summary.dt[label.map, test_label := i.label, on = "test.group"]
  }
  summary.dt[is.na(test_label), test_label := as.character(test.group)]
  
  summary.dt
}

plot_pooling_gain_heatmap <- function(summary.dt, grouping,
                                      algorithms = c("cv_glmnet", "ranger"),
                                      train_groups = c("other", "all"),
                                      use_log10 = FALSE,
                                      title = "",
                                      outfile = "pooling-gain-heatmap.png") {
  dt <- summary.dt[grouping == grouping]
  if (nrow(dt) == 0) stop("No rows for grouping: ", grouping)
  
  order.dt <- dt[algorithm == "cv_glmnet" & train.groups == "all",
                 .(order_val = mean(delta_mean)), by = test_label][order(order_val)]
  if (nrow(order.dt) > 0) {
    dt[, test_label := factor(test_label, levels = order.dt$test_label)]
  } else {
    dt[, test_label := factor(test_label, levels = unique(test_label))]
  }
  
  dt[, `:=`(
    task = factor(task, levels = unique(task)),
    algorithm = factor(algorithm, levels = algorithms),
    train.groups = factor(train.groups, levels = train_groups)
  )]
  
  delta_label <- if (use_log10) {
    "Delta log10(MSE):\n(train - same)\nBlue (-) = benefit;\nOrange (+) = harm"
  } else {
    "Delta MSE (train - same)\nBlue (-) = benefit;\nOrange (+) = harm"
  }
  
  dt[, delta_label_text := sprintf("%.2f", ifelse(abs(delta_mean) < 0.005, 0, delta_mean))]
  
  x_label <- if (grouping == "treatment") "Treatment" else "Site"
  
  gg <- ggplot(dt, aes(x = test_label, y = task, fill = delta_mean)) +
    geom_tile(color = "grey80", linewidth = 0.3) +
    geom_text(aes(label = delta_label_text), size = 2.2, color = "black") +
    facet_grid(algorithm ~ train.groups) +
    scale_fill_gradient2(
      low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
      name = delta_label
    ) +
    labs(title = title, x = x_label, y = "Task / experiment") +
    theme_bw() +
    theme(
      panel.spacing = grid::unit(0.4, "lines"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 9)
    )
  
  ggsave(outfile, gg, height = 6, width = 10, units = "in", dpi = 300)
  message("Wrote: ", outfile)
  invisible(gg)
}

## ---- Example runs (skips missing files) ----
runs <- list(
  list(
    file = "qsip_pc2_all_new-dim.C.between.sites.RData",
    task = "DIM + C",
    grouping = "site",
    label_col = "site"
  ),
  list(
    file = "qsip_pc2_all_new-dim.CN.combine.sites.RData",
    task = "DIM + CN",
    grouping = "site",
    label_col = "site"
  ),
  list(
    file = "qsip_pc2_all_new-dim.compare.sites.RData",
    task = "DIM all treatments",
    grouping = "site",
    label_col = "site"
  )
)

summary.list <- list()
for (run in runs) {
  if (!file.exists(run$file)) {
    message("Skip missing: ", run$file)
    next
  }
  load(run$file)
  summary.list[[run$file]] <- compute_pooling_gain_summary(
    bmr,
    task = run$task,
    grouping = run$grouping,
    label_col = run$label_col,
    algorithms = c("cv_glmnet", "ranger"),
    train_groups = c("other", "all"),
    use_log10 = T
  )
}

summary.dt <- rbindlist(summary.list, use.names = TRUE, fill = TRUE)
if (nrow(summary.dt) > 0) {
  plot_pooling_gain_heatmap(
    summary.dt,
    grouping = "site",
    algorithms = c("cv_glmnet", "ranger"),
    train_groups = c("other", "all"),
    use_log10 = T,
    title = "Pooling gain heatmap (site)",
    outfile = "qsip_pooling_gain_heatmap_site.png"
  )
}
