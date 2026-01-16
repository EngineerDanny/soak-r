library(data.table)
library(mlr3)
library(mlr3resampling)
library(ggplot2)

plot_pooling_gain <- function(bmr, comparison.RData, tit,
                              algorithms = c("cv_glmnet", "ranger"),
                              train_groups = c("other", "all"),
                              use_log10 = FALSE,
                              label_col = "site",
                              test_group_map = NULL,
                              colors = c(other = "#56B4E9", all = "#D55E00")) {
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
      delta = get(group_name) - same
    )]
  })
  delta.dt <- rbindlist(delta.list, use.names = TRUE)
  if (nrow(delta.dt) == 0) stop("No paired deltas found for train groups.")
  
  summary.dt <- delta.dt[, .(
    delta_mean = mean(delta),
    delta_sd = sd(delta),
    nfold = .N
  ), by = .(algorithm, test.group, train.groups)]
  summary.dt[, delta_se := delta_sd / sqrt(nfold)]
  summary.dt[, delta_lower := ifelse(nfold > 1,
                                     delta_mean + qt(0.025, nfold - 1) * delta_se,
                                     NA_real_)]
  summary.dt[, delta_upper := ifelse(nfold > 1,
                                     delta_mean + qt(0.975, nfold - 1) * delta_se,
                                     NA_real_)]
  
  label.map <- NULL
  if (!is.null(test_group_map)) {
    if (is.character(test_group_map) && !is.null(names(test_group_map))) {
      label.map <- data.table(
        test.group = names(test_group_map),
        label = unname(test_group_map)
      )
    } else if (is.data.table(test_group_map) &&
               all(c("test.group", "label") %in% names(test_group_map))) {
      label.map <- unique(test_group_map[, .(test.group, label)])
    } else {
      stop("test_group_map must be a named vector or data.table with test.group,label")
    }
  } else if (!is.null(label_col) && label_col %in% names(score.dt)) {
    label.map <- score.dt[, .(label = get(label_col)[1L]), by = test.group]
  }
  
  if (!is.null(label.map)) {
    summary.dt[label.map, label := i.label, on = "test.group"]
  }
  summary.dt[is.na(label), label := as.character(test.group)]
  summary.dt[, test.group.label := factor(label, levels = unique(label))]
  
  summary.dt[, train.groups := factor(train.groups, levels = train_groups)]
  summary.dt[, algorithm := factor(algorithm, levels = algorithms)]
  
  delta_label <- if (use_log10) {
    "Delta log10(MSE) (other/all - same)"
  } else {
    "Delta MSE (other/all - same)"
  }
  
  gg <- ggplot(summary.dt, aes(x = delta_mean, y = test.group.label, color = train.groups)) +
    geom_vline(xintercept = 0, color = "grey50", linewidth = 0.6) +
    geom_pointrange(aes(xmin = delta_lower, xmax = delta_upper),
                    position = position_dodge(width = 0.6),
                    linewidth = 0.6) +
    facet_wrap(~algorithm, nrow = 1) +
    scale_color_manual(values = colors) +
    labs(
      title = tit,
      x = delta_label,
      y = "Test group",
      color = "Train groups"
    ) +
    theme_bw() +
    theme(
      panel.spacing = grid::unit(0.6, "lines"),
      axis.text.y = element_text(size = 10)
    )
  
  out.file <- sub("RData$", "pooling-gain.png", comparison.RData)
  ggsave(out.file, gg, height = 2, width = 4, units = "in", dpi = 300)
  message("Wrote: ", out.file)
  
  invisible(list(summary = summary.dt, delta = delta.dt))
}

## ---- Example runs (one at a time) ----
## DIM +C between sites (test group = site)
load("qsip_pc2_all_new-dim.C.between.sites.RData")
plot_pooling_gain(
  bmr,
  comparison.RData = "qsip_pc2_all_new-dim.C.between.sites.RData",
  tit = "Pooling gain by site (DIM +C)",
  algorithms = c("cv_glmnet", "ranger"),
  train_groups = c("other", "all"),
  use_log10 = FALSE,
  label_col = "site"
)

## DIM +CN between sites (test group = site)
# load("qsip_pc2_all_new-dim.CN.combine.sites.RData")
# plot_pooling_gain(
#   bmr,
#   comparison.RData = "qsip_pc2_all_new-dim.CN.combine.sites.RData",
#   tit = "Pooling gain by site (DIM +CN)",
#   algorithms = c("cv_glmnet", "ranger"),
#   train_groups = c("other", "all"),
#   use_log10 = FALSE,
#   label_col = "site"
# )

## DIM compare sites (all treatments within dim; test group = site)
# load("qsip_pc2_all_new-dim.compare.sites.RData")
# plot_pooling_gain(
#   bmr,
#   comparison.RData = "qsip_pc2_all_new-dim.compare.sites.RData",
#   tit = "Pooling gain by site (DIM all treatments)",
#   algorithms = c("cv_glmnet", "ranger"),
#   train_groups = c("other", "all"),
#   use_log10 = FALSE,
#   label_col = "site"
# )

## DIM compare treatments (test group = treatment)
# load("qsip_pc2_all_new-dim.compare.treatments.RData")
# plot_pooling_gain(
#   bmr,
#   comparison.RData = "qsip_pc2_all_new-dim.compare.treatments.RData",
#   tit = "Pooling gain by treatment (DIM)",
#   algorithms = c("cv_glmnet", "ranger"),
#   train_groups = c("other", "all"),
#   use_log10 = FALSE,
#   label_col = "treatment"
# )
