suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mlr3)
  library(mlr3resampling)
  library(scales)
})

rdata.file <- "qsip_pc2_all_new-dim.C_CN_control.downsample.sizes0.compare.treatments.RData"
out.file <- "qsip_pc2_all_new-dim.C_CN_control.downsample.sizes0.compare.treatments.downsample-sizes.png"
out.file.byN <- "qsip_pc2_all_new-dim.C_CN_control.downsample.sizes0.compare.treatments.downsample-sizes-byN.png"
out.file.curve <- "qsip_pc2_all_new-dim.C_CN_control.downsample.sizes0.compare.treatments.downsample-learning-curves.png"
out.file.overview <- "qsip_pc2_all_new-dim.C_CN_control.downsample.sizes0.compare.treatments.downsample-overview.png"

algo.cols <- c(cv_glmnet = "#D55E00", featureless = "#0072B2")
subset.cols <- c(same = "#0072B2", other = "#009E73", all = "#E69F00")

base.theme <- theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "grey45"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.justification = "center",
    panel.spacing.x = grid::unit(0.35, "lines"),
    panel.spacing.y = grid::unit(0.2, "lines"),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.box = "vertical",
    legend.box.spacing = grid::unit(0.1, "lines"),
    legend.margin = margin(0, 0, 0, 0)
  )

save_plot <- function(path, plot_obj, width, height) {
  ggsave(path, plot_obj, width = width, height = height, dpi = 400, bg = "white")
}

if (!file.exists(rdata.file)) {
  stop(sprintf("Missing RData file: %s", rdata.file))
}

load(rdata.file)

if (!exists("bmr")) {
  stop("Expected object `bmr` in RData, but it was not found.")
}

score.obj <- mlr3resampling::score(bmr, mlr3::msr("regr.rmse"))
soak_sizes_score <- as.data.table(score.obj)

required.cols <- c("test.subset", "train.subsets", "n.train.groups", "algorithm", "regr.rmse")
missing.cols <- setdiff(required.cols, names(soak_sizes_score))
if (length(missing.cols) > 0L) {
  stop(sprintf("Missing required columns: %s", paste(missing.cols, collapse = ", ")))
}

soak_sizes_score <- soak_sizes_score[
  !is.na(test.subset) & !is.na(train.subsets) & !is.na(n.train.groups)
]
if (!"test.fold" %in% names(soak_sizes_score) && "iteration" %in% names(soak_sizes_score)) {
  soak_sizes_score[, test.fold := iteration]
}

algo.order <- c("cv_glmnet", "featureless")
soak_sizes_score[, algorithm := factor(as.character(algorithm), levels = unique(c(algo.order, as.character(algorithm))))]
soak_sizes_score[, train.subsets := factor(as.character(train.subsets), levels = c("same", "other", "all"))]

soak_sizes_score[, subset.N := paste(train.subsets, n.train.groups)]
levs <- soak_sizes_score[order(train.subsets, n.train.groups), unique(subset.N)]
soak_sizes_score[, subset.N.fac := factor(subset.N, levels = levs)]
levs.byN <- soak_sizes_score[order(n.train.groups, train.subsets), unique(subset.N)]
soak_sizes_score[, N.subset.fac := factor(subset.N, levels = levs.byN)]

if (all(unique(soak_sizes_score$test.subset) %in% c("control", "CN", "C"))) {
  soak_sizes_score[, test.subset := factor(test.subset, levels = c("control", "CN", "C"))]
}

overview.summary <- soak_sizes_score[, .(
  rmse_mean = mean(regr.rmse),
  rmse_sd = sd(regr.rmse)
), by = .(test.subset, algorithm, train.subsets)]

plot.overview <- ggplot(soak_sizes_score, aes(regr.rmse, train.subsets, color = algorithm)) +
  geom_point(
    alpha = 0.18,
    size = 1.1,
    shape = 16,
    position = position_jitter(height = 0.08, width = 0),
    show.legend = FALSE
  ) +
  geom_segment(
    data = overview.summary,
    aes(
      x = rmse_mean - rmse_sd,
      xend = rmse_mean + rmse_sd,
      y = train.subsets,
      yend = train.subsets
    ),
    linewidth = 0.9,
    position = position_dodge(width = 0.45)
  ) +
  geom_point(
    data = overview.summary,
    aes(x = rmse_mean, y = train.subsets),
    position = position_dodge(width = 0.45),
    shape = 21,
    fill = "white",
    size = 2.6,
    stroke = 0.9
  ) +
  facet_wrap("test.subset", nrow = 1, scales = "free_x", labeller = label_both) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_discrete(expand = expansion(mult = c(0.06, 0.06))) +
  scale_color_manual(values = algo.cols, drop = FALSE) +
  labs(
    title = "DIM C/CN/control downsampling overview",
    subtitle = "Fold-level points with mean +/- SD",
    x = "RMSE",
    y = "Train subset",
    color = "Algorithm"
  ) +
  base.theme

save_plot(out.file.overview, plot.overview, width = 9.8, height = 3.9)

sizes.summary <- soak_sizes_score[, .(
  rmse_mean = mean(regr.rmse),
  rmse_sd = sd(regr.rmse)
), by = .(test.subset, algorithm, subset.N.fac, N.subset.fac)]

plot.downsample <- ggplot(
  sizes.summary,
  aes(rmse_mean, subset.N.fac, color = algorithm, group = algorithm)
) +
  geom_segment(
    aes(
      x = rmse_mean - rmse_sd,
      xend = rmse_mean + rmse_sd,
      yend = subset.N.fac
    ),
    linewidth = 0.9,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    position = position_dodge(width = 0.6),
    shape = 21,
    fill = "white",
    size = 2.5,
    stroke = 0.9
  ) +
  facet_wrap("test.subset", labeller = label_both, scales = "free", nrow = 1) +
  scale_x_continuous(
    breaks = breaks_pretty(n = 3),
    labels = label_number(accuracy = 0.001)
  ) +
  scale_color_manual(values = algo.cols, drop = FALSE) +
  labs(
    title = "Downsampling effect by subset then size",
    subtitle = "Mean +/- SD RMSE over folds",
    x = "RMSE",
    y = "Train subset + n.train.groups",
    color = "Algorithm"
  ) +
  base.theme

save_plot(out.file, plot.downsample, width = 10.4, height = 3.9)

plot.downsample.byN <- ggplot(
  sizes.summary,
  aes(rmse_mean, N.subset.fac, color = algorithm, group = algorithm)
) +
  geom_segment(
    aes(
      x = rmse_mean - rmse_sd,
      xend = rmse_mean + rmse_sd,
      yend = N.subset.fac
    ),
    linewidth = 0.9,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    position = position_dodge(width = 0.6),
    shape = 21,
    fill = "white",
    size = 2.5,
    stroke = 0.9
  ) +
  facet_wrap("test.subset", labeller = label_both, scales = "free", nrow = 1) +
  scale_x_continuous(
    breaks = breaks_pretty(n = 3),
    labels = label_number(accuracy = 0.001)
  ) +
  scale_color_manual(values = algo.cols, drop = FALSE) +
  labs(
    title = "Downsampling effect by size then subset",
    subtitle = "Mean +/- SD RMSE over folds",
    x = "RMSE",
    y = "n.train.groups + train subset",
    color = "Algorithm"
  ) +
  base.theme

save_plot(out.file.byN, plot.downsample.byN, width = 10.4, height = 3.9)

curve.summary <- soak_sizes_score[, .(
  rmse_mean = mean(regr.rmse),
  rmse_sd = sd(regr.rmse)
), by = .(test.subset, algorithm, train.subsets, n.train.groups)]

plot.learning.curves <- ggplot(
  curve.summary,
  aes(n.train.groups, rmse_mean, color = train.subsets)
) +
  geom_errorbar(
    aes(
      ymin = rmse_mean - rmse_sd,
      ymax = rmse_mean + rmse_sd
    ),
    width = 0,
    alpha = 0.45,
    linewidth = 0.55
  ) +
  geom_line(
    aes(
      x = n.train.groups,
      y = rmse_mean,
      color = train.subsets,
      linetype = algorithm,
      group = interaction(train.subsets, algorithm)
    ),
    linewidth = 0.9
  ) +
  geom_point(
    aes(x = n.train.groups, y = rmse_mean),
    size = 2
  ) +
  facet_wrap("test.subset", nrow = 1, scales = "free_y", labeller = label_both) +
  scale_color_manual(values = subset.cols, drop = FALSE) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  labs(
    title = "Downsampling learning curves",
    subtitle = "Mean +/- SD RMSE over folds",
    x = "Number of training rows (n.train.groups)",
    y = "RMSE",
    color = "Train subset",
    linetype = "Algorithm"
  ) +
  base.theme

save_plot(out.file.curve, plot.learning.curves, width = 10.4, height = 3.8)

cat(sprintf("Saved plot: %s\n", out.file.overview))
cat(sprintf("Saved plot: %s\n", out.file))
cat("Ordered subset.N levels:\n")
print(levs)
cat(sprintf("Saved plot: %s\n", out.file.byN))
cat("Ordered subset.N levels by n.train.groups:\n")
print(levs.byN)
cat(sprintf("Saved plot: %s\n", out.file.curve))
