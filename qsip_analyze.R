library(data.table)
library(mlr3)
library(mlr3resampling)
library(ggplot2)
library(grid)

plot_soak_mse_figs <- function(bmr, comparison.RData, tit,
                               p.color = "#E69F00",
                               suffixes = c("all"),
                               label_col = "site") {
  
  # ---- 1) Fold-level score table ----
  score.dt <- as.data.table(mlr3resampling::score(bmr, mlr3::msr("regr.mse")))
  
  # Harmonize names
  if ("test.subset" %in% names(score.dt))   setnames(score.dt, "test.subset",  "test.group")
  if ("train.subsets" %in% names(score.dt)) setnames(score.dt, "train.subsets","train.groups")
  if (!"test.fold" %in% names(score.dt) && "iteration" %in% names(score.dt)) {
    setnames(score.dt, "iteration", "test.fold")
  }
  
  required_cols <- c("regr.mse","algorithm","test.group","train.groups","test.fold")
  miss <- setdiff(required_cols, names(score.dt))
  if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse=", "))
  
  # ---- 2) Rows per test.group ----
  if ("test" %in% names(score.dt)) {
    n.dt <- score.dt[, .(n.data = as.integer(median(vapply(test, length, integer(1))))),
                     by = test.group]
  } else {
    n.dt <- unique(score.dt[, .(test.group)])
    n.dt[, n.data := NA_integer_]
  }
  
  # ---- 3) Build test.group -> label mapping from score.dt ----
  if (!label_col %in% names(score.dt)) {
    stop("label_col='", label_col, "' not found. Available columns:\n",
         paste(names(score.dt), collapse=", "))
  }
  
  # One label per test.group
  label.map <- score.dt[, .(label = get(label_col)[1L]), by = test.group]
  
  # ---- 4) Add facet vars (robust to missing train.groups) ----
  add_facet_vars <- function(DT) {
    if ("algorithm" %in% names(DT)) {
      DT[algorithm == "regr.rpart", algorithm := "rpart"]
      DT[algorithm == "regr.ranger", algorithm := "ranger"]
      DT[algorithm %in% c("regr.kknn", "kknn"), algorithm := "knn"]
      algo_levels <- c("featureless", "knn", "rpart", "cv_glmnet", "ranger")
      algo_levels <- c(algo_levels, setdiff(unique(DT$algorithm), algo_levels))
      DT[, algorithm := factor(algorithm, levels = algo_levels)]
    }
    
    # Join label
    DT[label.map, label := i.label, on = "test.group"]
    DT[is.na(label), label := as.character(test.group)]
    DT[, `Test group` := gsub("_", "\n", paste0("\n", label))]
    
    # Only add Train\ngroups if train.groups exists in this DT
    if ("train.groups" %in% names(DT)) {
      DT[, `Train\ngroups` := paste0("\n", train.groups)]
    }
    
    # Rows annotation
    DT[n.dt, Rows := n.data, on = "test.group"]
    
    invisible(DT)
  }
  
  add_facet_vars(score.dt)
  
  # ---- 5) Summary stats ----
  score.stats <- score.dt[, .(
    regr.mse_mean = mean(regr.mse),
    regr.mse_sd   = sd(regr.mse),
    nfold         = .N
  ), by = .(train.groups, test.group, algorithm)]
  add_facet_vars(score.stats)
  
  # log scale for paired tests
  score.dt[, log10.regr.mse := log10(regr.mse)]
  
  # ---- 6) same reference line per algorithm (only cv_glmnet + ranger) ----
  same.algos <- c("cv_glmnet", "ranger")
  algo_colors <- c(cv_glmnet = p.color, ranger = "#009E73")
  same.dt <- score.stats[train.groups == "same" & algorithm %in% same.algos,
                         .(test.group, algorithm, same_mean = regr.mse_mean)]
  vline.dt <- copy(same.dt)
  add_facet_vars(vline.dt)
  
  # ---- 7) Paired t-tests: per algorithm other/all vs same (cv_glmnet + ranger only) ----
  score.wide.train <- dcast(
    score.dt[algorithm %in% same.algos],
    algorithm + test.fold + test.group ~ train.groups,
    value.var = "log10.regr.mse"
  )
  if (!("same" %in% names(score.wide.train))) stop("Need 'same' results for paired tests.")
  
  score.wide.train.compare <- melt(
    score.wide.train,
    measure.vars = intersect(c("all", "other"), names(score.wide.train)),
    variable.name = "train.groups",
    value.name = "value"
  )
  
  train.p.values <- score.wide.train.compare[
    !is.na(value) & !is.na(same), {
      if (.N < 2) {
        data.table(p.value = NA_real_, estimate = NA_real_)
      } else {
        test.res <- t.test(value, same, alternative = "two.sided", paired = TRUE)
        data.table(p.value = test.res$p.value, estimate = unname(test.res$estimate))
      }
    }, by = .(algorithm, train.groups, test.group)
  ][
    score.stats, on = c("train.groups", "test.group", "algorithm"), nomatch = 0L
  ][
    same.dt, on = c("test.group", "algorithm")
  ]
  train.p.values[, label_x := (regr.mse_mean + regr.mse_sd) * 1.05]
  add_facet_vars(train.p.values)
  
  # ---- 8) Plot config ----
  train.list <- list(all=c("same","other","all"))
  some <- function(DT, suffix) DT[train.groups %in% train.list[[suffix]]]
  n.test <- length(unique(score.dt$test.group))
  
  for (suffix in suffixes) {
    
    gg.stats <- ggplot()+
      ggtitle(tit)+
      theme_bw()+
      theme(
        panel.spacing=grid::unit(0, "lines"),
        axis.text.x=element_text(angle=30, hjust=1)
      )+
      geom_vline(aes(xintercept=same_mean, color=algorithm),
                 data=vline.dt, linewidth=0.6, show.legend=FALSE)+
      geom_text(aes(
        x = label_x, y = algorithm,
        label=ifelse(p.value<0.0001, "p<0.0001", sprintf("p=%.4f", p.value)),
        color=algorithm
      ),
      vjust=-0.8, size=2.0, show.legend=FALSE,
      data=some(train.p.values, suffix))+
      geom_point(aes(regr.mse_mean, algorithm), shape=1, data=some(score.stats, suffix))+
      geom_segment(aes(
        x = regr.mse_mean - regr.mse_sd, xend = regr.mse_mean + regr.mse_sd,
        y = algorithm, yend = algorithm
      ),
      linewidth=1, data=some(score.stats, suffix))+
      geom_segment(aes(
        x = same_mean, xend = regr.mse_mean,
        y = algorithm, yend = algorithm,
        color = algorithm
      ),
      data=some(train.p.values, suffix),
      show.legend=FALSE)+
      geom_blank(aes(regr.mse, algorithm), data=score.dt)+
      facet_grid(`Train\ngroups` ~ Rows + `Test group`,
                 scales="free_x",
                 labeller=label_both)+
      scale_y_discrete(drop = FALSE)+
      scale_color_manual(values = algo_colors, guide = "none")+
      scale_x_log10("Mean squared prediction error on test set\n(mean +/- SD over 10 folds, log scale)")
    
    out.stats <- sub("RData$", paste0(suffix, "-stats.png"), comparison.RData)
    png(out.stats, height=5, width=(n.test+1)*1.5, units="in", res=300)
    print(gg.stats)
    dev.off()
    message("Wrote: ", out.stats)
  }
  
  invisible(list(score=score.dt, stats=score.stats, p_train=train.p.values))
}

## DIM +CN between sites (test group = site)
## DIM compare sites (all treatments within dim; test group = site)
load("qsip_pc2_all_new-dim.compare.treatments.RData")
plot_soak_mse_figs(
   bmr,
   comparison.RData = "qsip_pc2_all_new-dim.compare.treatments.RData",
   tit = "DIM compare treatments (SOAK regression)",
   suffixes = c("all"),
   label_col = "treatment"
 )


if (F){
#load("qsip_pc2_all_new-dim.C.between.sites.RData")
#plot_soak_mse_figs(
#  bmr,
#  comparison.RData="qsip_pc2_all_new-dim.C.between.sites.RData",
#  tit="DIM +C between sites (SOAK regression)",
#  suffixes=c("all"),
#  label_col="site"
# )

## DIM +CN between sites (test group = site)
# load("qsip_pc2_all_new-dim.CN.combine.sites.RData")
# plot_soak_mse_figs(
#   bmr,
#   comparison.RData = "qsip_pc2_all_new-dim.CN.combine.sites.RData",
#   tit = "DIM +CN between sites (SOAK regression)",
#   suffixes = c("all"),
#   label_col = "site"
# )

## DIM compare sites (all treatments within dim; test group = site)
# load("qsip_pc2_all_new-dim.compare.sites.RData")
# plot_soak_mse_figs(
#   bmr,
#   comparison.RData = "qsip_pc2_all_new-dim.compare.sites.RData",
#   tit = "DIM compare sites (SOAK regression)",
#   suffixes = c("all"),
#   label_col = "site"
# )

## DIM compare treatments (test group = treatment)
# load("qsip_pc2_all_new-dim.compare.treatments.RData")
# plot_soak_mse_figs(
#   bmr,
#   comparison.RData = "qsip_pc2_all_new-dim.compare.treatments.RData",
#   tit = "DIM compare treatments (SOAK regression)",
#   suffixes = c("all"),
#   label_col = "treatment"
# )
}
