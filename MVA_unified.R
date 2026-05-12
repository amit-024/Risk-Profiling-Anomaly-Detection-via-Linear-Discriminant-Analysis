################################################################################
##  MULTIVARIATE ANALYSIS OF MILK TRANSPORTATION COST DATA
##
##  Variables
##    FuelType  -  1 = Gasoline, 2 = Diesel
##    X1        -  Fuel cost    (cents / mi)
##    X2        -  Repair cost  (cents / mi)
##    X3        -  Capital cost (cents / mi)
##
##  Required packages
##    MASS        -  LDA / QDA           (ships with base R)
##    MVN         -  Multivariate normality tests
##    psych       -  Factor analysis / KMO / Bartlett
##    heplots     -  Box's M test
##    nortest     -  Anderson-Darling, Lilliefors
##    car         -  Power transformations
##    ellipse     -  Confidence ellipses
##    GPArotation -  Oblique FA rotations
################################################################################


## ── 0. PACKAGES ───────────────────────────────────────────────────────────────

install.packages(c("MVN", "psych", "heplots", "nortest", "car",
                   "ellipse", "GPArotation"))

library(MASS)
library(MVN)
library(psych)
library(heplots)
library(nortest)
library(car)
library(ellipse)
library(GPArotation)


## ── 1. DATA INPUT ─────────────────────────────────────────────────────────────

setwd("C:/Users/Sanjhali.DESKTOP-O4PQA2H/Documents/MVA")

milk           <- read.table("milk_transportation.txt", header = FALSE,
                             col.names = c("FuelType", "X1", "X2", "X3"))
milk$FuelType  <- factor(milk$FuelType,
                         levels = c(1, 2),
                         labels = c("Gasoline", "Diesel"))

vars   <- c("X1", "X2", "X3")
groups <- levels(milk$FuelType)

cat("Dimensions:", nrow(milk), "x", ncol(milk), "\n")
print(head(milk))


## ── SHARED PALETTE ────────────────────────────────────────────────────────────
##  Used consistently across all three functions so every plot in the script
##  uses the same colours for the same groups.

.pal <- c(
  grp1    = "#4C9BE8",   # Gasoline / group 1
  grp2    = "#F4A261",   # Diesel   / group 2
  grp3    = "#66BB6A",
  grp4    = "#AB63FA",
  overall = "black",
  neutral = "grey40",
  heat_lo = "#2166AC",
  heat_hi = "#D6604D"
)

## Convenience: return a colour vector aligned to group membership
.grp_colours <- function(group_vec, groups) {
  pal <- unname(.pal[c("grp1", "grp2", "grp3", "grp4")])
  pal[as.integer(factor(group_vec, levels = groups))]
}


################################################################################
##  FUNCTION 1  run_analysis()
##  Sections: EDA · Univariate normality · Q-Q plots · Multivariate normality
##            Outlier detection (Mahalanobis)
################################################################################

run_analysis <- function(data,
                         group_var,
                         vars,
                         var_labels = NULL,
                         alpha      = 0.05) {

  ## ── parameter setup ─────────────────────────────────────────────────────────
  if (is.null(var_labels))
    var_labels <- setNames(vars, vars)

  groups      <- unique(data[[group_var]])
  n_groups    <- length(groups)
  p           <- length(vars)
  grp_pal     <- unname(.pal[paste0("grp", seq_len(n_groups))])
  col_heat    <- colorRampPalette(c(.pal["heat_lo"], "white", .pal["heat_hi"]))(100)

  ## helper: per-observation colour
  obs_col <- function() .grp_colours(data[[group_var]], groups)


  ## ── 2. EXPLORATORY DATA ANALYSIS ────────────────────────────────────────────

  ### 2.1  Overall summary
  cat("\n---- Overall Summary ----\n")
  print(summary(data))

  ### 2.2  Group descriptive statistics
  cat("\n---- Group Descriptive Statistics ----\n")
  for (g in groups) {
    cat("\nGroup:", g, "\n")
    sub   <- data[data[[group_var]] == g, vars]
    stats <- rbind(
      N    = sapply(sub, length),
      Mean = round(sapply(sub, mean), 3),
      SD   = round(sapply(sub, sd),   3),
      Min  = sapply(sub, min),
      Max  = sapply(sub, max)
    )
    print(stats)
  }

  ### 2.3  Boxplots (overall + per group, one panel per variable)
  par(mfrow = c(1, p), mar = c(5, 4, 4, 2))
  for (v in vars) {
    grp_fac <- factor(
      c(rep("Overall", nrow(data)), as.character(data[[group_var]])),
      levels = c("Overall", as.character(groups))
    )
    vals <- c(data[[v]], data[[v]])
    boxplot(vals ~ grp_fac,
            col    = c(.pal["neutral"], grp_pal)[seq_len(n_groups + 1)],
            main   = var_labels[v],
            xlab   = "",
            ylab   = "Cost (cents / mi)",
            border = "black")
  }
  mtext("Combined Boxplots of Transportation Costs",
        side = 3, line = -2, outer = TRUE, cex = 1.1, font = 2)
  par(mfrow = c(1, 1))

  ### 2.4  Density plots (overall + per group)
  par(mfrow = c(1, p), mar = c(5, 4, 4, 2))
  for (v in vars) {
    d_all    <- density(data[[v]])
    d_grps   <- lapply(groups, function(g) density(data[data[[group_var]] == g, v]))
    ylim_max <- max(d_all$y, sapply(d_grps, function(d) max(d$y)))

    plot(d_all,
         col  = .pal["overall"],
         lwd  = 2,
         main = var_labels[v],
         xlab = "Cost (cents / mi)",
         ylim = c(0, ylim_max))
    polygon(d_all, col = adjustcolor(.pal["overall"], 0.15), border = NA)

    for (i in seq_along(groups)) {
      lines(d_grps[[i]], col = grp_pal[i], lwd = 2)
      polygon(d_grps[[i]], col = adjustcolor(grp_pal[i], 0.30), border = NA)
    }

    legend("topright",
           legend = c("Overall", as.character(groups)),
           col    = c(.pal["overall"], grp_pal[seq_len(n_groups)]),
           lwd    = 2,
           bty    = "n",
           cex    = 0.8)
  }
  mtext("Combined Density Plots",
        side = 3, line = -2, outer = TRUE, cex = 1.1, font = 2)
  par(mfrow = c(1, 1))

  ### 2.5  Scatterplot matrix
  pairs(data[, vars],
        col    = obs_col(),
        pch    = 16,
        cex    = 1.2,
        labels = var_labels[vars],
        main   = paste0(
          "Scatterplot Matrix  (",
          paste(paste0(groups, " = ", grp_pal[seq_len(n_groups)]),
                collapse = ", "),
          ")"))

  ### 2.6  Pearson correlation matrix + heatmaps
  cor_mat    <- cor(data[, vars])
  cat("\n---- Pearson Correlation Matrix ----\n")
  print(round(cor_mat, 4))

  cor_grps <- list()
  for (g in groups) {
    cor_grps[[g]] <- cor(data[data[[group_var]] == g, vars],
                         use = "complete.obs")
    cat(paste0("\n---- Correlation Matrix: ", g, " ----\n"))
    print(round(cor_grps[[g]], 4))
  }

  par(mfrow = c(1, n_groups + 1))

  .draw_heatmap <- function(mat, title) {
    image(1:p, 1:p, mat[p:1, ],
          col  = col_heat,
          xaxt = "n", yaxt = "n",
          zlim = c(-1, 1),
          main = title)
    axis(1, at = 1:p, labels = vars)
    axis(2, at = 1:p, labels = rev(vars))
    for (i in 1:p) for (j in 1:p)
      text(i, p + 1 - j, round(mat[i, j], 2), cex = 1.2)
  }

  .draw_heatmap(cor_mat, "Correlation Heatmap")
  for (g in groups)
    .draw_heatmap(cor_grps[[g]], paste("Correlation Heatmap:", g))

  par(mfrow = c(1, 1))


  ## ── 3. UNIVARIATE NORMALITY ──────────────────────────────────────────────────

  .norm_tests <- function(x, indent = "  ") {
    sw <- shapiro.test(x)
    ks <- ks.test(x, "pnorm", mean = mean(x), sd = sd(x))
    ad <- ad.test(x)
    lf <- lillie.test(x)
    flag <- function(p) ifelse(p < alpha, "<-- NON-NORMAL", "")
    cat(sprintf("%sShapiro-Wilk        : W = %.4f, p = %.4f  %s\n",
                indent, sw$statistic, sw$p.value, flag(sw$p.value)))
    cat(sprintf("%sKolmogorov-Smirnov  : D = %.4f, p = %.4f  %s\n",
                indent, ks$statistic, ks$p.value, flag(ks$p.value)))
    cat(sprintf("%sAnderson-Darling    : A = %.4f, p = %.4f  %s\n",
                indent, ad$statistic, ad$p.value, flag(ad$p.value)))
    cat(sprintf("%sLilliefors          : D = %.4f, p = %.4f  %s\n",
                indent, lf$statistic, lf$p.value, flag(lf$p.value)))
  }

  cat("\n---- Univariate Normality Tests ----\n")
  cat("\n-- Overall --\n")
  for (v in vars) {
    cat(paste0("\nVariable: ", v, "\n"))
    .norm_tests(data[[v]], indent = "  ")
  }

  for (g in groups) {
    cat("\n-- Group:", g, "--\n")
    sub <- data[data[[group_var]] == g, vars]
    for (v in vars) {
      cat(paste0("\n  Variable: ", v, "\n"))
      .norm_tests(sub[[v]], indent = "    ")
    }
  }

  ### 3.1  Q-Q plots
  par(mfrow = c(n_groups + 1, p), mar = c(4, 4, 3, 1))

  for (v in vars) {
    qqnorm(data[[v]],
           main = paste("Q-Q:", v, "| Overall"),
           col  = .pal["neutral"],
           pch  = 16, cex = 0.9)
    qqline(data[[v]], col = .pal["heat_hi"], lwd = 1.5)
  }
  for (i in seq_along(groups)) {
    sub <- data[data[[group_var]] == groups[i], vars]
    for (v in vars) {
      qqnorm(sub[[v]],
             main = paste("Q-Q:", v, "|", groups[i]),
             col  = grp_pal[i],
             pch  = 16, cex = 0.9)
      qqline(sub[[v]], col = .pal["heat_hi"], lwd = 1.5)
    }
  }
  par(mfrow = c(1, 1))


  ## ── 4. MULTIVARIATE NORMALITY (MVN) ─────────────────────────────────────────

  .mvn_block <- function(df, label) {
    cat("\n--", label, "--\n")
    for (test in c("mardia", "hz", "royston", "doornik_hansen")) {
      cat(sprintf("\n%s\n",
                  switch(test,
                         mardia          = "Mardia Test",
                         hz              = "Henze-Zirkler Test",
                         royston         = "Royston Test",
                         doornik_hansen  = "Doornik-Hansen Test")))
      print(MVN::mvn(data = df, mvn_test = test)$multivariate_normality)
    }
  }

  cat("\n---- Multivariate Normality Tests ----\n")
  .mvn_block(data[, vars], "Overall")
  for (g in groups)
    .mvn_block(data[data[[group_var]] == g, vars], paste("Group:", g))


  ## ── 5. OUTLIER DETECTION (Mahalanobis D²) ───────────────────────────────────

  cat("\n---- Multivariate Outlier Detection (Mahalanobis D2) ----\n")

  X_mat      <- as.matrix(data[, vars])
  D2         <- mahalanobis(X_mat, colMeans(X_mat), cov(X_mat))
  data$D2    <- D2
  chi2_crit  <- qchisq(0.975, df = p)

  cat("Chi-square critical value (df =", p, ", 97.5%):",
      round(chi2_crit, 3), "\n")
  outlier_idx <- which(D2 > chi2_crit)
  cat("Outlier rows:", outlier_idx, "\n")
  if (length(outlier_idx) > 0)
    print(data[outlier_idx, ])

  ## Chi-square Q-Q plot
  qqplot(qchisq(ppoints(nrow(data)), df = p), sort(D2),
         main = "Chi-Square Q-Q Plot (Mahalanobis Distances)",
         xlab = paste0("Chi-Square Quantiles (df = ", p, ")"),
         ylab = "Squared Mahalanobis Distances",
         pch  = 16,
         col  = .pal["grp1"])
  abline(0, 1, col = .pal["heat_hi"], lwd = 2)

  ## Scatter with outliers flagged
  plot(data[[vars[1]]], data[[vars[p]]],
       col  = obs_col(),
       pch  = ifelse(D2 > chi2_crit, 8, 16),
       cex  = ifelse(D2 > chi2_crit, 1.6, 1.0),
       main = paste0("Outlier Detection: ", vars[1], " vs ", vars[p],
                     "  (* = outlier)"),
       xlab = var_labels[vars[1]],
       ylab = var_labels[vars[p]])
  text(data[[vars[1]]][outlier_idx],
       data[[vars[p]]][outlier_idx],
       labels = outlier_idx,
       pos    = 3,
       cex    = 0.8,
       col    = "black")
  legend("topright",
         legend = c(as.character(groups), "Outlier"),
         col    = c(grp_pal[seq_len(n_groups)], "black"),
         pch    = c(rep(16, n_groups), 8),
         bty    = "n")

  invisible(data)
}


################################################################################
##  FUNCTION 2  mv_analysis()
##  Sections: Box's M · Hotelling T² · CIs · MANOVA · LDA · QDA
################################################################################

mv_analysis <- function(data,
                        group_var,
                        response_vars,
                        alpha       = 0.05,
                        make_plots  = TRUE,
                        show_output = TRUE) {

  ## ── parameter setup ─────────────────────────────────────────────────────────
  groups  <- unique(data[[group_var]])
  if (length(groups) != 2)
    stop("group_var must have exactly 2 levels.")

  p       <- length(response_vars)
  g1      <- as.matrix(data[data[[group_var]] == groups[1], response_vars])
  g2      <- as.matrix(data[data[[group_var]] == groups[2], response_vars])
  n1      <- nrow(g1)
  n2      <- nrow(g2)

  grp_pal <- unname(.pal[c("grp1", "grp2")])

  ## helper: per-observation colour
  obs_col <- function() .grp_colours(data[[group_var]], groups)


  ## ── Box's M test ─────────────────────────────────────────────────────────────

  boxM_result <- heplots::boxM(data[, response_vars], data[[group_var]])

  if (show_output) {
    cat("\n---- Box's M Test ----\n")
    for (g in groups) {
      cat("\nGroup:", g, "\n")
      print(cov(data[data[[group_var]] == g, response_vars]))
    }
    print(boxM_result)
    cat(ifelse(boxM_result$p.value < alpha,
               "=> Reject H0: unequal covariance matrices  prefer QDA\n",
               "=> Fail to reject H0  LDA assumption holds\n"))
  }


  ## ── Basic statistics ─────────────────────────────────────────────────────────

  xbar1    <- colMeans(g1)
  xbar2    <- colMeans(g2)
  diff_vec <- xbar1 - xbar2
  S1       <- cov(g1)
  S2       <- cov(g2)
  Sp       <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)


  ## ── Hotelling T² ─────────────────────────────────────────────────────────────

  T2     <- as.numeric(
    (n1 * n2 / (n1 + n2)) * t(diff_vec) %*% solve(Sp) %*% diff_vec)
  F_stat <- T2 * (n1 + n2 - p - 1) / ((n1 + n2 - 2) * p)
  df1    <- p
  df2    <- n1 + n2 - p - 1
  p_val  <- pf(F_stat, df1, df2, lower.tail = FALSE)

  if (show_output) {
    cat("\n---- Hotelling T2 Test ----\n")
    cat(sprintf("T2 = %.4f\nF(%d, %d) = %.4f\np-value = %.6f\n",
                T2, df1, df2, F_stat, p_val))
  }


  ## ── Confidence intervals ─────────────────────────────────────────────────────

  F_crit   <- qf(1 - alpha, p, n1 + n2 - p - 1)
  T2_crit  <- p * (n1 + n2 - 2) / (n1 + n2 - p - 1) * F_crit

  ci_hotelling <- matrix(NA, nrow = p, ncol = 2,
                         dimnames = list(response_vars, c("lower", "upper")))
  ci_bon       <- matrix(NA, nrow = p, ncol = 2,
                         dimnames = list(response_vars, c("lower", "upper")))

  if (show_output) cat("\n---- Hotelling Simultaneous CIs ----\n")

  for (j in seq_along(response_vars)) {
    mg               <- sqrt(T2_crit * Sp[j, j] * (1/n1 + 1/n2))
    dj               <- diff_vec[j]
    ci_hotelling[j,] <- c(dj - mg, dj + mg)
    if (show_output)
      cat(sprintf("  %s : (%.3f , %.3f)\n",
                  response_vars[j], dj - mg, dj + mg))
  }

  if (show_output) cat("\n---- Bonferroni CIs ----\n")

  t_bon <- qt(1 - alpha / (2 * p), df = n1 + n2 - 2)

  for (j in seq_along(response_vars)) {
    se          <- sqrt(Sp[j, j] * (1/n1 + 1/n2))
    dj          <- diff_vec[j]
    ci_bon[j, ] <- c(dj - t_bon * se, dj + t_bon * se)
    if (show_output)
      cat(sprintf("  %s : (%.3f , %.3f)\n",
                  response_vars[j], dj - t_bon * se, dj + t_bon * se))
  }

  if (make_plots) {
    plot(seq_along(response_vars), diff_vec,
         xlim = c(0.5, p + 0.5),
         ylim = range(c(ci_hotelling, ci_bon)),
         xaxt = "n",
         pch  = 19,
         col  = "black",
         xlab = "Variable",
         ylab = paste0("Mean difference  (", groups[1], " - ", groups[2], ")"),
         main = "Group Differences with Confidence Intervals")
    axis(1, at = seq_along(response_vars), labels = response_vars)
    abline(h = 0, lty = 2, col = "grey")
    arrows(seq_along(response_vars),
           ci_hotelling[, 1], seq_along(response_vars), ci_hotelling[, 2],
           angle = 90, code = 3, length = 0.1,
           col = "blue", lwd = 2)
    arrows(seq_along(response_vars) + 0.1,
           ci_bon[, 1], seq_along(response_vars) + 0.1, ci_bon[, 2],
           angle = 90, code = 3, length = 0.1,
           col = "red", lwd = 2)
    legend("topright",
           legend = c("Hotelling T2", "Bonferroni"),
           col    = c("blue", "red"),
           lwd    = 2)
  }


  ## ── Overall mean CI and confidence ellipses ──────────────────────────────────

  cat("\n---- Overall Mean CI (Hotelling T2) ----\n")

  X_all        <- as.matrix(data[, response_vars])
  n            <- nrow(X_all)
  xbar_all     <- colMeans(X_all)
  S_all        <- cov(X_all)
  F_crit_all   <- qf(1 - alpha, p, n - p)
  T2_crit_all  <- p * (n - 1) / (n - p) * F_crit_all

  for (j in seq_along(response_vars)) {
    mg <- sqrt(T2_crit_all * S_all[j, j] / n)
    cat(sprintf("  %s : (%.3f , %.3f)\n",
                response_vars[j],
                xbar_all[j] - mg, xbar_all[j] + mg))
  }

  if (make_plots) {

    ## Single 2D ellipse (first two variables)
    ell_overall <- ellipse::ellipse(S_all[1:2, 1:2],
                                    centre = xbar_all[1:2],
                                    t      = sqrt(T2_crit_all / n))
    plot(X_all[, 1], X_all[, 2],
         pch  = 16,
         col  = "grey",
         xlab = response_vars[1],
         ylab = response_vars[2],
         main = "Overall Hotelling T2 Confidence Ellipse (2D)")
    lines(ell_overall, col = "blue", lwd = 2)
    points(xbar_all[1], xbar_all[2], pch = 19, col = "red")
    legend("topright",
           legend = c("Mean", "95% Ellipse"),
           col    = c("red", "blue"),
           lwd    = c(NA, 2),
           pch    = c(19, NA))

    ## Pairwise ellipses
    cat("\n---- Pairwise Hotelling T2 Confidence Ellipses ----\n")
    pairs_idx <- combn(seq_along(response_vars), 2)
    par(mfrow = c(1, ncol(pairs_idx)))

    for (k in seq_len(ncol(pairs_idx))) {
      i     <- pairs_idx[1, k]
      j     <- pairs_idx[2, k]
      mu_2d <- xbar_all[c(i, j)]
      S_2d  <- S_all[c(i, j), c(i, j)]
      ell   <- ellipse::ellipse(S_2d, centre = mu_2d,
                                t = sqrt(T2_crit_all / n))

      all_pts <- rbind(X_all[, c(i, j)], ell)
      plot(X_all[, i], X_all[, j],
           col  = "grey",
           pch  = 16,
           xlim = range(all_pts[, 1]),
           ylim = range(all_pts[, 2]),
           xlab = response_vars[i],
           ylab = response_vars[j],
           main = paste(response_vars[i], "vs", response_vars[j]))
      points(mu_2d[1], mu_2d[2], col = "red",  pch = 19, cex = 1.5)
      lines(ell,                  col = "blue", lwd = 2)
      legend("topright",
             legend = c("Data", "Mean", "95% CI Ellipse"),
             col    = c("grey", "red", "blue"),
             pch    = c(16, 19, NA),
             lwd    = c(NA, NA, 2))
    }
    par(mfrow = c(1, 1))
  }


  ## ── MANOVA ───────────────────────────────────────────────────────────────────

  form       <- as.formula(
    paste0("cbind(", paste(response_vars, collapse = ", "), ") ~ ", group_var))
  manova_fit <- manova(form, data = data)

  if (show_output) {
    cat("\n---- MANOVA ----\n")
    print(summary(manova_fit, test = "Wilks"))
    print(summary(manova_fit, test = "Pillai"))
  }


  ## ── LDA ──────────────────────────────────────────────────────────────────────

  disc_formula <- as.formula(
    paste(group_var, "~", paste(response_vars, collapse = " + ")))

  lda_fit  <- MASS::lda(disc_formula, data = data)
  lda_pred <- predict(lda_fit)

  conf_lda <- table(Actual    = data[[group_var]],
                    Predicted = lda_pred$class)
  cat("\n---- LDA ----\n")
  cat("Confusion Matrix (resubstitution):\n");  print(conf_lda)
  acc_lda  <- sum(diag(conf_lda)) / sum(conf_lda)
  cat(sprintf("APER = %.4f  |  Accuracy = %.2f%%\n",
              1 - acc_lda, acc_lda * 100))

  lda_cv  <- MASS::lda(disc_formula, data = data, CV = TRUE)
  conf_cv <- table(Actual    = data[[group_var]],
                   Predicted = lda_cv$class)
  cat("LOO-CV Confusion Matrix:\n");  print(conf_cv)
  acc_cv  <- sum(diag(conf_cv)) / sum(conf_cv)
  cat(sprintf("LOO-CV Error = %.4f  |  Accuracy = %.2f%%\n",
              1 - acc_cv, acc_cv * 100))

  if (make_plots) {

    ## Discriminant score density
    ld1_g1  <- lda_pred$x[data[[group_var]] == groups[1], 1]
    ld1_g2  <- lda_pred$x[data[[group_var]] == groups[2], 1]
    d1      <- density(ld1_g1)
    d2      <- density(ld1_g2)
    ylim_ld <- c(0, max(d1$y, d2$y) * 1.15)

    plot(d1,
         col  = grp_pal[1],
         lwd  = 2,
         main = "LDA Discriminant Score Distributions",
         xlab = "LD1",
         ylim = ylim_ld)
    polygon(d1, col = adjustcolor(grp_pal[1], 0.30), border = NA)
    lines(d2,   col = grp_pal[2], lwd = 2)
    polygon(d2, col = adjustcolor(grp_pal[2], 0.30), border = NA)
    rug(ld1_g1, col = grp_pal[1], side = 1)
    rug(ld1_g2, col = grp_pal[2], side = 3)
    legend("topright",
           legend = groups,
           fill   = adjustcolor(grp_pal, 0.40),
           border = grp_pal,
           bty    = "n")

    ## Classification scatter
    plot(data[[response_vars[1]]], data[[response_vars[2]]],
         col  = ifelse(data[[group_var]] == groups[1],
                       grp_pal[1], grp_pal[2]),
         pch  = ifelse(lda_pred$class == data[[group_var]], 16, 4),
         cex  = 1.2,
         main = "LDA Classification  (x = misclassified)",
         xlab = response_vars[1],
         ylab = response_vars[2])
    legend("topright",
           legend = c(groups[1], groups[2], "Misclassified"),
           col    = c(grp_pal, "black"),
           pch    = c(16, 16, 4),
           bty    = "n")
  }


  ## ── QDA ──────────────────────────────────────────────────────────────────────

  qda_fit  <- MASS::qda(disc_formula, data = data)
  qda_pred <- predict(qda_fit)

  conf_qda <- table(Actual    = data[[group_var]],
                    Predicted = qda_pred$class)
  cat("\n---- QDA ----\n")
  cat("Confusion Matrix (resubstitution):\n");  print(conf_qda)
  acc_qda  <- sum(diag(conf_qda)) / sum(conf_qda)
  cat(sprintf("APER = %.4f  |  Accuracy = %.2f%%\n",
              1 - acc_qda, acc_qda * 100))

  qda_cv  <- MASS::qda(disc_formula, data = data, CV = TRUE)
  conf_qcv <- table(Actual    = data[[group_var]],
                    Predicted = qda_cv$class)
  cat("LOO-CV Confusion Matrix:\n");  print(conf_qcv)
  acc_qcv  <- sum(diag(conf_qcv)) / sum(conf_qcv)
  cat(sprintf("LOO-CV Error = %.4f  |  Accuracy = %.2f%%\n",
              1 - acc_qcv, acc_qcv * 100))

  if (make_plots) {

    ## QDA uses posterior P(group1 | x) as discrimination axis
    qd1_g1  <- qda_pred$posterior[data[[group_var]] == groups[1], 1]
    qd1_g2  <- qda_pred$posterior[data[[group_var]] == groups[2], 1]
    d1      <- density(qd1_g1)
    d2      <- density(qd1_g2)
    ylim_qd <- c(0, max(d1$y, d2$y) * 1.15)

    plot(d1,
         col  = grp_pal[1],
         lwd  = 2,
         main = "QDA Classification Score Distributions",
         xlab = paste0("P(", groups[1], " | x)"),
         ylim = ylim_qd)
    polygon(d1, col = adjustcolor(grp_pal[1], 0.30), border = NA)
    lines(d2,   col = grp_pal[2], lwd = 2)
    polygon(d2, col = adjustcolor(grp_pal[2], 0.30), border = NA)
    rug(qd1_g1, col = grp_pal[1], side = 1)
    rug(qd1_g2, col = grp_pal[2], side = 3)
    legend("topright",
           legend = groups,
           fill   = adjustcolor(grp_pal, 0.40),
           border = grp_pal,
           bty    = "n")

    ## Classification scatter
    plot(data[[response_vars[1]]], data[[response_vars[2]]],
         col  = ifelse(data[[group_var]] == groups[1],
                       grp_pal[1], grp_pal[2]),
         pch  = ifelse(qda_pred$class == data[[group_var]], 16, 4),
         cex  = 1.2,
         main = "QDA Classification  (x = misclassified)",
         xlab = response_vars[1],
         ylab = response_vars[2])
    legend("topright",
           legend = c(groups[1], groups[2], "Misclassified"),
           col    = c(grp_pal, "black"),
           pch    = c(16, 16, 4),
           bty    = "n")
  }


  ## ── LDA vs QDA summary ───────────────────────────────────────────────────────

  cat("\n---- LDA vs QDA Comparison ----\n")
  comp <- data.frame(
    Method       = c("LDA Resub", "LDA LOO-CV", "QDA Resub", "QDA LOO-CV"),
    Error_Rate   = round(c(1 - acc_lda, 1 - acc_cv,
                           1 - acc_qda, 1 - acc_qcv), 4),
    Accuracy_Pct = round(c(acc_lda, acc_cv,
                           acc_qda, acc_qcv) * 100, 2)
  )
  print(comp)

  if (show_output) {
    cat("\nRecommendation:\n")
    cat(ifelse(boxM_result$p.value < alpha,
               "Use QDA (covariance matrices are unequal)\n",
               "Use LDA (covariance matrices are equal)\n"))
  }


  ## ── return ───────────────────────────────────────────────────────────────────

  return(invisible(list(
    boxM          = boxM_result,
    T2            = T2,
    p_value       = p_val,
    ci_hotelling  = ci_hotelling,
    ci_bonferroni = ci_bon,
    manova        = manova_fit,
    lda           = lda_fit,
    qda           = qda_fit,
    accuracy      = data.frame(
      Method   = c("LDA", "LDA_CV", "QDA", "QDA_CV"),
      Accuracy = c(acc_lda, acc_cv, acc_qda, acc_qcv)
    )
  )))
}


################################################################################
##  FUNCTION 3  mv_reg()
##  Multivariate multiple regression — used to impute detected outliers
################################################################################

mv_reg <- function(data,
                   group_var     = "FuelType",
                   group_value   = "Gasoline",
                   response_vars = c("X1", "X3"),
                   predictor     = "X2",
                   show_output   = TRUE) {

  sub_data    <- data[data[[group_var]] == group_value, ]
  formula_obj <- as.formula(
    paste0("cbind(", paste(response_vars, collapse = ", "), ") ~ ", predictor))
  mv_lm       <- lm(formula_obj, data = sub_data)

  if (show_output) {
    cat("\n---- Multivariate Multiple Regression ----\n")
    cat("Coefficients:\n")
    print(round(coef(mv_lm), 4))
    for (resp in response_vars) {
      cat(paste0("\nSummary for ", resp, ":\n"))
      print(summary(mv_lm)[[paste0("Response ", resp)]])
    }
    cat("\nMANOVA-style test (Pillai):\n")
    print(summary(manova(formula_obj, data = sub_data), test = "Pillai"))
  }

  return(invisible(mv_lm))
}


################################################################################
##  FUNCTION 4  pca_fa_analysis()
##  Sections: PCA · Parallel analysis · FA with four rotations
################################################################################

pca_fa_analysis <- function(data,
                            response_vars,
                            group_var   = NULL,
                            alpha       = 0.05,
                            n_factors   = NULL,
                            make_plots  = TRUE,
                            show_output = TRUE) {

  if (!requireNamespace("psych",       quietly = TRUE))
    stop("Package 'psych' required.")
  if (!requireNamespace("GPArotation", quietly = TRUE))
    stop("Package 'GPArotation' required.")

  ## ── parameter setup ─────────────────────────────────────────────────────────
  p        <- length(response_vars)
  X        <- as.matrix(data[, response_vars])
  groups   <- if (!is.null(group_var)) unique(data[[group_var]]) else NULL
  grp_pal  <- unname(.pal[paste0("grp", seq_len(max(1, length(groups))))])
  bar_cols <- grp_pal[seq_len(p)]

  ## helper: per-observation colour for score plots
  grp_col <- function() {
    if (is.null(group_var))
      rep(.pal["neutral"], nrow(data))
    else
      .grp_colours(data[[group_var]], groups)
  }


  ## ── PCA ──────────────────────────────────────────────────────────────────────

  cat("\n---- Principal Component Analysis ----\n")

  pca_fit <- prcomp(X, scale. = TRUE, center = TRUE)

  if (show_output) {
    print(summary(pca_fit))
    cat("\nLoadings (rotation matrix):\n")
    print(round(pca_fit$rotation, 4))
  }

  pca_var <- pca_fit$sdev^2
  pct_var <- pca_var / sum(pca_var) * 100
  cum_var <- cumsum(pct_var)
  n_pcs   <- length(pca_var)
  pc_labs <- paste0("PC", seq_len(n_pcs))

  if (make_plots) {

    ## Scree + cumulative variance
    par(mfrow = c(1, 2))

    bp <- barplot(pct_var,
                  names.arg = pc_labs,
                  col       = .pal["grp1"],
                  border    = "white",
                  main      = "Scree Plot",
                  ylab      = "% Variance Explained",
                  ylim      = c(0, max(pct_var) * 1.25))
    text(x      = bp,
         y      = pct_var + max(pct_var) * 0.05,
         labels = paste0(round(pct_var, 1), "%"),
         cex    = 0.9)

    plot(seq_len(n_pcs), cum_var,
         type  = "b",
         pch   = 16,
         col   = .pal["heat_hi"],
         main  = "Cumulative Variance Explained",
         xlab  = "Number of PCs",
         ylab  = "Cumulative %",
         ylim  = c(0, 100),
         xaxt  = "n")
    axis(1, at = seq_len(n_pcs), labels = pc_labs)
    abline(h = 80, lty = 2, col = "grey50")
    text(x = 1, y = 82, labels = "80% threshold",
         adj = 0, cex = 0.8, col = "grey50")

    par(mfrow = c(1, 1))

    ## Biplot
    biplot(pca_fit,
           col  = c(.pal["neutral"], .pal["heat_hi"]),
           main = "PCA Biplot")

    ## Score plot PC1 vs PC2
    scores <- as.data.frame(pca_fit$x)
    plot(scores$PC1, scores$PC2,
         col  = grp_col(),
         pch  = 16,
         cex  = 1.3,
         main = "PCA Score Plot: PC1 vs PC2",
         xlab = paste0("PC1 (", round(pct_var[1], 1), "%)"),
         ylab = paste0("PC2 (", round(pct_var[2], 1), "%)"))
    abline(h = 0, v = 0, lty = 2, col = "grey70")
    if (!is.null(group_var))
      legend("topright",
             legend = groups,
             col    = grp_pal[seq_along(groups)],
             pch    = 16,
             bty    = "n")

    ## Loading bar charts for PC1 and PC2
    par(mfrow = c(1, 2))
    for (k in 1:2) {
      barplot(pca_fit$rotation[, k],
              col       = bar_cols,
              main      = paste0("PC", k, " Loadings"),
              ylab      = "Loading",
              names.arg = response_vars,
              border    = "white",
              ylim      = c(-1, 1))
      abline(h = 0)
    }
    par(mfrow = c(1, 1))
  }


  ## ── Factor analysis ──────────────────────────────────────────────────────────

  cat("\n---- Exploratory Factor Analysis ----\n")

  ## Parallel analysis
  cat("\n-- Parallel Analysis --\n")
  fa_parallel  <- psych::fa.parallel(
    X,
    fm   = "ml",
    fa   = "fa",
    main = "Parallel Analysis Scree Plot",
    plot = make_plots
  )
  suggested_nf <- fa_parallel$nfact
  cat("Suggested number of factors:", suggested_nf, "\n")

  nf <- max(1L, as.integer(if (!is.null(n_factors)) n_factors else suggested_nf))
  cat(sprintf("Fitting models with nfactors = %d\n", nf))

  ## Inner helper: fit and report one FA model
  .fit_fa <- function(rotation_name, nfactors) {
    cat(sprintf("\n---- FA | rotation = %s ----\n", rotation_name))

    fit <- tryCatch(
      psych::fa(X, nfactors = nfactors, rotate = rotation_name, fm = "ml"),
      error = function(e) {
        cat("  Model failed:", conditionMessage(e), "\n")
        return(NULL)
      }
    )
    if (is.null(fit)) return(invisible(NULL))

    if (show_output) {
      cat("\nLoadings (cutoff = 0.3):\n")
      print(fit$loadings, cutoff = 0.3)
      cat("\nCommunalities:\n");  print(round(fit$communality,  4))
      cat("Uniquenesses:\n");    print(round(fit$uniquenesses, 4))

      rmsea_val <- tryCatch(round(as.numeric(fit$RMSEA[1]), 4),
                            error = function(e) NA_real_)
      cat("RMSEA:", ifelse(is.na(rmsea_val), "Not available", rmsea_val), "\n")

      var_tbl <- data.frame(
        Factor      = colnames(fit$loadings),
        SS_Loadings = round(colSums(fit$loadings^2), 4),
        Prop_Var    = round(colSums(fit$loadings^2) / p, 4),
        Cum_Var     = round(cumsum(colSums(fit$loadings^2) / p), 4)
      )
      cat("\nVariance explained:\n")
      print(var_tbl, row.names = FALSE)
    }

    if (make_plots) {
      nf_fit <- ncol(fit$loadings)
      par(mfrow = c(1, nf_fit))
      for (k in seq_len(nf_fit)) {
        barplot(fit$loadings[, k],
                col       = bar_cols,
                main      = sprintf("F%d Loadings  [%s]", k, rotation_name),
                ylab      = "Loading",
                names.arg = response_vars,
                border    = "white",
                ylim      = c(-1, 1))
        abline(h   = c(-0.3, 0, 0.3),
               lty = c(2, 1, 2),
               col = c("grey60", "black", "grey60"))
      }
      par(mfrow = c(1, 1))

      if (nf_fit >= 2)
        psych::fa.diagram(fit,
                          main = sprintf("Factor Path Diagram  [%s]",
                                         rotation_name))

      if (!is.null(group_var) && !is.null(fit$scores)) {
        fa_sc <- fit$scores[, 1]
        boxplot(fa_sc ~ data[[group_var]],
                col  = grp_pal[seq_along(groups)],
                main = sprintf("Factor 1 Scores by %s  [%s]",
                               group_var, rotation_name),
                xlab = group_var,
                ylab = "Factor 1 Score")
        stripchart(fa_sc ~ data[[group_var]],
                   vertical = TRUE,
                   method   = "jitter",
                   pch      = 16,
                   col      = adjustcolor("black", 0.4),
                   add      = TRUE)
      }
    }

    return(invisible(fit))
  }

  ## Run all four rotations
  rotations  <- c("none", "varimax", "quartimax", "promax")
  fa_results <- setNames(
    lapply(rotations, .fit_fa, nfactors = nf),
    rotations
  )

  ## Comparison table (Factor 1 loadings across rotations)
  cat("\n---- Rotation Comparison (Factor 1 loadings) ----\n")
  comp_df <- as.data.frame(
    do.call(cbind,
            lapply(rotations, function(rot) {
              fit <- fa_results[[rot]]
              if (is.null(fit)) rep(NA_real_, p) else round(fit$loadings[, 1], 4)
            }))
  )
  colnames(comp_df) <- rotations
  rownames(comp_df) <- response_vars
  print(comp_df)


  ## ── return ───────────────────────────────────────────────────────────────────

  return(invisible(list(
    pca          = pca_fit,
    pct_var      = pct_var,
    cum_var      = cum_var,
    parallel     = fa_parallel,
    fa_none      = fa_results[["none"]],
    fa_varimax   = fa_results[["varimax"]],
    fa_quartimax = fa_results[["quartimax"]],
    fa_promax    = fa_results[["promax"]]
  )))
}


################################################################################
##  ANALYSIS CALLS
################################################################################

## ── EDA + normality + outlier detection ──────────────────────────────────────

var_labels <- c(X1 = "Fuel Cost", X2 = "Repair Cost", X3 = "Capital Cost")

run_analysis(milk,          "FuelType", vars, var_labels)

milk_nooutlier <- milk[-c(9, 21, 41, 47, 56), ]
run_analysis(milk_nooutlier, "FuelType", vars, var_labels)


## ── Transformations ───────────────────────────────────────────────────────────

## Log
milk_log        <- milk[, c("FuelType", vars)]
milk_log[, vars] <- log(milk[, vars])
run_analysis(milk_log, "FuelType", vars, var_labels)

## Box-Cox
bc_model         <- car::powerTransform(milk[, vars], family = "bcPower")
summary(bc_model)
milk_bc          <- as.data.frame(car::bcPower(milk[, vars], bc_model$lambda))
colnames(milk_bc) <- vars
milk_bc$FuelType  <- milk$FuelType
run_analysis(milk_bc, "FuelType", vars, var_labels)

## Yeo-Johnson
yj_model         <- car::powerTransform(milk[, vars], family = "yjPower")
summary(yj_model)
milk_yj          <- as.data.frame(car::yjPower(milk[, vars], yj_model$lambda))
colnames(milk_yj) <- vars
milk_yj$FuelType  <- milk$FuelType
run_analysis(milk_yj, "FuelType", vars, var_labels)


## ── Post-transformation outlier removal ──────────────────────────────────────

milk_log_nooutlier <- milk_log[-c(4, 9), ]
milk_bc_nooutlier  <- milk_bc[-c(4, 9), ]
milk_yj_nooutlier  <- milk_yj[-c(4, 9), ]

run_analysis(milk_log_nooutlier, "FuelType", vars, var_labels)
run_analysis(milk_bc_nooutlier,  "FuelType", vars, var_labels)
run_analysis(milk_yj_nooutlier,  "FuelType", vars, var_labels)


## ── Outlier imputation via multivariate regression ───────────────────────────

milk_model  <- mv_reg(milk_nooutlier)
milk_pred   <- milk
milk_pred[c(9, 21, 41, 47, 56), c("X1", "X3")] <-
  predict(milk_model, newdata = milk[c(9, 21, 41, 47, 56), ])

log_model   <- mv_reg(milk_log_nooutlier)
log_pred    <- milk_log
log_pred[c(4, 9), c("X1", "X3")] <-
  predict(log_model, newdata = milk_log[c(4, 9), ])

bc_model2   <- mv_reg(milk_bc_nooutlier)
bc_pred     <- milk_bc
bc_pred[c(4, 9), c("X1", "X3")] <-
  predict(bc_model2, newdata = milk_bc[c(4, 9), ])

yj_model2   <- mv_reg(milk_yj_nooutlier)
yj_pred     <- milk_yj
yj_pred[c(4, 9), c("X1", "X3")] <-
  predict(yj_model2, newdata = milk_yj[c(4, 9), ])

run_analysis(milk_pred, "FuelType", vars, var_labels)
run_analysis(log_pred,  "FuelType", vars, var_labels)
run_analysis(bc_pred,   "FuelType", vars, var_labels)
run_analysis(yj_pred,   "FuelType", vars, var_labels)


## ── Hotelling T² · MANOVA · LDA · QDA ───────────────────────────────────────

mv_analysis(milk_pred, "FuelType", vars)
mv_analysis(log_pred,  "FuelType", vars)
mv_analysis(bc_pred,   "FuelType", vars)
mv_analysis(yj_pred,   "FuelType", vars)


## ── KMO and Bartlett pre-checks for PCA / FA ─────────────────────────────────

for (nm in c("milk", "milk_log", "milk_bc", "milk_yj",
             "milk_pred", "log_pred", "bc_pred", "yj_pred")) {
  df <- get(nm)
  cat(sprintf("\n---- KMO: %s ----\n", nm))
  print(psych::KMO(df[, vars]))
  cat(sprintf("\n---- Bartlett: %s ----\n", nm))
  print(psych::cortest.bartlett(cor(df[, vars]), n = nrow(df)))
}


## ── PCA + Factor Analysis ─────────────────────────────────────────────────────

for (nm in c("milk", "milk_log", "milk_bc", "milk_yj",
             "milk_pred", "log_pred", "bc_pred", "yj_pred")) {
  cat(sprintf("\n\n==== pca_fa_analysis: %s ====\n", nm))
  pca_fa_analysis(
    data          = get(nm),
    response_vars = vars,
    group_var     = "FuelType"
  )
}
