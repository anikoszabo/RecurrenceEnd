library(ggplot2)
library(dplyr)
library(tidyr)

#########################
# Simulation 1 figure
########################
load("Manuscript/Data/Simulation1_nsim500_boot100.RData")

mean_bias_long <- avebias_res1 |>
  rename(`Q75%_Coverage` = `Coverage_Q75%`,
         `Q50%_Coverage` = `Coverage_Q50%`,
         `Q25%_Coverage` = `Coverage_Q25%`) |>
  pivot_longer(starts_with("Q"), names_to = c("Quantile", ".value"),
               names_sep = "_") |>
  mutate(SD = sqrt(Variance),
         Coverage = 100 * Coverage,
         QuantileSF = as.numeric(substr(Quantile, 2, 3)), # for survival function
         QuantileDF = 100 - QuantileSF, # for distribution function
         QuantileDFpercent = paste0(QuantileDF, "%"))

# Bias Across Quantiles by Method and λ[r]
linetype_vals <- c("solid", "dashed", "dotdash", "twodash")
shape_vals <- c(16, 15, 17, 18)  # 18 = diamond
method_colors <- c(
  "NPMLE"     = "#E64B35",  # red
  "Full data" = "#4DBBD5",  # light blue
  "Naive"     = "#00A087",  # teal
  "Threshold" = "#3C5488"   # dark blue
)

ord <- c("Full data", "Naive", "Threshold", "NPMLE")

metric_label <- function(metric) {
  if_else(metric == "Coverage",
         "Coverage (%)",
         if_else(metric == "SD", "Empirical SD", metric))
}
add_reference <- function(p, metric) {
  if (metric == "Bias") p <- p + geom_hline(yintercept = 0, linetype = "dotted")
  if (metric == "Coverage") {
    p <- p + geom_hline(yintercept = 95, linetype = "dotted")
  }
  p
}

plot_study1 <- function(metric){
  d <- mean_bias_long
  d$value <- d[[metric]]
  p <- d |>
    ggplot(aes(x = QuantileDFpercent,
               y = value,
               color = Method,
               linetype = factor(p_c),
               shape = factor(p_c),
               group = interaction(Method, p_c)))
  p <- add_reference(p, metric) +
    geom_point(position = position_dodge(width = 0.3), size = 2.5, alpha = 0.9) +
    geom_line(position = position_dodge(width = 0.3)) +
    scale_linetype_manual(values = linetype_vals) +
    scale_shape_manual(values = shape_vals) +
    scale_color_manual(values = method_colors,
                     breaks = ord,
                     labels = c("Full data", "Naive", expression(tau[0.9]), "NPMLE")) +
    scale_y_continuous(limits = ~ range(c(0, .x))) +
    guides(color = guide_legend(override.aes = list(shape=NA, linewidth=1))) +
    facet_grid(
      EN ~ mu_d,
      labeller = label_bquote(
        rows = N(D) == .(EN),
        cols = mu[d] == .(mu_d)
      ),
      scales = "fixed"
    ) +
    theme_bw() +
    labs(x = "Evaluation percentile of D",
         y = metric_label(metric),
         color = "Method",
         linetype = expression(p[c]),
         shape = expression(p[c])) +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          strip.text  = element_text(size = 10),
          strip.background = element_rect(fill="white"))
    p
}

fig1_bias <- plot_study1("Bias")
fig1_sd <- plot_study1("SD")
fig1_rmse <- plot_study1("RMSE")
fig1_cov <- plot_study1("Coverage")

pdf("Manuscript/Figures/Simulation1_bias.pdf", width=7, height=5)
  fig1_bias
dev.off()


pdf("Manuscript/Figures/Simulation1_sd.pdf", width=7, height=5)
fig1_sd
dev.off()


pdf("Manuscript/Figures/Simulation1_rmse.pdf", width=7, height=5)
fig1_rmse
dev.off()


pdf("Manuscript/Figures/Simulation1_cov.pdf", width=7, height=5)
fig1_cov
dev.off()

#########################
# Simulation 2 figure
########################
load("Manuscript/Data/Simulation2_nsim500_boot100.RData")

ord <- c("Full data", "Naive", "Threshold",
         "NPMLE unadjusted", "NPMLE adjusted")

bias_long <- avebias_res2 |>
  rename(`Q75%_Coverage` = `Coverage_Q75%`,
         `Q50%_Coverage` = `Coverage_Q50%`,
         `Q25%_Coverage` = `Coverage_Q25%`) |>
  pivot_longer(starts_with("Q"), names_to = c("Quantile", ".value"),
               names_sep = "_") |>
  mutate(Method   = factor(Method, levels = ord),
         SD = sqrt(Variance),
         Coverage = 100 * Coverage,
         QuantileSF = as.numeric(substr(Quantile, 2, 3)), # for survival function
         QuantileDF = 100 - QuantileSF, # for distribution function
         QuantileDFpercent = paste0(QuantileDF, "%"))

bias_extra_long <- bias_long |>
  pivot_longer(c(Bias, SD, RMSE, Coverage),
               names_to = "Metric", values_to = "value") |>
  mutate(Metric = metric_label(Metric))

ref_data <- data.frame(
  Metric = metric_label(c("Bias", "Coverage")),
  Ref = c(0, 95)
)

method_colors2 <- c(
  "NPMLE unadjusted" = "#E64B35",
  "Full data"        = "#4DBBD5",
  "Naive"            = "#00A087",
  "Threshold"        = "#3C5488",
  "NPMLE adjusted"   = "#F39B7F"
)

fig2 <- bias_extra_long |>
  ggplot(aes(x = QuantileDFpercent, y = value, color = Method)) +
  facet_wrap(~Metric, nrow = 2, scales="free_y") +
  geom_hline(aes(yintercept = Ref), linetype = "dotted", data=ref_data) +
  geom_point(size = 2.5, alpha = 0.9) +
  geom_line(aes(group = Method)) +
  scale_linetype_manual(values = linetype_vals) +
  scale_color_manual(values = method_colors2,
                     limits = ord, breaks = ord,
                     labels = c("Full data", "Naive", expression(tau[0.9]),
                                "NPMLE unadjusted", "NPMLE adjusted")) +
  scale_y_continuous(limits = ~ range(c(0, .x))) +
  guides(color = guide_legend(override.aes = list(shape=NA, linewidth=1))) +
  theme_bw() +
  labs(x = "Evaluation percentile of D",
       y = "",
       color = "Method") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        strip.text  = element_text(size = 10),
        strip.background = element_rect(fill="white"))



pdf("Manuscript/Figures/Simulation2.pdf", width=7, height = 5)
fig2
dev.off()

######################
# AL data results
#####################
load("Manuscript/Data/HFdata.RData")

########### Random sample of 50 observations ###########

set.seed(42)

sample_ids <- HF |>
  distinct(patient.id) |>
  sample_n(50)  |>
  pull(patient.id)

# Step 1: Prepare original events
HF_sample <- HF |>
  filter(patient.id %in% sample_ids) |>
  # find length of followup
  group_by(patient.id) |>
  mutate(time = time / 365.24,
         max_time = max(time)) |>
  ungroup() |>
  # ensure order by max_time
  arrange(max_time) |>
  mutate(y_pos = dense_rank(max_time))

# Step 2: line segments
line_segments <- HF_sample |>
  filter(indicator == 0)

# Step 3: events
events <- HF_sample |>
  filter(indicator == 1)

fig3 <- ggplot(events, aes(x=time, y=y_pos)) +
  geom_segment(aes(yend = y_pos, xend = time), x=0, data=line_segments,
               color = "gray80", linewidth = 0.4) +
  #geom_point(color="forestgreen", alpha=0.8) +
  geom_point(color="black", alpha=0.5) +
  scale_x_continuous(
    breaks = seq(0, max(HF$time)),
    limits = c(-1, NA),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    labels = NULL,
    trans = "reverse") +
  labs(
    x = "Time before AL diagnosis date (years)",
    y = ""
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.title.position = "plot",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none"  # Remove legend
  )

pdf("Manuscript/Figures/HFsample.pdf",width = 7, height = 5)
fig3
dev.off()

########## Estimated time to onset #########

library(grid)
library(gridExtra)
library(gridBase)
library(RecurrenceEnd)

set.seed(42)

npmle <- estimate_end(
  Recur(time = time, id = patient.id, event = indicator) ~ age + race + sex,
  method = "NPMLE", bootCI = TRUE, data = HF)
naive_km <- estimate_end(
  Recur(time = time, id = patient.id, event = indicator) ~ age + race + sex,
  method = "naive", threshold = 0, bootCI = TRUE, data = HF)
data_driven_km <- estimate_end(
  Recur(time = time, id = patient.id, event = indicator) ~ age + race + sex,
  method = "quantile", quantile = 0.9, bootCI = TRUE, data = HF)

method_colors3 <- c(
  "NPMLE"     = "#E64B35",  # red
  "Naive"     = "#00A087",  # teal
  "data_driven" = "#7E6148"  # brown
)

# get medians for this and other methods
threshold_km <- estimate_end(
  Recur(time = time, id = patient.id, event = indicator) ~ age + race + sex,
  method = "threshold", threshold = 365, bootCI = TRUE, data = HF)
threshold_km2 <- estimate_end(
  Recur(time = time, id = patient.id, event = indicator) ~ age + race + sex,
  method = "threshold", threshold = 365*2, bootCI = TRUE, data = HF)
threshold_km3 <- estimate_end(
  Recur(time = time, id = patient.id, event = indicator) ~ age + race + sex,
  method = "threshold", threshold = 365*3, bootCI = TRUE, data = HF)

get_median_CI <- function(sf) {
  res <- median(sf)
  c("Median (years)" = sprintf("%.2f", res$quantile/365),
    "95% CI (years)" =
    sprintf("(%.2f, %.2f)", res$lower/365, res$upper/365))
}


tab <- as.data.frame(rbind(
  "Naive" = get_median_CI(naive_km),
  "Threshold = 1 year" = get_median_CI(threshold_km),
  "Threshold = 2 years" = get_median_CI(threshold_km2),
  "Threshold = 3 years" = get_median_CI(threshold_km3),
  "NPMLE" = get_median_CI(npmle)
))
tab$Method <- rownames(tab)

pdf("Manuscript/Figures/ALresults.pdf", width = 8, height = 6)
plot(npmle, do.points = FALSE, ylim = c(0, 1), xaxt = "n", las=1,
     xlab = "Time before AL diagnosis date (years)",
     ylab = "Probability of ongoing HF",
     main = "", bty="n",
     col = method_colors3["NPMLE"],
     conf.int = TRUE,
     conf.col = method_colors3["NPMLE"],
     conf.lty = 2,
    # conf.lwd = 1.5,
     lty = 1, lwd = 2)

lines(naive_km,
      col = method_colors3["Naive"],
      conf.int = TRUE,
      conf.col = method_colors3["Naive"],
      conf.lty = 2, #conf.lwd = 1.5,
      lty = 1, lwd = 2)

lines(data_driven_km,
      col = method_colors3["data_driven"],
      conf.int = TRUE,
      conf.col = method_colors3["data_driven"],
      conf.lty = 2, #conf.lwd = 1.5,
      lty = 1, lwd = 2)

# Add custom x-axis in years (convert days to years)
axis(1, at = seq(0, 4000, by = 365), labels = seq(0, 4000, by = 365) / 365)

legend("bottomleft",
       legend = c("NPMLE", "Naive", expression(tau[0.9] == 0.5 ~ "years")),
       col    = method_colors3[c("NPMLE", "Naive", "data_driven")],
       lty    = 1, lwd = 2, cex = 0.8, bty = "n")


tt <- ttheme_minimal(
  base_size = 9,
  core    = list(bg_params = list(fill = NA, col = NA),
                 fg_params = list(hjust = 0, x = 0.02)),
  colhead = list(bg_params = list(fill = NA, col = NA),
                 fg_params = list(fontface = 2, hjust = 0, x = 0.02))
)
tg <- tableGrob(tab[c("Method", "Median (years)", "95% CI (years)")],
                rows = NULL, theme = tt)

vps <- baseViewports()
pushViewport(vps$plot)
pushViewport(viewport(x = 0.98, y = 0.98, width = 0.45, height = 0.40,
                      just = c("right","top")))
grid.draw(tg)
dev.off()

########### Subgroup plots ###########
#### By gender ###########
set.seed(42)

method_colors2 <- c(
  "NPMLE.m"     = "#E64B35",  # red (male)
  "NPMLE.f"     = "#00A087"  # teal (male)
)


HF.m <- subset(HF, sex == "M")
npmle.m <- estimate_end(
  Recur(time = time, id = patient.id, event = indicator) ~ age + race,
  method = "NPMLE", bootCI = TRUE, data = HF.m)
HF.f <- subset(HF, sex == "F")
npmle.f <- estimate_end(
  Recur(time = time, id = patient.id, event = indicator) ~ age + race,
  method = "NPMLE", bootCI = TRUE, data = HF.f)

pdf("Manuscript/Figures/ALresultsBySex.pdf", width = 8, height = 6)
plot(npmle.m, do.points = FALSE, ylim = c(0, 1), xaxt = "n", las = 1,
     xlab = "Time before AL diagnosis date (years)",
     ylab = "Probability of ongoing HF",
     main = "", bty="n",
     col = method_colors2["NPMLE.m"],
     conf.int = TRUE,
     conf.col = method_colors2["NPMLE.m"],
     conf.lty = 2,
     #conf.lwd = 1.5,
     lty = 1, lwd = 2)

lines(npmle.f,
      col = method_colors2["NPMLE.f"],
      conf.int = TRUE,
      conf.col = method_colors2["NPMLE.f"],
      conf.lty = 2, #conf.lwd = 1.5,
      lty = 1, lwd = 2)


# Add custom x-axis in years (convert days to years)
axis(1, at = seq(0, 4000, by = 365), labels = seq(0, 4000, by = 365) / 365)

# Add legend
legend("topright",
       legend = c("NPMLE (Male)", "NPMLE (Female)"),
       col = method_colors2[c("NPMLE.m","NPMLE.f")],
       lty = 1, lwd = 2, cex = 0.8, bty = "n")
dev.off()

######### By age group ######
set.seed(42)

method_colors2b <- c(
  "NPMLE.ge65"     = "#E64B35",  # red (>=65)
  "NPMLE.lt65"     = "#00A087"  # teal (>=65)
)

HF.ge65 <- subset(HF, !is.na(age) & age >= 65)
HF.lt65 <- subset(HF, !is.na(age) & age < 65)

npmle.ge65 <- estimate_end(Recur(time = time, id = patient.id, event = indicator) ~ age + race,
                           method = "NPMLE", bootCI = TRUE, data = HF.ge65)
npmle.lt65 <- estimate_end(Recur(time = time, id = patient.id, event = indicator) ~ age + race,
                           method = "NPMLE", bootCI = TRUE, data = HF.lt65)

pdf("Manuscript/Figures/ALresultsByAge.pdf", width = 8, height = 6)
plot(npmle.ge65, do.points = FALSE, ylim = c(0, 1), xaxt = "n", las = 1,
     xlab = "Time before AL diagnosis date (years)",
     ylab = "Probability of ongoing HF",
     main = "", bty="n",
     col = method_colors2b["NPMLE.ge65"],
     conf.int = TRUE,
     conf.col = method_colors2b["NPMLE.ge65"],
     conf.lty = 2,
    # conf.lwd = 1.5,
     lty = 1, lwd = 2)

lines(npmle.lt65,
      col = method_colors2b["NPMLE.lt65"],
      conf.int = TRUE,
      conf.col = method_colors2b["NPMLE.lt65"],
      conf.lty = 2, #conf.lwd = 1.5,
      lty = 1, lwd = 2)

axis(1, at = seq(0, 4000, by = 365), labels = seq(0, 4000, by = 365) / 365)

legend("topright",
       legend = c("NPMLE (Age ≥ 65)", "NPMLE (Age < 65)"),
       col = method_colors2b[c("NPMLE.ge65","NPMLE.lt65")],
       lty = 1, lwd = 2, cex = 0.8, bty = "n")
dev.off()
