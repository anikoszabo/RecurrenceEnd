library(ggplot2)
library(dplyr)
library(tidyr)

#########################
# Simulation 1 figure
########################
load("Manuscript/Data/Simulation1.RData")

mean_bias_long <- avebias_res1 |>
  pivot_longer(cols = c(`Q25%`, `Q50%`, `Q75%`),
               names_to = "Quantile",
               values_to = "Bias") |>
  mutate(QuantileSF = as.numeric(substr(Quantile, 2, 3)), # for survival function
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

fig1 <- ggplot(mean_bias_long, aes(x = QuantileDFpercent, y = Bias,
                           color = Method,
                           linetype = factor(p_c),
                           shape = factor(p_c),
                           group = interaction(Method, p_c))) +
  geom_point(position = position_dodge(width = 0.3), size = 2.5, alpha = 0.9) +
  geom_line(position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_linetype_manual(values = linetype_vals) +
  scale_shape_manual(values = shape_vals) +
  scale_color_manual(values = method_colors,
                     breaks = ord,
                     labels = c("Full data", "Naive", expression(tau[0.9]), "NPMLE")) +
  facet_grid(
    EN ~ mu_d,
    labeller = label_bquote(
      rows = N(D) == .(EN),
      cols = mu[d] == .(mu_d)
    ),
    scales = "fixed"
  ) +
  theme_bw() +
  labs(x = "Quantile",
       y = "Bias",
       color = "Method",
       linetype = expression(p[c]),
       shape = expression(p[c])) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        strip.text  = element_text(size = 10))

pdf("Manuscript/Figures/Simulation1.pdf")
  fig1
dev.off()

#########################
# Simulation 2 figure
########################
load("Manuscript/Data/Simulation2.RData")

ord <- c("Full data", "Naive", "Threshold",
         "NPMLE unadjusted", "NPMLE adjusted")

bias_long <- bias_res2 |>
  pivot_longer(cols = c(`Q25%`, `Q50%`, `Q75%`),
               names_to = "Quantile",
               values_to = "Bias") |>
  dplyr::mutate(
    Method   = factor(Method, levels = ord),
    QuantileSF = as.numeric(substr(Quantile, 2, 3)), # for survival function
    QuantileDF = 100 - QuantileSF, # for distribution function
    QuantileDFpercent = paste0(QuantileDF, "%") |>
      factor(levels = c("25%","50%","75%"))
  )


method_colors2 <- c(
  "NPMLE unadjusted" = "#E64B35",
  "Full data"        = "#4DBBD5",
  "Naive"            = "#00A087",
  "Threshold"        = "#3C5488",
  "NPMLE adjusted"   = "#F39B7F"
)


fig2 <- ggplot(bias_long, aes(x = QuantileDFpercent, y = Bias, fill = Method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed", linewidth = 1) +
  labs(x = "Quantile", y = "Bias", fill = "Method") +
  theme_minimal() +
  scale_fill_manual(
    values = method_colors2,
    limits = ord, breaks = ord,
    labels = c("Full data", "Naive", expression(tau[0.9]),
               "NPMLE unadjusted", "NPMLE adjusted")
  )

pdf("Manuscript/Figures/Simulation2.pdf")
fig2
dev.off()
