# ------------------------------------------------------------------
# Figure for sensitivity analysis results
# ------------------------------------------------------------------

here::i_am("PROVIDE/sensitivity_figure.R")

library(ggplot2)
library(dplyr)
library(ggsci)

results <- readRDS(here::here("PROVIDE/results/sensititvity_results.Rds"))

# Combine point estimates with confidence intervals
ve_sens_plot_df <- results$ve_nai_sens %>%
  left_join(results$boot_est$ve_nai_sens_t0_365, by = "delta")

# Plot
sens_plot <- ggplot(ve_sens_plot_df, aes(x = delta, y = ve_nai * 100)) +
  geom_ribbon(aes(ymin = lower * 100, ymax = upper * 100), alpha = 0.2, fill = "#42B540FF") +
  geom_line(size = 1.2, color = "#42B540FF") +
  geom_point(size = 3, color = "#42B540FF") +
  geom_hline(yintercept = results$ve_fit$ve_ai * 100, linetype = "dashed", color = "gray2") +
  annotate(
    "text",
    x = 0.805,
    y = results$ve_fit$ve_ai * 100 + 4,
    hjust = 0,
    label = bquote(VE[AI] == .(round(results$ve_fit$ve_ai * 100, 0)) * "%" ~
                     "(" * .(round(results$boot_est$t0_365$ve_ai$lower * 100, 0)) * "%, " *
                     .(round(results$boot_est$t0_365$ve_ai$upper * 100, 0)) * "%)"),
    size = 5,
    color = "gray2"
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray2") +
  annotate(
    "text",
    x = 1.005,
    y = 0.8 * 100,
    hjust = 0,
    label = bquote(VE[NAI] == .(round(results$ve_nai_sens$ve_nai[results$ve_nai_sens$delta == 1] * 100, 0)) * "%" ~
                     "(" * .(round(results$boot_est$t0_365$ve_nai$lower * 100, 0)) * "%, " *
                     .(round(results$boot_est$t0_365$ve_nai$upper * 100, 0)) * "%)"),
    size = 5,
    color = "gray2"
  ) +
  scale_x_log10(
    breaks = pretty(ve_sens_plot_df$delta),  # or use c(0.1, 0.2, 0.5, 1, 2, 5, 10) for custom ticks
    labels = scales::label_number()          # or scales::label_math() if you want expressions
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    labels = scales::label_percent(scale = 1)
  ) +
  labs(
    #title = expression("Sensitivity Analysis of " * VE[NAI]),
    x = expression(delta),
    y = expression(VE[NAI]~"(95 \u0025 CI)")
  ) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none")

ggsave(here::here("PROVIDE/results/sens_plot_provide.jpg"), sens_plot, width = 12, height = 8, dpi = 300)
