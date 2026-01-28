here::i_am("R/make_all_figures.R")

source(here::here("R/make_sankey_plot.R"))
# devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(ggplot2)
library(dplyr)
library(patchwork)

#----------- Overall figure ----------------# 

sankey_plot_overall <- make_sankey_plot()

ggsave(
  filename = here::here("figs/main_sankey.pdf"),
  plot = sankey_plot_overall,
  device = "pdf",
  width = 4.5,
  height = 5,
  units = "in"
)

#----------- VE_Symp ----------------# 

sankey_plot_ve_symp <- make_sankey_plot(
  gray_out_nodes = c(
    "Uninfected_Placebo",
    "Asymptomatic_Placebo"                 
  ),
  gray_out_flows = c(
    "Asymptomatic Prevented", 
    "Never Infected",
    "Always Asymptomatic"
  ),
  ps_label_x_adjust = c(
    "Both Prevented" = -0.16,
    "Symptoms Prevented" = 0.11
  ),
  ps_label_y_adjust = c(
    "Both Prevented" = 0.10,
    "Symptoms Prevented" = -0.01)
)

sankey_plot_ve_symp <-  sankey_plot_ve_symp +
  annotate("text", x = 1.1, y = 0.2, label = "}", size = 20) + 
  annotate("text", x = 1.24, y = 0.18, label = expression(bold(VE[Sym])), size = 3.75)

ggsave(
  filename = here::here("figs/ve_sym_sankey.pdf"),
  plot = sankey_plot_ve_symp,
  device = "pdf",
  width = 4.5,
  height = 5,
  units = "in"
)

#----------- VE_Inf ----------------# 

sankey_plot_ve_inf <- make_sankey_plot(
  gray_out_nodes = c(
    "Uninfected_Placebo"                 
  ),
  gray_out_flows = c(
    "Never Infected"
  ),
  ps_label_x_adjust = c(
		"Both Prevented" = -0.04,
		"Always Asymptomatic" = 0.15
	),
  ps_label_y_adjust = c("Both Prevented" = 0.10)
)

sankey_plot_ve_inf <- sankey_plot_ve_inf +
  annotate("text", x = 1.065, y = 0.12, label = "}", size = 7) + 
  annotate("text", x = 1.08, y = -0.165, label = "}", size = 10.5) +
  geom_curve(
    aes(x = 1.08, y = 0.115, xend = 1.2, yend = 0.02),
    arrow = arrow(length = unit(0.08, "in")),
    curvature = -0.3,   # amount of curve
    linewidth = 0.5,   # thickness
    color = "black"
  ) + 
  geom_curve(
    aes(x = 1.10, y = -0.165, xend = 1.2, yend = -0.06),
    arrow = arrow(length = unit(0.08, "in")),
    curvature = 0.3,   # amount of curve
    linewidth = 0.5,   # thickness
    color = "black"
  ) +
  annotate("text", x = 1.2, y = -0.03, label = expression(bold(VE[Inf])), size = 3.75)

ggsave(
  filename = here::here("figs/ve_inf_sankey.pdf"),
  plot = sankey_plot_ve_inf,
  device = "pdf",
  width = 4.5,
  height = 5,
  units = "in"
)

#----------- VE_AI ----------------# 

sankey_plot_ve_ai <- make_sankey_plot(
  gray_out_nodes = c(
    "Uninfected_Placebo",
    "Symptomatic_Vaccine"              
  ),
  gray_out_flows = c(
    "Never Infected",
    "Always Symptomatic",
    "Both Prevented"
  ),
  ps_label_x_adjust = c(
		"Both Prevented" = -0.16
	),
  ps_label_y_adjust = c("Both Prevented" = 0.10)
)

ggsave(
  filename = here::here("figs/ve_ai_sankey.pdf"),
  plot = sankey_plot_ve_ai,
  device = "pdf",
  width = 4.5,
  height = 5,
  units = "in"
)

#----------- VE_NAI ----------------# 

sankey_plot_ve_nai <- make_sankey_plot(
  gray_out_nodes = c(
    "Uninfected_Placebo",
    "Symptomatic_Vaccine",            
    "Symptomatic_Placebo"              
  ),
  gray_out_flows = c(
    "Never Infected",
    "Always Symptomatic",
    "Both Prevented",
    "Symptoms Prevented"
  ),
  ps_label_x_adjust = c(
		"Both Prevented" = -0.16,
    "Asymptomatic Prevented" = 0.09
	),
  ps_label_y_adjust = c(
    "Both Prevented" = 0.10,
    "Asymptomatic Prevented" = -0.03
  )
)

sankey_plot_ve_nai <- sankey_plot_ve_nai +
  annotate("text", x = 1.08, y = -0.169, label = "}", size = 11.5) +
  annotate("text", x = 1.192, y = -0.179, label = expression(bold(VE["NAI"])), size = 3.75)

ggsave(
  filename = here::here("figs/ve_nai_sankey.pdf"),
  plot = sankey_plot_ve_nai,
  device = "pdf",
  width = 4.5,
  height = 5,
  units = "in"
)

#-------- combined layout -----------#

# left_side <- (sankey_plot_ve_symp + sankey_plot_ve_inf) / (sankey_plot_ve_ai + sankey_plot_ve_nai)



