make_sankey_plot <- function(
  # order = "Always Symptomatic", "Symptoms Prevented", "Both Prevented", 
  #         "Always Asymptomatic", "Asymptomatic Prevented", "Never Infected")
  ps_proportions = c(0.22, 0.14, 0.05, 0.14, 0.10, 0.2),
  # flow colors
  flow_colors = c(
	  "Both Prevented" = "#ff684c", 
	  "Asymptomatic Prevented" = "#ffda66", 
	  "Never Infected" = "#309143", 
	  "Symptoms Prevented" = "#e03531", 
	  "Always Asymptomatic" = "#f0bd27", 
	  "Always Symptomatic" = "#b60a1c" 
	),
	node_colors = c(
	  "Symptomatic_Placebo" = "#b60a1c",
	  "Asymptomatic_Placebo" = "#e39802", # (grayed out)
	  "Uninfected_Placebo"   = "#309143", # (grayed out)
	  "Symptomatic_Vaccine"  = "#b60a1c",
	  "Asymptomatic_Vaccine" = "#e39802",
	  "Uninfected_Vaccine"   = "#309143"
	),
	ps_label_x_adjust = c("Both Prevented" = -0.26),
	ps_label_y_adjust = c("Both Prevented" = 0.135),
	label_colors = c(
		"Both Prevented" = "black", 
		"Asymptomatic Prevented" = "black",
		"Never Infected" = "black",
		"Symptoms Prevented" = "black",
		"Always Asymptomatic" = "black",   
		"Always Symptomatic" = "black"
	),
  gray_out_nodes = NULL, # names(node_colors) gives options
  gray_out_flows = NULL, #names(flow_colors) gives options
  gray_color_nodes_and_flows = "#DBE2E9",
  gray_color_flow_labels = "gray70",
  ...
){

  if(!is.null(gray_out_nodes)){
    node_colors[gray_out_nodes] <- gray_color_nodes_and_flows
  }
  if(!is.null(gray_out_flows)){
    flow_colors[gray_out_flows] <- gray_color_nodes_and_flows
    label_colors[gray_out_flows] <- gray_color_flow_labels
  }

	flow_data <- data.frame(
	  Placebo = c("Symptomatic", "Symptomatic", "Symptomatic", "Asymptomatic", "Asymptomatic", "Uninfected"),
	  Vaccine = c("Symptomatic", "Asymptomatic", "Uninfected", "Asymptomatic", "Uninfected", "Uninfected"),
	  value = ps_proportions,
	  label = c("Always Symptomatic", "Symptoms Prevented", "Both Prevented", 
	            "Always Asymptomatic", "Asymptomatic Prevented", "Never Infected"),
	  flow_id = 1:6
	)

	# to control vertical order (bottom to top)
	node_levels <- c("Uninfected", "Asymptomatic", "Symptomatic")
	flow_data$Placebo <- factor(flow_data$Placebo, levels = node_levels)
	flow_data$Vaccine <- factor(flow_data$Vaccine, levels = node_levels)

	# Convert to long format with ggsankey
	sankey_long <- flow_data %>%
	  make_long(Placebo, Vaccine, value = value)

	sankey_long$node <- factor(sankey_long$node, levels = node_levels)
	sankey_long$next_node <- factor(sankey_long$next_node, levels = node_levels)

	sankey_long <- sankey_long %>%
  mutate(node_side = paste(node, x, sep = "_"))

  sankey_plot <- ggplot(
		sankey_long, aes(
			x = x, 
			next_x = next_x, 
			node = node, 
			next_node = next_node, 
			value = value,
			fill = node_side
		)) +
	  geom_sankey(flow.alpha = 0.7, flow.fill = rep(flow_colors, each = 600)) +
	  geom_sankey_text(
	    aes(label = node),
	    angle = 90,         # Turn text sideways
	    fontface = "bold",  # Bold font
	    color = "white",    # Text color
	    size = 3.5           # Text size
	  ) + 
	  scale_fill_manual(values = node_colors) +
	  theme_sankey(base_size = 10) +
	  scale_x_discrete(position = "top") +
	  theme(
	    axis.text.x.top = element_text(
	      color = "black", face = "bold", size = 12, 
	      margin = margin(b = -10)
	    ),
	    axis.ticks.x.top = element_blank(),
	    axis.line.x.top = element_blank(),
	    axis.text.x.bottom = element_blank(),
	    axis.ticks.x.bottom = element_blank(),
	    axis.line.x.bottom = element_blank(),
	    legend.position = "none",
	    plot.margin = margin(t = -10, r = -40, b = -10, l = -55) 
	  ) + 
	  xlab("") +
	  ylab("")

  plot_data <- ggplot_build(sankey_plot)
	flow_polygons <- plot_data$data[[1]]

	flow_labels <- flow_polygons %>%
	  group_by(group) %>%
	  summarize(
	    x = mean(x),
	    y = mean(y),
	    .groups = "drop"
	  ) %>%
	  mutate(
	    label = c(
	      "Both Prevented", "Asymptomatic Prevented", "Never Infected", 
	      "Symptoms Prevented", "Always Asymptomatic", "Always Symptomatic"
	    )
	  )

	ps_x_adjust <- names(ps_label_x_adjust)
	for(ps in ps_x_adjust){
		flow_labels$x[flow_labels$label == ps] <- flow_labels$x[flow_labels$label == ps] + ps_label_x_adjust[ps]
	}
  ps_y_adjust <- names(ps_label_y_adjust)
	for(ps in ps_y_adjust){
		flow_labels$y[flow_labels$label == ps] <- flow_labels$y[flow_labels$label == ps] + ps_label_y_adjust[ps]
	}

	sankey_plot_with_labels <- sankey_plot +
	  geom_text(
	    data = flow_labels,
	    aes(x = x, y = y, label = label),
	    fontface = "bold",
	    size = 2.7,
	    color = label_colors,
	    inherit.aes = FALSE
	  )

	return(sankey_plot_with_labels)
}