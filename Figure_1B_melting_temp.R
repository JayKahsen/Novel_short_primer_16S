###############################################################################
message("1 # SETUP")
###############################################################################

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)))
  }
  getwd()
}

source(file.path(get_script_dir(), "_globalStuff.R"))

script_title <- "figure_1b_melting_temp"
set_output(script_title)

script_description <- "
Creates the melting temperature comparison plot for the novel short primer manuscript.
Uses data_tables/melting_temperatures.csv and saves the figure to output_plot.
"

set_output(script_title)
plot_title <- "melting temperature"

###############################################################################
message("25 # END SETUP")
###############################################################################

###############################################################################
message("29 # PLOT SETTINGS")
###############################################################################

strip_text_size <- 14

common_theme <- theme(
  legend.position = "none",
  plot.margin = margin(20, 10, 10, 10),
  axis.ticks.x = element_blank(),
  axis.text.x = element_blank(),
  legend.title = element_blank(),
  axis.text.y = element_text(size = 14, face = "bold"),
  plot.background = element_rect(fill = color_palette["general_background"]),
  axis.title.y = element_text(size = 14, face = "bold"),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  strip.text = element_text(size = 14),
  plot.title = element_blank()
)

gPlot <- function(p) {
  p <- p +
    scale_color_manual(values = color_palette) +
    guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
    common_theme
  print(p)
  p
}

###############################################################################
message("59 # MAIN")
###############################################################################

df <- read.csv("data_tables/melting_temperatures.csv") %>%
  mutate(Direction = if_else(Direction == "Forward", "515 Forward Primers", Direction)) %>%
  mutate(Direction = if_else(Direction == "Reverse", "806 Reverse Primers", Direction)) %>%
  mutate(plot_group = paste0(Direction, Group)) %>%
  ungroup()

annotate_df <- df %>%
  group_by(Direction, Group) %>%
  summarize(mean_melting_temp = mean(MeltTemp), .groups = "drop") %>%
  mutate(metric = "Mean Tm")

###############################################################################
message(paste("74 #", plot_title))
###############################################################################

unique_x_labels <- unique(df$Group)
outline_color_x <- "black"
backgrounds_x <- lapply(
  unique_x_labels,
  function(label) element_rect(color = outline_color_x, fill = color_palette[[label]])
)

p <- ggplot(df, aes(x = Direction, y = MeltTemp)) +
  geom_boxplot(aes(color = Direction), outlier.shape = NA) +
  geom_point(
    aes(color = Direction, shape = Direction),
    position = position_jitter(width = 0.2, height = 0),
    alpha = 0.5,
    size = 3
  ) +
  facet_wrap2(
    ~Group,
    strip = strip_themed(
      background_x = backgrounds_x,
      text_x = elem_list_text(
        color = c("white", rep("black", length(unique_x_labels) - 1)),
        face = rep("bold", length(unique_x_labels)),
        size = rep(strip_text_size, length(unique_x_labels))
      )
    )
  ) +
  ylab("Theoretical Melting Temperature (Tm) °C") +
  xlab("") +
  geom_text(
    data = annotate_df,
    aes(y = 64, label = paste0(sprintf("%.1f", mean_melting_temp), " °C"), color = metric),
    hjust = 0.5,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette)

gPlot(p)

ggsave(file.path(output_plot, paste0(script_title, ".png")), width = 8, height = 6.5)

###############################################################################
message(paste("114 # FINISHED", script_title))
###############################################################################
