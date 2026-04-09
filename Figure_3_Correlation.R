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

script_title <- "Figure_3_Correlation"
script_description <- "
Creates the Figure 3 correlation plot comparing Standard and Truncated pooled primer
sets across the active taxonomic levels in matrix_names.csv.
"

set_output(script_title)

###############################################################################
message("24 # END SETUP")
###############################################################################

###############################################################################
message("28 # PLOT SETTINGS")
###############################################################################

big_constant_number <- 1000000
default_text_size <- 11
strip_text_size <- 14
axis_title_text_size <- 14
new_text_size <- 4
v_adjust <- 1.4
g_text_color <- "seagreen"
g_x_position <- 0.1

common_theme <- theme(
  plot.title = element_blank(),
  legend.title = element_blank(),
  plot.background = element_rect(fill = color_palette["general_background"]),
  legend.background = element_rect(fill = color_palette["general_background"]),
  legend.text = element_text(size = 14, face = "bold"),
  axis.title = element_text(size = 12, face = "bold"),
  axis.title.x = element_text(size = axis_title_text_size, margin = margin(t = 10)),
  axis.title.y = element_text(size = axis_title_text_size, margin = margin(r = 10)),
  axis.text = element_markdown(face = "bold", size = rel(1.2)),
  plot.caption = element_text(hjust = 0.5)
)

gPlot <- function(p) {
  p <- p +
    common_theme +
    scale_fill_manual(values = color_palette) +
    scale_color_manual(values = color_palette) +
    scale_x_log10(labels = log_labels_bold, breaks = log_breaks) +
    scale_y_log10(labels = log_labels_bold, breaks = log_breaks) +
    labs(
      title = "",
      x = "Reads from Standard Primer Set Amplification",
      y = "Reads from Truncated Primer Set Amplification"
    )
  p
}

###############################################################################
message("69 # HELPER FUNCTIONS")
###############################################################################

build_correlation_rows <- function(target_row) {
  taxa_levs <- matrix_names[target_row, "taxa_levs"]
  primer_set <- matrix_names[target_row, "primer_sets"]

  if (primer_set == "Ind") {
    return(NULL)
  }

  matrix_df <- read.csv(matrix_names[target_row, "file_path"], check.names = FALSE, row.names = 1) %>%
    t() %>%
    as.data.frame() %>%
    xPlode_sample_name() %>%
    imPlode_sample_name()

  feature_names <- names(matrix_df)

  df2 <- matrix_df %>%
    xPlode_sample_name() %>%
    pivot_longer(cols = all_of(feature_names), names_to = "ASV", values_to = "counts") %>%
    group_by(Type, ASV) %>%
    mutate(ASV_counts = sum(counts)) %>%
    filter(ASV_counts > 0) %>%
    group_by(Type, ASV, Primer) %>%
    summarize(Primer_counts = mean(counts), .groups = "drop_last") %>%
    group_by(Type) %>%
    mutate(mean_type_counts = sum(Primer_counts)) %>%
    ungroup() %>%
    pivot_wider(names_from = Primer, values_from = Primer_counts, values_fill = 0) %>%
    group_by(Type) %>%
    mutate(
      unique_standard = sum(TruncatedA == 0),
      unique_Truncated = sum(Standard == 0),
      unique_feature = n_distinct(ASV),
      shared_feature = sum(Standard != 0 & TruncatedA != 0)
    ) %>%
    ungroup() %>%
    mutate(
      taxa = ASV,
      taxa_level = matrix_names$taxa_plural[target_row]
    )

  if (matrix_names[target_row, "taxa_levs"] == "ASV") {
    df2 <- df2 %>%
      select(-taxa) %>%
      left_join(ASV_taxa %>% select(ASV, taxa), by = "ASV")
  }

  df2 %>%
    mutate(
      Domain = sapply(taxa, extract_domain),
      more = case_when(
        TruncatedA > Standard ~ "Truncated",
        TruncatedA < Standard ~ "Standard",
        TRUE ~ "equal"
      )
    )
}

###############################################################################
message("125 # BUILD DATA")
###############################################################################

active_rows <- which(matrix_names$taxa_levs != "Family")

correlation_input <- bind_rows(lapply(active_rows, build_correlation_rows))
write.csv(correlation_input, file.path(output_data, "correlation.csv"), row.names = FALSE)

correlation_df <- correlation_input %>%
  group_by(taxa_level, Type) %>%
  summarise(
    correlation = cor(Standard, TruncatedA),
    p_value = cor.test(Standard, TruncatedA)$p.value,
    unique_standard_counts = sum(if_else(TruncatedA == 0, Standard, 0)),
    unique_truncated_counts = sum(if_else(Standard == 0, TruncatedA, 0)),
    total_standard_counts = sum(Standard),
    total_truncated_counts = sum(TruncatedA),
    unique_Truncated = sum(Standard == 0),
    unique_Standard = sum(TruncatedA == 0),
    unique_feature = n_distinct(ASV),
    shared_feature = sum(Standard != 0 & TruncatedA != 0),
    .groups = "drop"
  ) %>%
  mutate(correlation_Label = paste0("R^2: ", round(correlation^2, 3)))

###############################################################################
message("150 # STRIP STYLING")
###############################################################################

unique_x_labels <- rev(unique(correlation_input$taxa_level))
outline_color_x <- "black"
s <- paste0(
  'backgrounds_x <- list(element_rect(color=outline_color_x,fill = color_palette["',
  paste0(unique_x_labels, collapse = '"]),element_rect(color=outline_color_x,fill = color_palette["'),
  '"]))'
)
eval(parse(text = s))

unique_y_labels <- unique(correlation_input$Type)
outline_color_y <- "black"
s <- paste0(
  'backgrounds_y <- list(element_rect(color=outline_color_y,fill = color_palette["',
  paste0(unique_y_labels, collapse = '"]),element_rect(color=outline_color_y,fill = color_palette["'),
  '"]))'
)
eval(parse(text = s))

###############################################################################
message("168 # correlation")
###############################################################################

p <- ggplot(correlation_input, aes(x = Standard, y = TruncatedA)) +
  geom_point(aes(color = more)) +
  geom_point(data = correlation_input[correlation_input$Domain == "Archaea", ], aes(color = Domain), shape = 1, stroke = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
  geom_richtext(
    data = correlation_df,
    aes(x = g_x_position, y = Inf, label = paste("Standard Only Features:", unique_Standard)),
    hjust = 0,
    vjust = v_adjust * 1,
    color = color_palette["St"],
    fontface = "bold",
    size = new_text_size
  ) +
  geom_richtext(
    data = correlation_df,
    aes(x = g_x_position, y = Inf, label = paste("Truncated Only Features:", unique_Truncated)),
    hjust = 0,
    vjust = v_adjust * 2,
    color = color_palette["Tr"],
    fontface = "bold",
    size = new_text_size
  ) +
  geom_richtext(
    data = correlation_df,
    aes(x = g_x_position, y = Inf, label = paste("Shared Features:", shared_feature)),
    hjust = 0,
    vjust = v_adjust * 3,
    color = g_text_color,
    fontface = "bold",
    size = new_text_size
  ) +
  geom_richtext(
    data = correlation_df,
    aes(x = g_x_position, y = Inf, label = paste("Total Features:", unique_feature)),
    hjust = 0,
    vjust = v_adjust * 4,
    color = g_text_color,
    fontface = "bold",
    size = new_text_size
  ) +
  geom_richtext(
    data = correlation_df,
    aes(x = Inf, y = 0, label = sprintf("R<sup>2</sup>: %.3f", correlation^2)),
    hjust = 1.1,
    vjust = -1,
    color = "orangered",
    fontface = "bold",
    size = new_text_size
  ) +
  facet_grid2(
    Type ~ taxa_level,
    strip = strip_themed(
      background_x = backgrounds_x,
      background_y = backgrounds_y,
      text_x = elem_list_text(
        color = c("white", rep("black", length(unique_x_labels) - 1)),
        face = rep("bold", length(unique_x_labels)),
        size = rep(strip_text_size, length(unique_x_labels))
      ),
      text_y = elem_list_text(
        face = rep("bold", length(unique_y_labels)),
        size = rep(strip_text_size, length(unique_y_labels))
      )
    )
  ) +
  xlab("Number of Reads (Standard Primer Set)") +
  ylab("Number of Reads (Truncated Primer Set)") +
  ggtitle("")

gPlot(p)

ggsave(file.path(output_plot, "correlation.png"), width = 15, height = 10)

###############################################################################
message(paste("249 # FINISHED", script_title))
###############################################################################
