###############################################################################
message("1 # SETUP")
###############################################################################
rm(list = ls())

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

script_title <- "global"
script_description <- "
Shared setup for the Novel Short Primer 16S workflow.
Loads the local helper file, package libraries, data tables,
shared orders, palettes, and common plotting helpers used by the analysis scripts.
"



set.seed(20240725)
options(scipen = 999)
options(repos = c(CRAN = "https://cloud.r-project.org"))
.libPaths(unique(c("C:/_Repositories/R_libs", .libPaths())))

wd <- get_script_dir()
project_folder <- basename(wd)
setwd(wd)

source(file.path(wd, "helperJ.R"))

required_libraries <- c(
  "tidyverse",
  "vegan",
  "gridExtra",
  "compositions",
  "ggh4x",
  "ggtext",
  "patchwork",
  "readxl",
  "broom",
  "scales",
  "openxlsx",
  "svglite",
  "insect"
)

load_libraries(required_libraries)
qPrint(project_folder)

###############################################################################
message("36 # END SETUP")
###############################################################################

###############################################################################
message("40 # LOAD FILES")
###############################################################################

qLoad("data_tables/meta.csv")
qLoad("data_tables/matrix_names.csv")
qLoad("data_tables/ASV_taxa.csv")

if (exists("meta")) {
  meta_names <- names(meta)
  qPrint(meta_names)
}

###############################################################################
message("54 # END LOAD FILES")
###############################################################################

###############################################################################
message("58 # VARIABLES")
###############################################################################

sample_type_order <- c("Feces", "Skin", "Soil", "Wastewater")
sample_types <- sample_type_order

select_pools <- c("Standard", "TruncatedA")
primer_sets <- c("Ind", "Pool")

Ind <- c(
  "StandardO", "TruncatedFO", "V1.2-12.CC", "V2.2-13.CA", "V3.3-14.TC",
  "V4.4-15.TA", "V5.5-16.xA", "V6.6-15.xC", "V7.7-17.xA", "V8.8-17.xC"
)
Pool <- c("Standard", "TruncatedA")

all_primers <- c(
  "Standard", "StandardO", "TruncatedA", "TruncatedF", "TruncatedFO",
  "TruncatedR", "V1.2-12.CC", "V2.2-13.CA", "V3.3-14.TC", "V4.4-15.TA",
  "V5.5-16.xA", "V6.6-15.xC", "V7.7-17.xA", "V8.8-17.xC"
)

taxa_levels <- c("Phylum", "Family", "Genus", "ASV")
taxa_plural <- c("Phyla", "Families", "Genera", "ASVs")
names(taxa_plural) <- taxa_levels

if (exists("meta")) {
  meta_data <- names(meta)
  meta_primer <- meta %>%
    select(Primer, primer_label) %>%
    distinct()

  primer_order <- meta_primer$Primer
  primer_label_order <- unique(meta_primer$primer_label)

  qPrint(primer_order)
  qPrint(primer_label_order)
}

set_output()

###############################################################################
message("101 # END VARIABLES")
###############################################################################

###############################################################################
message("105 # PALETTES")
###############################################################################

color_palette <- c(
  "515 Forward Primers" = "#339943",
  "806 Reverse Primers" = "#8B1A1A",
  "Tr" = "#E3B22B",
  "St" = "#666666",
  "Truncated" = "#E3B22B",
  "Standard" = "#666666",
  "V1" = "#8B1A1A",
  "V2" = "#339943",
  "V3" = "#FF4593",
  "V4" = "#206699",
  "V5" = "#00FF00",
  "V6" = "#EE7788",
  "V7" = "#7FFFD4",
  "V8" = "#DD88DD",
  "Archaea" = "#FF4040",
  "Bacteria" = "#000080",
  "general_background" = "#FAFCFF",
  "Fecal" = "lavender",
  "Skin" = "#E7FFE1",
  "Soil" = "tan",
  "Water" = "aliceblue",
  "Feces" = "lavender",
  "Wastewater" = "aliceblue",
  "sig" = "red",
  "ns" = "black",
  "ASVs" = "black",
  "Genera" = "#E3ABFF",
  "Phyla" = "skyblue1",
  "ASV" = "black",
  "Genus" = "#E3ABFF",
  "Phylum" = "skyblue1",
  "Mean Tm" = "firebrick2",
  "white" = "white"
)

shape_palette <- c(
  "515 Forward Primers" = 19,
  "806 Reverse Primers" = 17
)

palette_color <- color_palette
palette_shape <- shape_palette
palette_label <- c(
  "Tr" = "Truncated",
  "St" = "Standard"
)

###############################################################################
message("149 # END PALETTES")
###############################################################################

###############################################################################
message("153 # THEMES")
###############################################################################

common_theme <- theme(axis.text.x = element_text(angle = 90)) +
  theme(
    plot.margin = margin(5.5, 5.5, 45, 5.5),
    plot.background = element_rect(fill = color_palette["general_background"]),
    legend.background = element_rect(fill = color_palette["general_background"]),
    legend.position = c(0.5, -0.09),
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    caption.title = element_text(hjust = 0.5)
  )

gPlot <- function(p) {
  p <- p +
    scale_color_manual(values = color_palette) +
    guides(color = "none") +
    scale_x_log10(labels = log_labels, breaks = log_breaks) +
    common_theme +
    labs(
      title = taxa_levs,
      x = "Mean Reads",
      y = "",
      caption = "mean differential abunance"
    )
  print(p)
  p
}

###############################################################################
message("184 # END THEMES")
###############################################################################

###############################################################################
message(paste("188 # FINISHED", script_title))
###############################################################################
