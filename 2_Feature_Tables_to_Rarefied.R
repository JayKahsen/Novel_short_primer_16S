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

script_title <- "feature_tables_to_rarefied"
script_description <- "
Builds matrix_names.csv and the raw, rarefied, and filtered matrix layers
from the separate QIIME feature tables. This workflow follows the original
published logic by reading each taxonomic rank table directly rather than
collapsing higher ranks from ASVs.
"

filter_by_counts <- 0
min_relative_abundance <- 1 / 100000
minimum_counts <- c(D4 = 150000, D5 = 100000, D6 = 240000, D7 = 240000, E1 = 50000, E2 = 300000)

raw_matrix_dir <- "data_tables/raw_matrix/"
rarefied_matrix_dir <- "data_tables/rarefied_matrix/"
filtered_matrix_dir <- "data_tables/filtered_matrix/"

dir.create(raw_matrix_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rarefied_matrix_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(filtered_matrix_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
message("32 # END SETUP")
###############################################################################

###############################################################################
message("36 # MATRIX NAMES")
###############################################################################

matrix_names <- expand.grid(
  taxa_levs = taxa_levels,
  primer_sets = primer_sets,
  stringsAsFactors = FALSE
) %>%
  mutate(
    taxa_plural = taxa_plural[taxa_levs],
    raw_path = paste0(raw_matrix_dir, "raw_matrix_", taxa_levs, ".csv"),
    rarefied_path = paste0(rarefied_matrix_dir, "rarefied_matrix_", taxa_levs, "_", primer_sets, ".csv"),
    filtered_path = paste0(filtered_matrix_dir, "filtered_matrix_", taxa_levs, "_", primer_sets, ".csv"),
    file_path = filtered_path
  ) %>%
  arrange(match(taxa_levs, taxa_levels), primer_sets)

write.csv(matrix_names, "data_tables/matrix_names.csv", row.names = FALSE)

###############################################################################
message("52 # END MATRIX NAMES")
###############################################################################

###############################################################################
message("56 # RAW MATRICES")
###############################################################################

sample_reads_df2 <- NULL

for (taxa_levs in taxa_levels) {
  qPrint(taxa_levs)

  if (taxa_levs == "ASV") {
    raw_matrix <- read.delim(
      "data_tables/original/FeatureTable_PostRemoval.tsv",
      check.names = FALSE,
      row.names = 1,
      skip = 1
    ) %>%
      t() %>%
      as.data.frame()
  } else {
    raw_matrix <- read.delim(
      paste0("data_tables/original/Table_", taxa_levs, ".tsv"),
      check.names = FALSE,
      row.names = 1,
      skip = 1
    ) %>%
      t() %>%
      as.data.frame()
  }

  raw_matrix1 <- raw_matrix %>%
    tibble::rownames_to_column(var = "sample_ID") %>%
    left_join(meta, by = "sample_ID") %>%
    select(any_of(meta_names), everything()) %>%
    filter(!is.na(sample_name)) %>%
    imPlode_sample_name()

  sample_reads <- rowSums(raw_matrix1) %>%
    as.data.frame() %>%
    setNames("Reads") %>%
    mutate(taxa_level = paste0(taxa_levs, "_Reads")) %>%
    xPlode_sample_name()

  sample_reads_df2 <- bind_rows(sample_reads_df2, sample_reads)

  raw_matrix1 %>%
    t() %>%
    as.data.frame() %>%
    write.csv(paste0(raw_matrix_dir, "raw_matrix_", taxa_levs, ".csv"))
}

sample_reads_df <- sample_reads_df2 %>%
  pivot_wider(names_from = taxa_level, values_from = Reads) %>%
  imPlode_sample_name()

write.csv(sample_reads_df, "data_tables/sample_reads.csv")

###############################################################################
message("103 # END RAW MATRICES")
###############################################################################

###############################################################################
message("107 # RAREFIED MATRICES")
###############################################################################

for (r in seq_len(nrow(matrix_names))) {
  taxa_levs <- matrix_names[r, "taxa_levs"]
  primer_set <- matrix_names[r, "primer_sets"]
  p_title <- paste0(script_title, "_", taxa_levs, "_", primer_set)

  if (primer_set == "Pool") {
    selected <- Pool
    sample_types <- c("Skin", "Soil")
  }

  if (primer_set == "Ind") {
    selected <- Ind
    sample_types <- sample_type_order
  }

  qPrint(primer_set)
  qPrint(selected)

  matrix_df <- read.csv(matrix_names[r, "raw_path"], check.names = FALSE, row.names = 1) %>%
    t() %>%
    as.data.frame()

  matrix_df1 <- matrix_df
  matrix_df1$sample_count <- rowSums(matrix_df)

  primer_set_df <- matrix_df1 %>%
    xPlode_sample_name() %>%
    filter(Type %in% sample_types) %>%
    filter(Primer %in% selected) %>%
    ungroup()

  rarefied_df <- NULL

  for (t in sample_types) {
    qPrint(t)
    qPrint(p_title)

    type_df1 <- primer_set_df %>%
      filter(Type == t)

    min_sample_count <- minimum_counts[unique(type_df1$Exp)]
    qPrint(min_sample_count)

    type_df <- type_df1 %>%
      filter(sample_count > min_sample_count) %>%
      select(-sample_count) %>%
      imPlode_sample_name() %>%
      as.matrix()

    set.seed(20240725)

    rarefied_type_df <- rrarefy(type_df, sample = min_sample_count) %>%
      as.data.frame()

    rarefied_df <- bind_rows(rarefied_df, rarefied_type_df)
  }

  rarefied_df %>%
    t() %>%
    as.data.frame() %>%
    write.csv(paste0(rarefied_matrix_dir, "rarefied_matrix_", taxa_levs, "_", primer_set, ".csv"))
}

###############################################################################
message("173 # END RAREFIED MATRICES")
###############################################################################

###############################################################################
message("177 # FILTERED MATRICES")
###############################################################################

for (r in seq_len(nrow(matrix_names))) {
  taxa_levs <- matrix_names[r, "taxa_levs"]
  primer_set <- matrix_names[r, "primer_sets"]
  p_title <- paste0(script_title, "_", taxa_levs, "_", primer_set)

  qPrint(p_title)

  matrix_df <- read.csv(matrix_names[r, "rarefied_path"], check.names = FALSE, row.names = 1) %>%
    t() %>%
    as.data.frame()

  feature_names <- names(matrix_df)

  df <- matrix_df %>%
    xPlode_sample_name() %>%
    pivot_longer(cols = all_of(feature_names), names_to = "feature", values_to = "counts") %>%
    group_by(Type) %>%
    mutate(total_counts = sum(counts)) %>%
    group_by(feature, Type) %>%
    mutate(feature_counts = sum(counts)) %>%
    ungroup() %>%
    mutate(rel_abund = feature_counts / total_counts) %>%
    filter(rel_abund > min_relative_abundance) %>%
    select(-c(total_counts, feature_counts, rel_abund)) %>%
    ungroup() %>%
    pivot_wider(names_from = feature, values_from = counts, values_fill = list(counts = 0)) %>%
    imPlode_sample_name() %>%
    t() %>%
    as.data.frame()

  df <- df[, colSums(df, na.rm = TRUE) > filter_by_counts, drop = FALSE]

  write.csv(df, paste0(filtered_matrix_dir, "filtered_matrix_", taxa_levs, "_", primer_set, ".csv"))
}

###############################################################################
message("222 # END FILTERED MATRICES")
###############################################################################

###############################################################################
message(paste("226 # FINISHED", script_title))
###############################################################################
