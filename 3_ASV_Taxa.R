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

script_title <- "asv_taxa"
script_description <- "
Creates ASV_taxa.csv from QIIME taxonomy output and sequence fasta input.
Also summarizes observed ASV sequence lengths in the Ind and Pool ASV matrices.
"

set_output(script_title)

###############################################################################
message("25 # END SETUP")
###############################################################################

###############################################################################
message("29 # BUILD ASV TAXA")
###############################################################################

sequence_data <- readFASTA("data_tables/sequences.fasta") %>%
  dna2char() %>%
  as.data.frame() %>%
  `names<-`("Sequence") %>%
  rownames_to_column(var = "ASV")

ASV_taxa <- read.delim("data_tables/original/taxa_data.tsv") %>%
  slice(-1) %>%
  setNames(c("ASV", "taxa", "Confidence")) %>%
  select(ASV, taxa) %>%
  mutate(
    Domain = str_extract(taxa, "d__([^;\\s]*)"),
    Phylum = str_extract(taxa, "p__([^;\\s]*)"),
    Class = str_extract(taxa, "c__([^;\\s]*)"),
    Order = str_extract(taxa, "o__([^;\\s]*)"),
    Family = str_extract(taxa, "f__([^;\\s]*)"),
    Genus = str_extract(taxa, "g__([^;\\s]*)"),
    Species = str_extract(taxa, "s__([^;\\s]*)")
  ) %>%
  mutate(
    default = case_when(
      !is.na(Species) ~ Species,
      !is.na(Genus) ~ Genus,
      !is.na(Family) ~ Family,
      !is.na(Order) ~ Order,
      !is.na(Class) ~ Class,
      !is.na(Phylum) ~ Phylum,
      !is.na(Domain) ~ Domain,
      TRUE ~ ASV
    )
  ) %>%
  mutate(
    taxa = case_when(is.na(taxa) ~ ASV, TRUE ~ taxa),
    Domain = case_when(is.na(Domain) ~ taxa, TRUE ~ Domain),
    Phylum = case_when(is.na(Phylum) ~ Domain, TRUE ~ Phylum),
    Class = case_when(is.na(Class) ~ Phylum, TRUE ~ Class),
    Order = case_when(is.na(Order) ~ Class, TRUE ~ Order),
    Family = case_when(is.na(Family) ~ Order, TRUE ~ Family),
    Genus = case_when(is.na(Genus) ~ Family, TRUE ~ Genus),
    Species = case_when(is.na(Species) ~ Genus, TRUE ~ Species)
  ) %>%
  left_join(sequence_data, by = "ASV") %>%
  mutate(sequence_length = nchar(Sequence))

write.csv(ASV_taxa, "data_tables/ASV_taxa.csv", row.names = FALSE)

###############################################################################
message("79 # SEQUENCE LENGTH SUMMARIES")
###############################################################################

sequence_length_plot <- ggplot(ASV_taxa, aes(x = sequence_length)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "black") +
  labs(
    title = "Sequence Length Distribution",
    x = "Sequence Length (bp)",
    y = "count"
  ) +
  theme_minimal()

ggsave(file.path(output_plot, "sequence_length_distribution.png"), sequence_length_plot, width = 8, height = 6)

length_table <- ASV_taxa %>%
  count(sequence_length, name = "count") %>%
  arrange(desc(count))

write.csv(length_table, file.path(output_data, "sequence_length_distribution.csv"), row.names = FALSE)

###############################################################################
message("98 # ASV LENGTH BY PRIMER SET")
###############################################################################

summarize_asv_lengths <- function(target_row) {
  matrix_df <- read.csv(matrix_names[target_row, "file_path"], check.names = FALSE, row.names = 1) %>%
    t() %>%
    as.data.frame()

  feature_names <- names(matrix_df)

  length_df <- matrix_df %>%
    xPlode_sample_name() %>%
    pivot_longer(cols = any_of(feature_names), names_to = "ASV", values_to = "counts") %>%
    filter(counts > 0) %>%
    left_join(ASV_taxa, by = "ASV") %>%
    select(sample_name, Type, primer_label, sequence_length, ASV) %>%
    group_by(primer_label, Type) %>%
    count(sequence_length, name = "counts") %>%
    mutate(expected_length = ifelse(sequence_length %in% c("252", "253"), "252+253", as.character(sequence_length))) %>%
    group_by(Type, primer_label, expected_length) %>%
    summarize(counts = sum(counts), .groups = "drop") %>%
    rename(sequence_length = expected_length) %>%
    group_by(primer_label, Type) %>%
    mutate(
      total_counts = sum(counts),
      percent = round(100 * counts / total_counts, 2)
    ) %>%
    select(-counts, -total_counts) %>%
    arrange(desc(percent))

  length_df
}

pool_row <- which(matrix_names$taxa_levs == "ASV" & matrix_names$primer_sets == "Pool")
ind_row <- which(matrix_names$taxa_levs == "ASV" & matrix_names$primer_sets == "Ind")

pool_lengths <- summarize_asv_lengths(pool_row)
ind_lengths <- summarize_asv_lengths(ind_row)

pool_plot <- ggplot(pool_lengths, aes(x = sequence_length, y = percent)) +
  geom_col(fill = "steelblue", color = "black") +
  facet_grid(Type ~ primer_label) +
  labs(
    title = "Sequence Length Distribution",
    x = "Sequence Length (bp)",
    y = "percent"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

ind_plot <- ggplot(ind_lengths, aes(x = sequence_length, y = percent)) +
  geom_col(fill = "steelblue", color = "black") +
  facet_grid(Type ~ primer_label) +
  labs(
    title = "Sequence Length Distribution",
    x = "Sequence Length (bp)",
    y = "percent"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

ggsave(file.path(output_plot, "distribution_pool.png"), pool_plot, width = 8, height = 6)
ggsave(file.path(output_plot, "distribution_ind.png"), ind_plot, width = 12, height = 9)

write.csv(
  pool_lengths %>% pivot_wider(names_from = primer_label, values_from = percent),
  file.path(output_data, "unique_ASVs_length_percent_table_pool.csv"),
  row.names = FALSE
)

write.csv(
  ind_lengths %>% pivot_wider(names_from = primer_label, values_from = percent),
  file.path(output_data, "unique_ASVs_length_percent_table_ind.csv"),
  row.names = FALSE
)

###############################################################################
message(paste("175 # FINISHED", script_title))
###############################################################################
