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

script_title <- "meta_data"
script_description <- "
Builds data_tables/meta.csv from the sample IDs present in the original phylum table.
The merged file name and total QIIME reads are derived directly from the phylum table
sample IDs and counts.
"

meta_columns <- c(
  "sample_name", "sample_ID", "Exp", "Type", "file_name",
  "Rpl", "Primer", "primer_label", "Group", "Forward_Primer",
  "Reverse_Primer", "cName", "QIIME_reads"
)

primer_names <- c(
  "StandardO", "TruncatedFO", "Standard", "TruncatedR", "TruncatedF",
  "TruncatedA", "V1.2-12.CC", "V2.2-13.CA", "V3.3-14.TC", "V4.4-15.TA",
  "V5.5-16.xA", "V6.6-15.xC", "V7.7-17.xA", "V8.8-17.xC"
)
primer_labels <- c("St", "Tr", "St", "Tr", "Tr", "Tr", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8")
primer_groups <- c(rep("Pool", 6), rep("Ind", 8))

recode_labels <- setNames(primer_labels, primer_names)
recode_groups <- setNames(primer_groups, primer_names)
primer_to_file <- c(
  "StandardO" = "Full",
  "TruncatedFO" = "Pool",
  "Standard" = "515F_806R",
  "TruncatedF" = "515F_CasR",
  "TruncatedR" = "CasF_806R",
  "TruncatedA" = "CasF_CasR",
  "V1.2-12.CC" = "CC_2_12",
  "V2.2-13.CA" = "CA_2_13",
  "V3.3-14.TC" = "TC_3_14",
  "V4.4-15.TA" = "TA_4_15",
  "V5.5-16.xA" = "A_5_16",
  "V6.6-15.xC" = "C_6_15",
  "V7.7-17.xA" = "A_7_17",
  "V8.8-17.xC" = "C_8_17"
)

###############################################################################
message("35 # END SETUP")
###############################################################################

###############################################################################
message("39 # LOAD SAMPLE IDS")
###############################################################################

sample_id_df <- read.delim(
  "data_tables/original/Table_Phylum.tsv",
  check.names = FALSE,
  row.names = 1,
  skip = 1
) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample_ID") %>%
  mutate(QIIME_reads = rowSums(across(-sample_ID, as.numeric), na.rm = TRUE)) %>%
  select(sample_ID, QIIME_reads)

###############################################################################
message("54 # END LOAD SAMPLE IDS")
###############################################################################

###############################################################################
message("58 # BUILD META")
###############################################################################

meta_df <- sample_id_df %>%
  separate(sample_ID, into = c("Exp", "Type", "Primer", "Rpl"), sep = "_", remove = FALSE) %>%
  mutate(
    Type = case_when(
      Type == "Fecal" ~ "Feces",
      Type == "Water" ~ "Wastewater",
      TRUE ~ Type
    ),
    primer_label = recode(Primer, !!!recode_labels),
    Group = recode(Primer, !!!recode_groups),
    Forward_Primer = case_when(
      primer_label == "St" ~ "515F",
      Primer %in% c("TruncatedA", "TruncatedF", "TruncatedFO") ~ "tr515F",
      Primer %in% c("TruncatedR") ~ "515F",
      TRUE ~ paste0(primer_label, "-tr515F")
    ),
    Reverse_Primer = case_when(
      Primer %in% c("TruncatedA", "TruncatedR") ~ "tr515R",
      TRUE ~ "806R"
    ),
    sample_name = paste(Exp, Type, Forward_Primer, Reverse_Primer, Rpl, sep = "_"),
    cName = paste(Exp, Type, primer_label, sep = "_"),
    file_name = paste(Exp, Type, recode(Primer, !!!primer_to_file), Rpl, sep = "_"),
    file_name = paste0(file_name, ".fastq")
  ) %>%
  select(any_of(meta_columns)) %>%
  arrange(sample_ID)

write.csv(meta_df, "data_tables/meta.csv", row.names = FALSE)

###############################################################################
message(paste("88 # FINISHED", script_title))
###############################################################################
