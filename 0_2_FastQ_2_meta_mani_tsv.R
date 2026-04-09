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

script_title <- "fastq_manifest_and_metadata"
script_description <- "
Creates the QIIME manifest.tsv and metadata.tsv files from the repository metadata table.
The script assumes merged files already exist and converts the local metadata
into the tab-delimited QIIME input files.
"

linux_path <- paste0("/home/rushgmcf/Jay/", project_folder, "/demultiplexed_seqs/")

###############################################################################
message("21 # END SETUP")
###############################################################################

###############################################################################
message("25 # BUILD MANIFEST")
###############################################################################

manifest_source <- read.csv("data_tables/meta.csv", check.names = FALSE)

manifest_df <- manifest_source %>%
  transmute(
    `sample-id` = sample_ID,
    `absolute-filepath` = paste0(linux_path, file_name)
) %>%
  arrange(`sample-id`) %>%
  ungroup()

manifest_df %>%
  select(`sample-id`, `absolute-filepath`) %>%
  write.table("manifest.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

###############################################################################
message("38 # END BUILD MANIFEST")
###############################################################################

###############################################################################
message("42 # BUILD METADATA")
###############################################################################

metadata_df <- manifest_source %>%
  select(sample_name, cName, Type, Primer)

metadata_type_row <- c("#q2:types", "categorical", "categorical", "categorical")
metadata_out <- rbind(metadata_type_row, metadata_df)

write.table(metadata_out, "metadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

###############################################################################
message("53 # END BUILD METADATA")
###############################################################################

###############################################################################
message(paste("57 # FINISHED", script_title))
###############################################################################
