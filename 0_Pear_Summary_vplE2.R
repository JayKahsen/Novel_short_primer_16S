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

script_title <- "pear_summary"
set_output(script_title)

script_description <- "
Creates file-name, repooling, and aggregate merge summaries from Pear_Summary.txt.
Requires data_tables/Pear_Summary.txt from the read-merging workflow.
"

###############################################################################
message("24 # END SETUP")
###############################################################################

###############################################################################
message("28 # MAIN")
###############################################################################

if (!file.exists("data_tables/Pear_Summary.txt")) {
  stop("Required input missing: data_tables/Pear_Summary.txt")
}

pear_summary <- read.delim("data_tables/Pear_Summary.txt", sep = "\t", header = TRUE) %>%
  select(-X) %>%
  mutate(file_name = gsub("_.*", "", SampleName)) %>%
  mutate(SampleName = file_name) %>%
  separate(SampleName, into = c("Exp", "Type", "Primer_F", "Primer_R", "Rpl"), sep = "-") %>%
  mutate(
    file_name = gsub("-", "_", file_name),
    file_name = paste0(file_name, ".fastq"),
    Primers = paste0(Primer_F, ".", Primer_R),
    cName = paste0(Exp, "_", Type, "_", Primers),
    sample_name = paste0(Exp, "_", Type, "_", Primers, "_", Rpl)
  ) %>%
  select(file_name, sample_name, cName, everything())

write.csv(pear_summary, "data_tables/file_names.csv", row.names = FALSE)
write.csv(pear_summary, paste0(output_plot,"Pear_Summary.csv"), row.names = FALSE)

###############################################################################
message("47 # REPOOLING")
###############################################################################

fine_adj <- 0.07

repooling <- pear_summary %>%
  mutate(
    count_ratio = mean(Raw_Count) / Raw_Count,
    weight = 1,
    initial_volume = 16,
    additional_water = 20,
    repooling_ratio = count_ratio * weight * initial_volume / (initial_volume + additional_water),
    repooling_initial = repooling_ratio * fine_adj / min(repooling_ratio),
    repooling_volume = round(pmin(pmax(repooling_initial, 2), 15), 2)
  )

write.csv(repooling, paste0(output_plot,"repooling.csv"), row.names = FALSE)

###############################################################################
message("64 # MERGE SUMMARY")
###############################################################################

aggregate_merge_percentage <- pear_summary %>%
  group_by(Exp, Type, Primers) %>%
  summarize(
    aggregate_Raw_Count = sum(Raw_Count),
    aggregate_Merged_Count = sum(Merged_Count),
    Merging_Percentage = 100 * aggregate_Merged_Count / aggregate_Raw_Count,
    .groups = "drop"
  )

write.csv(aggregate_merge_percentage, paste0(output_plot,"aggregated_merge_percentage.csv"), row.names = FALSE)

###############################################################################
message(paste("77 # FINISHED", script_title))
###############################################################################
