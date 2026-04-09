################################################################################
# LOCAL HELPER
################################################################################

load_libraries <- function(required_libraries, quietly = TRUE) {
  stopifnot(is.character(required_libraries))

  for (lib in required_libraries) {
    if (!requireNamespace(lib, quietly = TRUE)) {
      install.packages(lib, dependencies = TRUE)
    }

    suppressPackageStartupMessages(
      library(lib, character.only = TRUE, quietly = quietly, warn.conflicts = FALSE)
    )
  }

  message(paste0("Loaded libraries: ", paste(required_libraries, collapse = ", ")))
  invisible(required_libraries)
}

################################################################################

qPrint <- function(var = NULL) {
  if (is.null(var)) {
    env <- parent.frame()
    objects <- ls(envir = env)
    vars <- objects[sapply(objects, function(x) !is.function(get(x, envir = env)))]

    if (length(vars) > 0) {
      last_var_name <- tail(vars, 1)
      var_name <- last_var_name
      var <- get(var_name, envir = env)
    } else {
      env <- .GlobalEnv
      objects <- ls(envir = env)
      vars <- objects[sapply(objects, function(x) !is.function(get(x, envir = env)))]

      if (length(vars) > 0) {
        last_var_name <- tail(vars, 1)
        var_name <- last_var_name
        var <- get(var_name, envir = env)
      } else {
        cat("No variables found in the global environment.\n")
        return(invisible(NULL))
      }
    }
  } else {
    var_name <- deparse(substitute(var))
  }

  value_str <- if (is.vector(var)) paste(var, collapse = ", ") else toString(var)
  value_str <- gsub("/", "", value_str)
  output_string <- sprintf("%s = %s", var_name, value_str)

  cat(paste0("\n\t", output_string, "\t\n"))
  invisible(output_string)
}

################################################################################

qLoad <- function(file_path) {
  if (!file.exists(file_path)) {
    message("\t", file_path, " does not exist.")
    return(invisible(NULL))
  }

  file_name <- tools::file_path_sans_ext(basename(file_path))

  if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
    data <- read.csv(file_path, check.names = FALSE)
  } else if (grepl("\\.xlsx$", file_path, ignore.case = TRUE)) {
    data <- readxl::read_excel(file_path)
  } else {
    stop("Unsupported file type. Please provide a CSV or XLSX file.")
  }

  assign(file_name, data, envir = .GlobalEnv)
  message("\tFile successfully loaded:\t", basename(file_path))
  invisible(data)
}

################################################################################

set_output <- function(subfolder = NULL) {
  base_data <- "output_data"
  base_plot <- "output_plot"

  output_data <- if (is.null(subfolder)) {
    paste0(base_data, "/")
  } else {
    paste0(file.path(base_data, subfolder), "/")
  }

  output_plot <- if (is.null(subfolder)) {
    paste0(base_plot, "/")
  } else {
    paste0(file.path(base_plot, subfolder), "/")
  }

  dir.create(output_data, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_plot, recursive = TRUE, showWarnings = FALSE)

  assign("output_data", output_data, envir = parent.frame())
  assign("output_plot", output_plot, envir = parent.frame())

  invisible(list(output_data = output_data, output_plot = output_plot))
}

################################################################################

xPlode_sample_name <- function(df) {
  message("xPloding")
  rows_df0 <- nrow(df)
  bad_cols <- intersect(colnames(df), meta_names)
  if (length(bad_cols) > 0) {
    stop(paste("duplicate columns found\n", paste(bad_cols, collapse = ", ")))
  }

  df_out <- df %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample_name") %>%
    left_join(meta, by = "sample_name") %>%
    select(any_of(meta_names), everything()) %>%
    ungroup()

  if (nrow(df_out) != rows_df0) {
    stop("meta has duplicate rows")
  }

  df_out
}

################################################################################

imPlode_sample_name <- function(df) {
  message("imPloding")
  df %>%
    as.data.frame() %>%
    tibble::column_to_rownames(var = "sample_name") %>%
    select(-any_of(meta_names))
}

################################################################################

qShannon <- function(x, s = "off") {
  if (s == "on") {
    print("shannon vector")
    print(x)
  }
  x <- x[x != 0]
  if (length(x) == 0) {
    return(NA_real_)
  }
  total <- sum(x)
  -sum(x / total * log(x / total))
}

qEvenness <- function(x, s = "off") {
  if (s == "on") {
    print("evenness vector")
    print(x)
  }
  x <- x[x != 0]
  if (length(x) == 0) {
    return(NA_real_)
  }
  if (length(x) == 1) {
    return(1)
  }
  total <- sum(x)
  -sum(x / total * log(x / total)) / log(length(x))
}

qRichness <- function(x) {
  sum(x > 0)
}

################################################################################

extract_domain <- function(taxa) {
  if (grepl("archaea", taxa, ignore.case = TRUE)) {
    "Archaea"
  } else if (grepl("bacteria", taxa, ignore.case = TRUE)) {
    "Bacteria"
  } else if (grepl("eukaryote", taxa, ignore.case = TRUE)) {
    "Eukaryote"
  } else {
    "unknown"
  }
}

extract_propi <- function(taxa) {
  if (exists("taxa_levs", inherits = TRUE) && taxa_levs == "Phylum") {
    if (grepl("p__Actinobacteriota", taxa, ignore.case = TRUE)) "Propi" else "no"
  } else {
    if (grepl("g__Cutibacterium", taxa, ignore.case = TRUE)) "Propi" else "no"
  }
}

extract_lachno <- function(taxa) {
  if (grepl("f__Lachnospiraceae", taxa, ignore.case = TRUE)) "Lachno" else "no"
}

extract_muri <- function(taxa) {
  if (grepl("f__Muribaculaceae", taxa, ignore.case = TRUE)) "Muri" else "no"
}

################################################################################

log_breaks <- c(0, 10^seq(0, 12, by = 1))
log10_breaks <- seq(0, 12, by = 1)
log2_breaks <- seq(0, 12, by = 1)
fold_change_breaks <- seq(-12, 12, by = 1)

log_labels <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) {
      return(NA)
    }
    if (val == 1) {
      return(expression(10^0))
    }
    as.expression(bquote(10^.(round(log10(val)))))
  })
}

log_labels_bold <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) {
      return(NA)
    }
    exponent <- round(log10(val))
    paste0("<b>10<sup>", exponent, "</sup></b>")
  })
}

fold_change_labels_bold <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) {
      return(NA)
    }
    paste0("<b>2<sup>", val, "</sup></b>")
  })
}

################################################################################

feature_labels <- function(group = "feature", df = df_matrix, averaging_group = "sample_name") {
  feature_names <- names(df)
  if (!("feature" %in% group)) {
    group <- c(group, "feature")
  }

  generate_new_feature_label <- function(feature) {
    feature_split <- strsplit(feature, ";")[[1]]
    num_segments <- length(feature_split)

    idx <- num_segments
    while (
      idx > 1 &&
      (grepl("__$", feature_split[idx]) ||
         grepl("uncultured", feature_split[idx], ignore.case = TRUE) ||
         grepl("unclassified", feature_split[idx], ignore.case = TRUE))
    ) {
      idx <- idx - 1
    }

    start_idx <- idx

    if (idx < num_segments) {
      next_seg <- feature_split[idx + 1]
      if (grepl("uncultured", next_seg, ignore.case = TRUE) ||
          grepl("unclassified", next_seg, ignore.case = TRUE)) {
        start_idx <- idx
        num_segments <- idx + 1
      } else {
        num_segments <- idx
      }
    } else {
      num_segments <- idx
    }

    paste(feature_split[start_idx:num_segments], collapse = ";")
  }

  df_processed <- df %>%
    as.data.frame() %>%
    xPlode_sample_name() %>%
    pivot_longer(cols = all_of(feature_names), names_to = "feature", values_to = "counts")

  if (exists("taxa_levs", inherits = TRUE) && taxa_levs == "ASV") {
    df_processed <- df_processed %>%
      mutate(ASV = feature) %>%
      left_join(ASV_taxa %>% select(ASV, taxa), by = "ASV")
  } else {
    df_processed <- df_processed %>%
      mutate(taxa = feature)
  }

  df_processed %>%
    mutate(
      Domain = sapply(taxa, extract_domain),
      feature_label_original = sapply(taxa, generate_new_feature_label),
      feature_label = if_else(feature_label_original == "__", "unclassified", feature_label_original),
      feature_label = sub("^p__", "", feature_label)
    ) %>%
    group_by(feature) %>%
    mutate(feature_counts = sum(counts)) %>%
    ungroup() %>%
    mutate(total_counts = sum(counts)) %>%
    mutate(total_rel_abun = feature_counts / total_counts) %>%
    group_by(sample_name) %>%
    mutate(sample_counts = sum(counts)) %>%
    ungroup() %>%
    mutate(rel_abun = counts / sample_counts) %>%
    group_by(across(all_of(c(group, averaging_group))), feature_label, feature_label_original, Domain, taxa, total_rel_abun) %>%
    summarize(
      average_rel_abun = mean(rel_abun),
      mean_feature_counts = mean(feature_counts),
      mean_total_counts = mean(total_counts),
      mean_total_rel_abun = mean(total_rel_abun),
      summarized_rows = n(),
      .groups = "drop"
    ) %>%
    group_by(across(all_of(group)), feature_label, feature_label_original, Domain, taxa, total_rel_abun) %>%
    summarize(
      group_mean_rel_abun = mean(average_rel_abun),
      group_mean_feature_counts = mean(mean_feature_counts),
      group_mean_total_counts = mean(mean_total_counts),
      group_mean_total_rel_abun = mean(mean_total_rel_abun),
      mean_prevous_summarized_rows = mean(summarized_rows),
      summarized_rows = n(),
      .groups = "drop"
    ) %>%
    ungroup() %>%
    mutate(
      rel_abun_label = paste0(signif(100 * group_mean_rel_abun, 2), "%"),
      total_rel_abun_label = paste0(signif(100 * total_rel_abun, 2), "%"),
      feature_label2 = paste0("<span style='color:", palette_color[Domain], "'><i>", feature_label, "</i> ", rel_abun_label, "</span>"),
      feature_label3 = paste0(feature_label, " ", signif(100 * group_mean_rel_abun, 2), "%"),
      total_feature_label = paste0("<span style='color:", palette_color[Domain], "'><i>", feature_label, "</i> ", total_rel_abun_label, "</span>")
    ) %>%
    arrange(desc(group_mean_rel_abun)) %>%
    mutate(plot_order_test = factor(feature_label3, levels = rev(unique(feature_label3)))) %>%
    mutate(plot_order = factor(feature_label2, levels = rev(unique(feature_label2))))
}

################################################################################

Kruskal_Mann_Whitney_Test <- function(df = df_matrix, testing_group = "y_axis_group", category_group = "test_across", collapse_data_used = FALSE) {
  dfx <- df %>%
    xPlode_sample_name() %>%
    ungroup()

  kw_group_results <- data.frame()

  unique_categories <- dfx %>%
    distinct(across(all_of(category_group))) %>%
    split(seq(nrow(.)))

  for (category in unique_categories) {
    category_df <- as.data.frame(category)
    category_str <- paste(category_df, collapse = "_")
    category_vals <- as.list(category_df[1, , drop = TRUE])

    keep_rows <- rep(TRUE, nrow(dfx))
    for (i in seq_along(category_group)) {
      keep_rows <- keep_rows & dfx[[category_group[i]]] == category_vals[[i]]
    }

    category_data1 <- dfx[keep_rows, , drop = FALSE]
    rownames(category_data1) <- NULL

    category_columns <- as.list(category_df)
    names(category_columns) <- category_group

    if (nrow(category_data1) < 2) {
      warning(paste("Skipping category", category_str, "due to insufficient rows"))
      kw_group_results <- bind_rows(
        kw_group_results,
        tibble(
          feature = NA,
          estimate = NA,
          statistic = NA,
          p.value = 1,
          testing = testing_group,
          test_type = NA,
          adjusted_pvalue = 1,
          significance = "99",
          adjusted_significance = "99",
          adj_check = NA,
          !!!category_columns
        )
      )
      next
    }

    category_data <- category_data1 %>%
      imPlode_sample_name() %>%
      select(where(~ sum(.) != 0)) %>%
      t() %>%
      as.data.frame() %>%
      mutate(across(everything(), ~ . / sum(.))) %>%
      t() %>%
      as.data.frame() %>%
      xPlode_sample_name()

    check <- category_data %>%
      select(all_of(testing_group)) %>%
      distinct() %>%
      pull()

    cat("Number of unique categories in", testing_group, "for", category_str, ":", length(check), "\n")

    if (length(check) < 2) {
      warning(paste("Skipping category", category_str, "due to insufficient unique categories"))
      kw_group_results <- bind_rows(
        kw_group_results,
        tibble(
          feature = NA,
          estimate = NA,
          statistic = NA,
          p.value = 1,
          testing = testing_group,
          test_type = NA,
          adjusted_pvalue = 1,
          significance = "ns",
          adjusted_significance = "ns",
          adj_check = NA,
          !!!category_columns
        )
      )
      next
    }

    test_type <- if (length(check) == 2) "Welch’s t-test" else "Kruskal-Wallis"
    test_func <- if (length(check) == 2) wilcox.test else kruskal.test

    for (group in testing_group) {
      group_var <- category_data[[group]]

      if (length(unique(group_var)) > 1) {
        feature_data <- category_data %>% select(-any_of(meta_names))

        test_results <- feature_data %>%
          map_df(function(.x) {
            tryCatch(
              broom::tidy(test_func(.x ~ group_var)),
              error = function(e) {
                tibble(
                  estimate = NA,
                  statistic = NA,
                  p.value = 1
                )
              }
            )
          }, .id = "feature") %>%
          mutate(
            data_used = if (collapse_data_used) {
              map_chr(.$feature, function(f) {
                vals <- feature_data[[f]]
                grouped_vals <- split(vals, group_var)
                paste(
                  names(grouped_vals),
                  purrr::map(grouped_vals, ~ paste0(signif(.x, 4), collapse = ", ")),
                  sep = ": ",
                  collapse = " | "
                )
              })
            } else {
              NA_character_
            }
          ) %>%
          mutate(
            testing = group,
            test_type = test_type,
            adjusted_pvalue = p.adjust(p.value, method = "BH"),
            significance = case_when(
              p.value < 0.001 ~ "***",
              p.value < 0.01 ~ "**",
              p.value < 0.05 ~ "*",
              p.value < 0.1 ~ ".",
              TRUE ~ "ns"
            ),
            adjusted_significance = case_when(
              adjusted_pvalue < 0.001 ~ "***",
              adjusted_pvalue < 0.01 ~ "**",
              adjusted_pvalue < 0.05 ~ "*",
              adjusted_pvalue < 0.1 ~ ".",
              TRUE ~ "ns"
            ),
            adj_check = adjusted_pvalue / p.value,
            !!!category_columns
          )

        kw_group_results <- bind_rows(kw_group_results, test_results)
      } else {
        warning(paste("Skipping group", group, "in category", category_str, "because it has only one unique category."))
      }
    }
  }

  kw_group_results
}
