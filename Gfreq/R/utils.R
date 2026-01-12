#' Check that all required arguments are provided
#'
#' This function checks that all elements of a named list are non-NULL.
#' @param args Named list of required arguments to check.
#'   Each element must be non-NULL.
#'
#' @return Nothing. The function stops with an error if at least one element is NULL.
#' @export
check_missing_args <- function(args) {
  null_entries <- names(args)[vapply(args, is.null, logical(1))]

  if (length(null_entries) > 0) {
    msg <- paste(
      "The following required argument is missing:",
      paste(sprintf("\n- %s", null_entries), collapse = "\n")
    )
    stop(msg, call. = FALSE)
  }
}


#' Check file path and stop if it does not exist
#'
#' @param name Name of the file (for display purposes)
#' @param path Path to the file
#' @export
check_file_exists <- function(name, path) {
  if (!file.exists(path)) {
    stop(sprintf("File or folder does not exist:\n- %s: %s", name, path), call. = FALSE)
  }
}

#' Check required and optional input file paths
#'
#' @param required_files Named list of required file/folder paths
#' @param optional_files Named list of optional file/folder paths (default: NULL)
#' @export
check_input_files <- function(required_files, optional_files = NULL) {
  all_files <- required_files
  if (!is.null(optional_files)) {
    optional_files <- optional_files[!vapply(optional_files, is.null, logical(1))]
    all_files <- c(required_files, optional_files)
  }

  for (name in names(all_files)) {
    check_file_exists(name, all_files[[name]])
  }
}


#' Convert VCF genotype to numeric
#'
#' Converts a genotype string (e.g. "0/1" or "0|1") into a numeric representation.
#' @param x A character string representing a genotype in the order REF/ALT (e.g. "0/0", "0/1", "1/1")
#' @return A numeric value: 1 (homozygous ref), 0.5 (heterozygous), 0 (homozygous alt), or NA
#' @export
convert_gt_to_numeric <- function(x) {
  if (x %in% c("0/0", "0|0")) return(1)
  else if (x %in% c("0/1", "1/0", "0|1", "1|0")) return(0.5)
  else if (x %in% c("1/1", "1|1")) return(0)
  else return(NA)
}

#' Extract genotype matrix from VCF and convert to numeric
#'
#' @param vcf_path Path to VCF file
#' @param lib_names_corresp_path Path to a table matching sample library names to sample genotype names (default = NULL)
#' @return Numeric genotyping matrix with values 1 (hom. REF), 0.5 (het), and 0 (hom. ALT) (SNPs x samples)
#' @export
vcf_to_numeric_matrix <- function(vcf_path, lib_names_corresp_path = NULL) {
  writeLines("\n\nImporting vcf genotyping matrix:\n")
  vcf <- read.vcfR(vcf_path)
  gt_matrix <- extract.gt(vcf, element = "GT")
  numeric_matrix <- apply(gt_matrix, c(1, 2), convert_gt_to_numeric)

  if (!is.null(lib_names_corresp_path)) {
    lib_names <- read.table(lib_names_corresp_path, header = TRUE)
    rownames(lib_names) <- lib_names$Library_name
    colnames(numeric_matrix) <- lib_names[colnames(numeric_matrix), "Genotype"]
  }

  return(numeric_matrix)
}

#' Reorient REF/ALT-coded genotype matrix to major/minor per SNP
#'
#' @param numeric_matrix Numeric matrix/data.frame of genotypes (rows = SNPs, cols = samples):
#'   1 = hom. REF, 0.5 = het, 0 = hom. ALT.
#'
#' @return Numeric matrix with rows reoriented to major/minor:
#'   1 = hom. major, 0.5 = het, 0 = hom. minor.
#' @export
num_gt_matrix_to_majmin <- function(numeric_matrix) {
  f_ref <- rowMeans(numeric_matrix, na.rm = TRUE)
  inv <- !is.na(f_ref) & (f_ref < 0.5)
  numeric_matrix[inv, ] <- 1 - numeric_matrix[inv, ]
  numeric_matrix
}

#' Build a major/minor-oriented numeric genotype matrix from a VCF
#'
#' @param vcf_path Path to VCF file.
#' @param lib_names_corresp_path Optional path to a two-column mapping table
#'   (Library_name -> Genotype) for renaming samples.
#'
#' @return Numeric genotyping matrix (SNPs x samples) oriented to major/minor:
#'   1 = hom. major, 0.5 = het, 0 = hom. minor.
#' @export
vcf_to_majmin_num_gt_matrix <- function(vcf_path, lib_names_corresp_path = NULL) {
  numeric_matrix <- vcf_to_numeric_matrix(
    vcf_path = vcf_path,
    lib_names_corresp_path = lib_names_corresp_path
  )
  majmin_num_matrix<-num_gt_matrix_to_majmin(numeric_matrix)
  return(majmin_num_matrix)
}

#' Check that required columns are present in a data frame
#'
#' Verifies that all expected column names are present in the input data frame.
#'
#' @param df A data frame to check.
#' @param required_cols A character vector of required column names.
#'
#' @return Nothing.
#' @export
check_columns <- function(df, required_cols) {
  df_name <- deparse(substitute(df))
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf("Missing required column(s) in '%s': %s",
                 df_name, paste(missing, collapse = ", ")))
  }
}

#' Set row names from CHROM and POS columns and remove them
#'
#' Combines the `CHROM` and `POS` columns to form row names of the form `CHROM_POS`,
#' and removes those two columns from the data frame.
#'
#' @param df A data frame with columns named `CHROM` and `POS`.
#'
#' @return A data frame with updated row names and without the `CHROM` and `POS` columns.
#' @export
set_rownames_from_chrom_pos <- function(df) {
  check_columns(df, c("CHROM", "POS"))

  rownames(df) <- paste(df$CHROM, df$POS, sep = "_")
  df <- df[, !(colnames(df) %in% c("CHROM", "POS"))]

  return(df)
}


#' Melt a Genotype Frequency Table to Long Format
#'
#' @param genotype_freqs Dataframe containing genotype frequencies, with a "Mixture" column and at least two other columns representing components.
#' @param freq_col_name A character string specifying the name to use for the frequency column in the output.
#'
#' @return A dataframe in long format, with columns Mixture, Component, and a frequency column named according to "freq_col_name".
#' @export
melt_genotype_freqs<-function(genotype_freqs, freq_col_name){
  genotype_freqs_melt<-melt(genotype_freqs, id.vars = "Mixture")
  colnames(genotype_freqs_melt)[2:3]<-c("Component", freq_col_name)
  return(genotype_freqs_melt)
}

#' Filter SNPs based on minimum sequencing depth
#'
#' This function removes SNPs (rows) whose minimum sequencing depth across samples
#' is below a specified threshold (`min_depth`). It applies the same filtering
#' to both the main depth matrix and any additional data frames provided.
#'
#' @param depths Read depths data.frame (SNPs x sample). Row names must correspond to SNP identifiers.
#' @param extra_files Either a single data.frame or a list of data.frames containing
#'   additional SNP-related information (e.g. allele frequencies) with the same
#'   row names as `depths`.
#' @param min_depth Minimum read depth threshold (numeric value). SNPs with at least one sample showing a depth < threshold in the depths data.frame we be removed.
#'
#' @return A list of filtered data.frames, starting with the filtered `depths`
#'   followed by the filtered elements of `extra_files`, in the same order.
#' @export
filter_lowdepth_snps<-function(depths, extra_files, min_depth){
  if (is.data.frame(extra_files)) {
    extra_files <- list(extra_files)
  }

  snp_to_filter<-rownames(depths)[apply(depths, 1, min) < min_depth]
  n_filtered <- length(snp_to_filter)

  filtered_dfs <- lapply(c(list(depths), extra_files), function(df) {
    df[!rownames(df) %in% snp_to_filter, ]
  })

  if (n_filtered == 0) {
    writeLines(sprintf("No SNPs filtered (all depths ≥ %d).", min_depth))
  } else {
    writeLines(sprintf("\n\nINFO: Filtered %d SNPs with min depth < %d.", n_filtered, min_depth))
  }

  return(filtered_dfs)
}

#' Check that multiple SNP tables share identical SNP sets
#'
#' @param dfs [list] A list of data frames or matrices to check. Each element
#'   must have rownames corresponding to SNP identifiers (e.g.: all_freqs or read_depths dataframes)
#' @param context [character] Optional string describing the context of the check
#'
#' @return Invisibly returns `TRUE` if all tables contain identical SNP sets.
#'   Stops execution with an error message if any mismatch is detected.
#' @export
check_identical_snps <- function(dfs, context = NULL) {
  rn_list <- lapply(dfs, rownames)

  ref <- rn_list[[1]]
  for (i in seq_along(rn_list)[-1]) {
    if (!setequal(ref, rn_list[[i]])) {
      stop(paste0(
        "Mismatch between SNPs detected among input files",
        if (!is.null(context)) paste0(" (", context, ")"), "."
      ))
    }
  }

  invisible(TRUE)
}

#' Load and format all input data for genotype frequency estimation.
#'
#' This function imports, formats, and filters input data for genotype frequency
#' estimation from mixtures (and optionally from component libraries).
#' It reads the genotyping matrix from a VCF file, loads allele frequency and
#' depth tables, optionally processes expected frequencies, and applies SNP
#' filtering based on a minimum depth threshold.
#'
#' The function automatically adapts to the presence or absence of component
#' input files and expected frequency tables. It returns a named list containing
#' the formatted data, ready to be passed to downstream estimation functions.
#'
#' @param genotyping_vcf_path [character] Path to the VCF file used to build the genotyping matrix.
#' @param lib_names_corresp_path [character] Path to the correspondence table mapping library names to samples.
#' @param allele_freqs_mixtures_path [character] Path to the allele frequency table for mixtures.
#' @param snp_depths_mixtures_path [character] Path to the SNP depth table for mixtures.
#' @param expected_freqs_mixtures_path [character|NULL] Optional path to the expected genotype frequencies for mixtures.
#' @param allele_freqs_components_path [character|NULL] Optional path to the allele frequency table for component libraries.
#' @param snp_depths_components_path [character|NULL] Optional path to the SNP depth table for component libraries.
#' @param expected_freqs_components_path [character|NULL] Optional path to the expected genotype frequencies for component libraries.
#' @param min_depth [numeric] Minimum SNP depth threshold used for filtering.
#'
#' @return A named list containing the following elements:
#'   - `genotyping_matrix`: numeric matrix built from the VCF file.
#'   - `allele_freqs_mixtures`: formatted allele frequency table for mixtures.
#'   - `snp_depths_mixtures`: formatted SNP depth table for mixtures.
#'   - `expected_freqs_mixtures_melt`: melted expected genotype frequencies for mixtures (if provided).
#'   - `allele_freqs_components`: formatted allele frequency table for components (if provided).
#'   - `snp_depths_components`: formatted SNP depth table for components (if provided).
#'   - `expected_freqs_components_melt`: melted expected genotype frequencies for components (if provided).
#' @export
load_inputs <- function(
    genotyping_vcf_path,
    lib_names_corresp_path,
    allele_freqs_mixtures_path,
    snp_depths_mixtures_path,
    expected_freqs_mixtures_path = NULL,
    allele_freqs_components_path = NULL,
    snp_depths_components_path = NULL,
    expected_freqs_components_path = NULL,
    min_depth = 0
) {
  genotyping_matrix <- vcf_to_numeric_matrix(
    vcf_path = genotyping_vcf_path,
    lib_names_corresp_path = lib_names_corresp_path
  )

  allele_freqs_mixtures <- read.table(allele_freqs_mixtures_path, header = TRUE)
  snp_depths_mixtures <- read.table(snp_depths_mixtures_path, header = TRUE)

  allele_freqs_mixtures <- set_rownames_from_chrom_pos(allele_freqs_mixtures)
  snp_depths_mixtures <- set_rownames_from_chrom_pos(snp_depths_mixtures)

  check_identical_snps(dfs=list(allele_freqs_mixtures, snp_depths_mixtures), context = "between mixtures all_freqs and read_depths")

  expected_freqs_mixtures_melt <- NULL
  if (!is.null(expected_freqs_mixtures_path)) {
    expected_freqs_mixtures <- read.table(expected_freqs_mixtures_path, header = TRUE)
    expected_freqs_mixtures_melt <- melt_genotype_freqs(expected_freqs_mixtures, "ExpFreq")
  }

  allele_freqs_components <- NULL
  snp_depths_components <- NULL
  expected_freqs_components_melt <- NULL

  if (!is.null(allele_freqs_components_path)) {
    allele_freqs_components <- read.table(allele_freqs_components_path, header = TRUE)
    snp_depths_components <- read.table(snp_depths_components_path, header = TRUE)

    allele_freqs_components <- set_rownames_from_chrom_pos(allele_freqs_components)
    snp_depths_components <- set_rownames_from_chrom_pos(snp_depths_components)

    check_identical_snps(dfs=list(allele_freqs_components, snp_depths_components), context = "between components all_freqs and read_depths")
    check_identical_snps(dfs=list(allele_freqs_components, allele_freqs_mixtures), context = "between components and mixtures all_freqs")


    if (!is.null(expected_freqs_components_path)) {
      expected_freqs_components <- read.table(expected_freqs_components_path, header = TRUE)
      expected_freqs_components_melt <- melt_genotype_freqs(expected_freqs_components, "ExpFreq")
    }

    c(
      snp_depths_mixtures,
      allele_freqs_mixtures,
      snp_depths_components,
      allele_freqs_components
    ) %<-% filter_lowdepth_snps(
      depths = snp_depths_mixtures,
      extra_files = list(
        allele_freqs_mixtures,
        snp_depths_components,
        allele_freqs_components
      ),
      min_depth = min_depth
    )
  } else {
    c(snp_depths_mixtures, allele_freqs_mixtures) %<-% filter_lowdepth_snps(
      depths = snp_depths_mixtures,
      extra_files = allele_freqs_mixtures,
      min_depth = min_depth
    )
  }

  return(list(
    genotyping_matrix = genotyping_matrix,
    allele_freqs_mixtures = allele_freqs_mixtures,
    snp_depths_mixtures = snp_depths_mixtures,
    expected_freqs_mixtures_melt = expected_freqs_mixtures_melt,
    allele_freqs_components = allele_freqs_components,
    snp_depths_components = snp_depths_components,
    expected_freqs_components_melt = expected_freqs_components_melt
  ))
}


# ----

#' Estimate genotype frequencies from genotyping matrix and allele frequencies
#'
#' @param genotyping_matrix Genotyping matrix with values 1 (hom. REF), 0.5 (het), and 0 (hom. ALT) (SNPs x components)
#' @param allele_freqs Allele frequency matrix (SNPs x mixtures)
#' @param snp_depths Read depths matrix (SNPs x mixtures) (default = NULL)
#' @param expected_freqs_melt Melted dataframe of expected genotype frequencies (default = NULL)
#' @param est_freq_col_name Column name for estimated frequencies (default = "EstFreq")
#' @param out_file Path to the output file to write the estimated genotype frequencies (default = NULL)
#' @return A melted dataframe of the estimated genotype frequencies
#' @export
estimate_genotype_freqs <- function(genotyping_matrix, allele_freqs, snp_depths = NULL, expected_freqs_melt = NULL, est_freq_col_name = "EstFreq", out_file = NULL) {
  nb_mixtures <- ncol(allele_freqs)
  nb_components <- ncol(genotyping_matrix)
  genotype_freqs <- data.frame(matrix(nrow = nb_components, ncol = nb_mixtures))
  colnames(genotype_freqs) <- colnames(allele_freqs)
  rownames(genotype_freqs) <- colnames(genotyping_matrix)

  if (!is.null(snp_depths) && !setequal(rownames(allele_freqs), rownames(snp_depths))) {
    stop("Mismatch between SNPs in allele_freqs and snp_depths - they must be identical.")
  }

  common_snps <- intersect(rownames(genotyping_matrix), rownames(allele_freqs))

  message(sprintf("Estimating genotype frequencies with %d SNPs.", length(common_snps)))
  
  allele_freqs <- allele_freqs[common_snps,]
  genotyping_matrix <- genotyping_matrix[common_snps,]
  snp_depths <- if (!is.null(snp_depths)) snp_depths[common_snps,] else NULL

  for (mixture_name in colnames(allele_freqs)) {
    weights <- if (!is.null(snp_depths)) snp_depths[, mixture_name] else NULL
    est_freqs <- coefficients(lm(allele_freqs[, mixture_name] ~ genotyping_matrix - 1, weights = weights))
    genotype_freqs[, mixture_name] <- round(est_freqs, 4)
  }

  genotype_freqs$Component <- rownames(genotype_freqs)
  genotype_freqs_melt <- melt(genotype_freqs, id.vars = "Component")
  colnames(genotype_freqs_melt)[2:3] <- c("Mixture", est_freq_col_name)

  if (!is.null(expected_freqs_melt)) {
    genotype_freqs_melt <- merge(genotype_freqs_melt, expected_freqs_melt, by = c("Component", "Mixture"))
  }

  if (!is.null(out_file)){
    write.table(genotype_freqs_melt, file = out_file, row.names = FALSE, quote = F, sep ="\t")
  }

  return(genotype_freqs_melt)
}



# -----

#' Force alphabetical ordering of a column based on labels
#'
#' @param x A vector (character or factor).
#'
#' @return A factor with levels sorted alphabetically.
#' @export
factor_alphabetical <- function(x) {
  x_chr <- as.character(x)
  factor(x_chr, levels = sort(unique(x_chr)))
}

#' Convert a vector to a factor with an explicit level order
#'
#' @param x A vector (character or factor).
#' @param levels Character vector defining the desired order.
#'
#' @return A factor with levels ordered as provided.
#' @export
factor_with_levels <- function(x, levels) {
  x_chr <- as.character(x)

  missing_levels <- setdiff(unique(x_chr), levels)
  if (length(missing_levels) > 0) {
    stop(
      sprintf(
        "Provided levels are missing value(s): %s",
        paste(missing_levels, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  factor(x_chr, levels = levels)
}

#' Plot boxplots for a long-format dataset
#'
#' @param df Long-format dataframe.
#' @param x_col String. Column name used for the x-axis (grouping variable).
#' @param y_col String. Column name used for the y-axis (numeric variable).
#' @param x_lab String. Label for the x-axis.
#' @param y_lab String. Label for the y-axis.
#' @param title String. Plot title.
#' @param subtitle Optional string. Plot subtitle (default = NULL).
#' @param hline Optional numeric value for a horizontal reference line (default = NULL).
#' @param show_x_labels Logical. Whether to display x-axis labels (default = TRUE).
#' @param out_file Optional path to save the plot as a PNG. If NULL, the plot is printed.
#'
#' @return A ggplot object (printed or saved).
#' @export
plot_boxplots <- function(
    df,
    x_col,
    y_col,
    x_lab,
    y_lab,
    title,
    subtitle = NULL,
    hline = NULL,
    show_x_labels = TRUE,
    out_file = NULL
) {
  
  check_columns(df, c(x_col, y_col))
  
  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_boxplot(outlier.size = 0.5, color = "grey30", fill = "skyblue") +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid.major.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      x = x_lab,
      y = y_lab,
      title = title,
      subtitle = subtitle
    )
  
  if (!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = "dashed", color = "orange") +
      annotate("text", x = -Inf, y = hline, label = hline,
               hjust = 1.45, vjust = +0.4, color="orange", cex=3.2) +
      coord_cartesian(clip="off")
  }
  
  if (!isTRUE(show_x_labels)) {
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  } else {
    p <- p +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
      )
  }
  
  if (!is.null(out_file)) {
    ggsave(out_file, p, width = 10, height = 8, dpi = 300, bg = "white")
    writeLines(c("", "Output file saved in:", out_file, ""))
  } else {
    print(p)
  }
}

#' Plot boxplots across varying conditions
#'
#' @param df_melt Data frame in long format
#' @param x_col Column name for x-axis (string)
#' @param y_col Column name for y-axis (string)
#' @param fill_col Column name for fill (string)
#' @param ref_col Optional column to use for horizontal dashed lines (default = NULL)
#' @param y_lab Label for y-axis
#' @param x_lab Label for x-axis
#' @param fill_lab Legend title for fill
#' @param show_zero Logical (default = FALSE). If TRUE, adds a horizontal dashed line at y = 0 to help visualize deviations around zero.
#'
#' @return A ggplot object
#' @export
plot_condition_effect_boxplots <- function(df_melt, x_col, y_col, fill_col, ref_col = NULL, y_lab = "", x_lab = "", fill_lab = "Group", show_zero = F) {
  # Order legend and conditions
  df_melt[[fill_col]] <- factor_alphabetical(df_melt[[fill_col]])

  condition_levels <- attr(df_melt, "condition_levels")

  if (!is.null(condition_levels)) {
    df_melt[[x_col]] <- factor_with_levels(df_melt[[x_col]], levels = condition_levels)
  } else {
      df_melt[[x_col]] <- factor_alphabetical(df_melt[[x_col]])
  }

  # Define palette
  fill_levels <- levels(df_melt[[fill_col]])
  color_palette <- scales::hue_pal()(length(fill_levels))
  names(color_palette) <- fill_levels

  p <- ggplot(df_melt, aes(x = .data[[x_col]], y = .data[[y_col]], fill = .data[[fill_col]]))+
    scale_x_discrete(expand = expansion(mult = c(0, 0)))


  # Optional reference lines and labels
  if (!is.null(ref_col)) {
    ref_vals <- as.numeric(as.character(df_melt[[ref_col]]))
    p <- p +
      geom_hline(
        aes(
          yintercept = ref_vals,
          color = .data[[fill_col]]),
          linetype = "dashed"
        ) +
      geom_text(
        aes(
          x = 0,
          y = ref_vals,
          label = .data[[ref_col]],
          color = .data[[fill_col]]
        ),
        hjust = 0,
        vjust = -0.5,
        size = 3
      ) +
      scale_color_manual(values = color_palette, guide = "none")
  }
  if (isTRUE(show_zero)) {
    p <- p +
      geom_hline(yintercept = 0, color = "grey40", linetype = "dashed", linewidth = 0.6)
  }

  p <- p +
    geom_boxplot() +
      scale_fill_manual(values = color_palette) +
      xlab(x_lab) +
      ylab(y_lab) +
      labs(fill = fill_lab) +
      theme_minimal(base_size = 16)

  return(p)
}


#' Plot boxplots of (1 or more sets of) estimated frequencies by expected frequencies
#'
#' @param freqs_df Dataframe with at least one estimated frequency column and a "ExpFreq" column (providing the expected frequencies)
#' @param variable_name Label for x-axis (below all conditions)
#' @param out_file Path to the output file to save the plot (default = NULL)
#' @return A ggplot object
#' @export
compare_boxplots_est_freq <- function(freqs_df, variable_name, out_file = NULL) {
  df_melt <- melt(freqs_df[, -c(1:2)], id.vars = "ExpFreq")
  colnames(df_melt)[2:3] <- c("Condition", "EstFreq")
  attr(df_melt, "condition_levels") <- attr(freqs_df, "condition_levels")


  p <- plot_condition_effect_boxplots(
    df_melt   = df_melt,
    x_col     = "Condition",
    y_col     = "EstFreq",
    fill_col  = "ExpFreq",
    ref_col   = "ExpFreq",
    y_lab     = "Estimated frequency",
    x_lab     = variable_name,
    fill_lab  = "Expected frequency"
  )

  if (!is.null(out_file)) {
    n_conditions <- length(unique(df_melt$Condition))
    n_expFreq <- length(unique(df_melt$ExpFreq))
    ggsave(out_file, plot = p, width = min(49, 2+0.5*n_expFreq*n_conditions), height = 4, bg = "white")
  }

  return(p)
}



#' Compute mean and SD per expected frequency
#'
#' @param freqs_df Dataframe with an "ExpFreq" colomn (providing the expected frequencies) and 1 or more estimated frequencies column(s)
#' @param out_dir Path the the output directory to write output files (default = NULL)
#' @param suffix Suffix to add in the output file names (default = "")
#' @return A list with two dataframes: expected frequencies mean and standard deviation
#' @export
compute_stats_per_exp_freq <- function(freqs_df, out_dir = NULL, suffix = "") {
  output <- list()
  df_melt <- melt(freqs_df[, -c(1:2)], id.vars = "ExpFreq")
  mean_df <- aggregate(value ~ ExpFreq + variable, data = df_melt, function(x) round(mean(x), 4))
  sd_df   <- aggregate(value ~ ExpFreq + variable, data = df_melt, function(x) round(sd(x), 4))

  output$mean <- dcast(mean_df, ExpFreq ~ variable)
  output$sd <- dcast(sd_df, ExpFreq ~ variable)

  if (! is.null(out_dir)){
    names(output$mean)[-1] <- paste0(names(output$mean)[-1], "_mean")
    names(output$sd)[-1] <- paste0(names(output$sd)[-1], "_SD")
    write.table(output$mean, glue("{out_dir}/est_geno_freqs_mean{suffix}.tsv"), row.names = FALSE, quote = F, sep ="\t")
    write.table(output$sd, glue("{out_dir}/est_geno_freqs_sd{suffix}.tsv"), row.names = FALSE, quote = F, sep ="\t")
  }

  return(output)
}

#' Add linear model annotation to an existing plot
#'
#' @param fit A fitted lm object
#' @param col Text color (default = "darkred")
#' @param cex Text size (default = 1)
#' @param x_offset Relative x offset from the left (default = 0.05)
#' @param y_offset Relative y offset from the top (default = 0.05)
#'
#' @return NULL
add_lm_annotation <- function(
  fit,
  col = "darkred",
  cex = 1,
  x_offset = 0.05,
  y_offset = 0.05
) {
  stopifnot(inherits(fit, "lm"))

  coeffs <- coef(fit)
  r2 <- summary(fit)$r.squared

  usr <- par("usr")
  x_pos <- usr[1] + x_offset * (usr[2] - usr[1])
  y_pos <- usr[4] - y_offset * (usr[4] - usr[3])

  text(
    x_pos,
    y_pos,
    labels = bquote(
      y == .(round(coeffs[1], 3)) + .(round(coeffs[2], 3)) * x
    ),
    adj = c(0, 1),
    col = col,
    cex = cex
  )

  text(
    x_pos,
    y_pos - 0.08 * (usr[4] - usr[3]),
    labels = bquote(
      R^2 == .(round(r2, 3))
    ),
    adj = c(0, 1),
    col = col,
    cex = cex
  )

  invisible(NULL)
}


#' Plot estimated vs expected frequencies with regression
#'
#' @param freqs_df Dataframe with an "ExpFreq" colomn (providing the expected frequencies) and 1 or more estimated frequencies column(s)
#' @param extra_freqs_df Dataframe with additional points to plot (they won't be used to compute the regression) (default = NULL)
#' @param out_file Path to the output file to save the plot (default = NULL)
#' @export
plot_correlation_with_exp_freq <- function(freqs_df, extra_freqs_df = NULL, out_file = NULL) {

  if (! is.null(out_file)) { pdf(out_file, width = 11.11, height = 8.33) }

  for (condition in colnames(freqs_df)[!colnames(freqs_df) %in% c("Component", "Mixture", "ExpFreq")]) {
    par(
      cex.axis = 1.2,
      cex.lab = 1.4,
      cex.main = 1.5
    )
    plot(freqs_df[, condition] ~ freqs_df$ExpFreq, pch = 16, las = 2, col = alpha("darkred", 0.3),
         ylim = c(0, 1.1), xlim = c(0, 1.1), ylab = "Estimated frequencies", xlab = "Expected frequencies", main = condition)

    if (!is.null(extra_freqs_df)) {
      points(extra_freqs_df$ExpFreq, extra_freqs_df$EstFreq, pch = 17, col = alpha("black", 0.5))
    }

    abline(a = 0, b = 1, col = "black")
    fit <- lm(freqs_df[, condition] ~ freqs_df$ExpFreq)
    coeffs <- coef(fit)

    if (length(unique(freqs_df$ExpFreq)) > 1) {
      abline(fit, col = "darkred")
      add_lm_annotation(fit)
    }
  }

  if (! is.null(out_file)) { dev.off() }
}

#' Compute mean bias and significance groups for one condition
#'
#' @param freqs_df A data.frame with columns Component, ExpFreq, and an error column
#' @param error_col String. Name of the error column to analyze (must start with "error_")
#' @param condition String. Condition name used to label significance group column
#'
#' @return A list with:
#'   - components: data.frame with mean bias per Component and significance group
#'   - exp_freqs: data.frame with mean bias per ExpFreq
#' @export
compute_bias_from_errors <- function(freqs_df, error_col, condition) {
  stopifnot(startsWith(error_col, "error_"))

  formula_comp <- as.formula(glue("{error_col} ~ Component"))
  formula_exp  <- as.formula(glue("{error_col} ~ ExpFreq"))

  comp_mean <- aggregate(formula_comp, freqs_df, function(x) round(mean(x), 4))
  expf_mean <- aggregate(formula_exp, freqs_df, function(x) round(mean(x), 4))

  aov_res <- aov(formula_comp, freqs_df)
  tukey_res <- TukeyHSD(aov_res)
  sig_groups <- multcompLetters4(aov_res, tukey_res)[[1]]$Letters

  bias_col <- sub("^error_", "bias_", error_col)
  colnames(comp_mean)[2] <- bias_col
  colnames(expf_mean)[2] <- bias_col

  comp_mean[[glue("{condition}_significance_group")]] <- sig_groups[comp_mean$Component]

  return(list(
    components = comp_mean,
    exp_freqs = expf_mean
  ))
}


#' Plot grouped error boxplots by component and expected frequency
#'
#' @param errors_long_df Long-format dataframe with columns: Component, ExpFreq, Error, Condition
#' @param out_dir Optional output directory to save plots (default = NULL)
#' @param suffix Optional suffix to append to output file names (default = "")
#'
#' @return NULL. Plots are printed if no output directory is given.
#' @export
plot_error_boxplots <- function(
    errors_long_df,
    variable_name,
    out_dir = NULL,
    suffix = ""
) {
  p_component <- plot_condition_effect_boxplots(
    df_melt = errors_long_df,
    x_col = "Condition",
    y_col = "Error",
    fill_col = "Component",
    y_lab = "Error (estimated - expected)",
    x_lab = variable_name,
    fill_lab = "Component",
    show_zero = T
  )

  p_expfreq <- plot_condition_effect_boxplots(
    df_melt = errors_long_df,
    x_col = "Condition",
    y_col = "Error",
    fill_col = "ExpFreq",
    y_lab = "Error (estimated - expected)",
    x_lab = variable_name,
    fill_lab = "Expected frequency",
    show_zero = T
  )

  if (!is.null(out_dir)) {
    n_conditions <- length(unique(errors_long_df$Condition))
    n_components <- length(unique(errors_long_df$Component))
    n_expFreq <- length(unique(errors_long_df$ExpFreq))
    ggsave(glue("{out_dir}/error_by_component{suffix}.png"), p_component, width = min(49, 3+0.5*n_components*n_conditions), height = 4, dpi = 300, bg = "white")
    ggsave(glue("{out_dir}/error_by_exp_freq{suffix}.png"), p_expfreq, width = min(49, 3+0.5*n_expFreq*n_conditions), height = 4, dpi = 300, bg = "white")
  } else {
    print(p_component)
    print(p_expfreq)
  }
}


#' Compute estimation bias and plot grouped error boxplots
#'
#' @param freqs_df Dataframe with an "ExpFreq" colomn (providing the expected frequencies), a "Component" column, and 1 or more estimated frequencies column(s)
#' @param variable_name Label for the x-axis (e.g. name of the strategy variable)
#' @param out_dir Path the the output directory to write output files (default = NULL)
#' @param suffix Suffix to add in the output file names (default = "")
#' @return A list with bias dataframes per component and per expected frequency
#' @export
compute_bias <- function(freqs_df, variable_name = "", out_dir = NULL, suffix = "") {
  est_freq_cols <- setdiff(colnames(freqs_df), c("Component", "Mixture", "ExpFreq"))
  components_bias <- data.frame(Component = unique(freqs_df$Component))
  exp_freqs_bias <- data.frame(ExpFreq = unique(freqs_df$ExpFreq))
  errors_df<-freqs_df
  for (condition in est_freq_cols) {
    error_col <- glue("error_{condition}")

    errors_df[[error_col]] <- errors_df[[condition]] - errors_df[["ExpFreq"]]

    bias <- compute_bias_from_errors(errors_df, error_col, condition)
    components_bias <- merge(components_bias, bias$components)
    exp_freqs_bias  <- merge(exp_freqs_bias, bias$exp_freqs)

  }

  errors_long_df <- errors_df %>%
    tidyr::pivot_longer(
      cols = tidyselect::all_of(paste0("error_", est_freq_cols)),
      names_to = "Condition",
      names_prefix = "error_",
      values_to = "Error"
    )

  attr(errors_long_df, "condition_levels") <- attr(freqs_df, "condition_levels")

  plot_error_boxplots(errors_long_df, variable_name, out_dir, suffix)

  if (! is.null(out_dir)) {
    write.table(components_bias, glue("{out_dir}/bias_per_component{suffix}.tsv"), row.names = FALSE, quote = F, sep ="\t")
    write.table(exp_freqs_bias, glue("{out_dir}/bias_per_exp_freq{suffix}.tsv"), row.names = FALSE, quote = F, sep ="\t")
  }
  return(list(components = components_bias, exp_freqs = exp_freqs_bias))
}


#' Compare multiple sets of estimated genotype frequencies
#'
#' This function compares estimated genotype frequencies from different strategies by:
#' - plotting estimated vs. expected frequencies as boxplots
#' - computing summary statistics per expected frequency
#' - plotting correlation with expected frequencies
#' - computing estimation bias
#'
#' @param freqs_to_compare Dataframe with a column "ExpFreq" (providing the expected frequencies) and one column per strategy providing the estimated frequencies
#' @param variable_name Label for the x-axis (e.g. name of the strategy variable)
#' @param suffix Suffix to add in the output file names
#' @param out_dir Path to the output directory to write result files (plots and tables)
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{stats}{A list of dataframes from \code{compute_stats_per_exp_freq()}: mean and standard deviation per expected frequency}
#'   \item{bias}{A list of dataframes from \code{compute_bias()}: bias per component and per expected frequency}
#' }
#' @export
compare_different_est_freqs<-function(freqs_to_compare, variable_name, out_dir, suffix){
  compare_boxplots_est_freq(
    freqs_df = freqs_to_compare,
    variable_name = variable_name,
    out_file = glue("{out_dir}/est_geno_freqs_boxplot{suffix}.png")
  )

  stats <- compute_stats_per_exp_freq(
    freqs_df = freqs_to_compare,
    out_dir = out_dir,
    suffix = suffix
  )

  plot_correlation_with_exp_freq(
    freqs_df = freqs_to_compare,
    out_file = glue("{out_dir}/est_exp_freqs_scatterplot{suffix}.pdf")
  )

  bias <- compute_bias(
    freqs_df = freqs_to_compare,
    variable_name = variable_name,
    out_dir = out_dir,
    suffix = suffix
  )

  return(list(stats = stats, bias = bias))

}


#' Evaluate the effect of weight vector on genotype frequency estimation
#'
#' This function compares estimated genotype frequencies using:
#' - read depth as weight
#' - no weighting
#' Then compares both results with plots, bias, and summary stats.
#'
#' @param genotyping_matrix Genotyping matrix with values 1 (hom. REF), 0.5 (het), and 0 (hom. ALT) (SNPs x components)
#' @param allele_freqs Allele frequency matrix (SNPs x mixtures)
#' @param snp_depths Read depths matrix (SNPs x mixtures)
#' @param expected_freqs_melt Melted dataframe of expected genotype frequencies
#' @param out_dir Path the the output directory to write output files
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{stats}{A list of dataframes from \code{compute_stats_per_exp_freq()}: mean and standard deviation per expected frequency}
#'   \item{bias}{A list of dataframes from \code{compute_bias()}: bias per component and per expected frequency}
#' }
#' @export
evaluate_weight_vector_effect<-function(genotyping_matrix, allele_freqs, snp_depths, expected_freqs_melt, out_dir){
  writeLines("\n\nCompare the effect of the weight vector on genotype frequency estimation:\n")

  freqs_depth_weight<-estimate_genotype_freqs(
    genotyping_matrix = genotyping_matrix,
    allele_freqs = allele_freqs,
    snp_depths = snp_depths,
    expected_freqs_melt = expected_freqs_melt,
    est_freq_col_name = "read_depth"
  )
  freqs_no_weight<-estimate_genotype_freqs(
    genotyping_matrix = genotyping_matrix,
    allele_freqs = allele_freqs,
    expected_freqs_melt = expected_freqs_melt,
    est_freq_col_name = "none"
  )
  freqs_to_compare<-merge(freqs_depth_weight, freqs_no_weight)

  attr(freqs_to_compare, "condition_levels") <- c("none", "read_depth")

  stats<-compare_different_est_freqs(
    freqs_to_compare = freqs_to_compare,
    variable_name = "Weight vector",
    out_dir = out_dir,
    suffix = "_weight_effect"
  )
  return(stats)
}


#' Compare the effect of random SNP subsampling on genotype frequency estimation
#'
#' @param genotyping_matrix Genotyping matrix with values 1 (hom. REF), 0.5 (het), and 0 (hom. ALT) (SNPs x components)
#' @param allele_freqs Allele frequency matrix (SNPs x mixtures)
#' @param step_size Number of SNPs to increment by when subsampling
#' @param snp_depths Read depths matrix (SNPs x mixtures) (default = NULL)
#' @param expected_freqs_melt Melted dataframe of expected genotype frequencies
#' @return Dataframe of estimated frequencies per condition
#' @export
estimate_geno_freqs_snps_subsampling <- function(genotyping_matrix, allele_freqs, step_size, snp_depths = NULL, expected_freqs_melt) {

  common_snps <- intersect(rownames(genotyping_matrix), rownames(allele_freqs))
  total_nb_snps <- length(common_snps)
  nb_snps_to_sample <- unique(c(
    seq(from = max(5, step_size), to = total_nb_snps, by = step_size),
    total_nb_snps
  ))

  condition_levels <- paste0(nb_snps_to_sample, "_SNPs")

  subsampled_genotype_freqs <- data.frame(Mixture = colnames(allele_freqs))

  for (nb_snps in nb_snps_to_sample) {
    condition <- glue("{nb_snps}_SNPs")
    snps <- sample(common_snps, nb_snps)
    if (! is.null(snp_depths)){
      snp_depths_subset<-snp_depths[snps,]
    } else {
      snp_depths_subset<-NULL
    }

    genotype_freqs <- estimate_genotype_freqs(
      genotyping_matrix = genotyping_matrix[snps,],
      allele_freqs = allele_freqs[snps,],
      snp_depths = snp_depths_subset,
      expected_freqs_melt = expected_freqs_melt,
      est_freq_col_name = condition
      )

    subsampled_genotype_freqs <- merge(subsampled_genotype_freqs, genotype_freqs)
  }

  attr(subsampled_genotype_freqs, "condition_levels") <- condition_levels

  return(subsampled_genotype_freqs)
}

#' Estimate genotype frequency errors under SNP subsampling with replicates
#' @param genotyping_matrix Genotyping matrix with values 1 (hom. REF), 0.5 (het), and 0 (hom. ALT) (SNPs x components)
#' @param allele_freqs Allele frequency matrix (SNPs x mixtures)
#' @param step_size Number of SNPs to increment by when subsampling
#' @param snp_depths Read depths matrix (SNPs x mixtures) (default = NULL)
#' @param nb_reps Integer. Number of random subsampling replicates per SNP subset size (default = 5).
#' @param expected_freqs_melt Melted dataframe of expected genotype frequencies
#' @return A dataframe of the estimation errors
#' @export
estimate_freq_errors_subsampling <- function(genotyping_matrix, allele_freqs, step_size, nb_reps=5, snp_depths = NULL, expected_freqs_melt) {
  
  common_snps <- intersect(rownames(genotyping_matrix), rownames(allele_freqs))
  total_nb_snps <- length(common_snps)
  nb_snps_to_sample <- unique(c(
    seq(from = max(5, step_size), to = total_nb_snps, by = step_size),
    total_nb_snps
  ))
  
  res <- vector("list", length(nb_snps_to_sample) * nb_reps)
  k <- 0L
  
  for (nb_snps in nb_snps_to_sample) {
    for (rep in 1:nb_reps){
      k <- k + 1L
      snps <- sample(common_snps, nb_snps)
      if (! is.null(snp_depths)){
        snp_depths_subset<-snp_depths[snps,]
      } else {
        snp_depths_subset<-NULL
      }
      
      genotype_freqs <- suppressMessages(
        estimate_genotype_freqs(
          genotyping_matrix = genotyping_matrix[snps,],
          allele_freqs = allele_freqs[snps,],
          snp_depths = snp_depths_subset,
          expected_freqs_melt = expected_freqs_melt
        )
      )
      errors <- genotype_freqs$EstFreq - genotype_freqs$ExpFreq
      
      res[[k]] <- data.frame(
        nb_SNPs = nb_snps,
        replicate = rep,
        error = errors
      )
    }
  }
  errors <- do.call(rbind, res)
  errors$nb_SNPs <- factor(errors$nb_SNPs, levels = nb_snps_to_sample)
  return(errors)
}

#' Evaluate the effect of SNP subsampling on genotype frequency estimation
#'
#' This function evaluates the impact of subsampling different numbers of SNPs on the estimation
#' of genotype frequencies. For each SNP subset size (provided in `nb_snps_vector`), the function:
#' - Estimates genotype frequencies
#' - Merges all estimates for comparison
#' - Produces summary statistics, boxplots, correlation plots, and bias analyses
#'
#' @param genotyping_matrix Genotyping matrix with values 1 (hom. REF), 0.5 (het), and 0 (hom. ALT) (SNPs x components)
#' @param allele_freqs Allele frequency matrix (SNPs x mixtures)
#' @param snp_depths Read depths matrix (SNPs x mixtures)
#' @param expected_freqs_melt Melted dataframe of expected genotype frequencies
#' @param step_size Number of SNPs to increment by when subsampling
#' @param nb_reps Integer. Number of random subsampling replicates per SNP subset size (default = 5).
#' @param out_dir Path to the output directory to save plots and result files
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{stats}{A list of dataframes from \code{compute_stats_per_exp_freq()}: mean and standard deviation per expected frequency}
#'   \item{bias}{A list of dataframes from \code{compute_bias()}: bias per component and per expected frequency}
#' }
#' @export
evaluate_subsampling_effect<-function(genotyping_matrix, allele_freqs, snp_depths, expected_freqs_melt, step_size, nb_reps=5, out_dir){
  writeLines("\n\nCompare the effect of random SNP subsampling on genotype frequency estimation:\n")

  freqs_to_compare_subsampling<-estimate_geno_freqs_snps_subsampling(
    genotyping_matrix = genotyping_matrix,
    allele_freqs = allele_freqs,
    snp_depths = snp_depths,
    expected_freqs_melt = expected_freqs_melt,
    step_size = step_size
  )

  compare_different_est_freqs(
    freqs_to_compare = freqs_to_compare_subsampling,
    variable_name = "Number of SNPs",
    out_dir = out_dir,
    suffix = "_subsampling_effect"
  )

  errors<-estimate_freq_errors_subsampling(
    genotyping_matrix = genotyping_matrix,
    allele_freqs = allele_freqs,
    snp_depths = snp_depths,
    expected_freqs_melt = expected_freqs_melt,
    step_size = step_size,
    nb_reps = nb_reps
  )
  
  plot_boxplots(
    df = errors,
    x_col = "nb_SNPs",
    y_col = "error",
    x_lab = "Number of SNPs",
    y_lab = "Error (estimated - expected)",
    title = "Estimation error under SNP subsampling",
    subtitle = glue("n = {nb_reps} replicates per SNP count"),
    show_x_labels = TRUE,
    out_file = glue("{out_dir}/subsampling_freq_errors_boxplot.png")
  )
}


#' Plot histogram of Minor Allele Frequencies (MAF) from a VCF object
#'
#' This function extracts the minor allele frequencies (MAF) from a `vcfR` object
#' and plots a histogram of their distribution.
#'
#' @param vcf A `vcfR` object, as read by [vcfR::read.vcfR()].
#' @param x_lim_values A numeric vector of length 2 indicating the minimum and maximum x-axis limits (e.g., `c(0, 0.5)`).
#'
#' @return A histogram plot of the MAF distribution. No return value; the function is used for its side effect (plot).
#' @export
plot_MAF_hist<-function(vcf, x_lim_values, out_dir = NULL){
  min_x_lim<-x_lim_values[1]
  max_x_lim<-x_lim_values[2]
  maf<-maf(vcf)
  maf <- as.data.frame(cbind(row.names(maf), maf[,"Frequency"]), stringsAsFactors = F)
  colnames(maf)<-c("CHR_POS", "MAF")
  maf$MAF<-as.numeric(maf$MAF)


  if (!is.null(out_dir)){ png(glue("{out_dir}/MAF_hist.png"), width = 800, height = 600) }
  par(
    cex.axis = 1.2,
    cex.lab = 1.4,
    cex.main = 1.5
  )
  hist(maf$MAF,
       ylab = "Number of SNPs",
       xlab = "Minor Allele Frequency",
       main = "",
       breaks = seq(0, 1, by = 0.1),
       xlim = c(min_x_lim, max_x_lim),
       xaxt = "n")
  axis(1, at = seq(min_x_lim, max_x_lim, by = 0.1))
  if (!is.null(out_dir)){
    dev.off()
    write.table(maf, glue("{out_dir}/MAFs.tsv"), row.names = FALSE, quote = F, sep ="\t")
    writeLines(c(
      "",
      "Output files saved in:",
      glue("{out_dir}/MAF_hist.png"),
      glue("{out_dir}/MAFs.tsv"),
      ""
    ))
  } else {
      return(maf)
    }

}


#' Plot marker set intersections as an UpSet and save to PNG.
#'
#' @param vcf A numeric genotyping matrix (rows = SNPs, columns = samples) where 1 = homozygous for the major allele, 0 = homozygous for the minor allele. Typically produced by [vcf_to_majmin_numeric_matrix()].
#' @param out_file Output PNG file path.
#' @return An UpSet plot of the marker set intersections. No return value; the function is used for its side effect (plot).
#' @export
plot_marker_set_intersections <- function(numeric_matrix, out_file = NULL) {

  df <- as.data.frame(numeric_matrix, stringsAsFactors = FALSE)

  non_binary_rows <- apply(df, 1, function(row) {
    vals <- suppressWarnings(as.numeric(row))
    any(!(vals %in% c(0, 1)), na.rm = TRUE)
  })

  nb_removed <- sum(non_binary_rows)
  if (nb_removed > 0) {
    warning(sprintf("[%s] Removed %d marker(s) presenting genotyping values other than 0/1.",
                    as.character(sys.call()[[1]]), nb_removed),
            call. = FALSE)
    df <- df[!non_binary_rows, , drop = FALSE]
  }

  df[] <- lapply(df, function(x) ifelse(is.na(x), 0L, as.integer(x)))

  if (!is.null(out_file)) { png(out_file, width = 1600, height = 1000, res = 150) }

  p<-upset(
    df,
    nsets = ncol(df),
    nintersects = 40,
    order.by = "freq",
    keep.order = TRUE,
    text.scale = 1.8
  )
  print(p)

  if (!is.null(out_file)) {
    dev.off()
    writeLines(c("", "Output file saved in:", out_file, ""))
  }
}


#' Merge multiple depth files by CHROM and POS
#'
#' @param depth_list Path to a text file listing depth file paths (one per line)
#' Listed files must be read depths matrices (SNPs x samples) with CHROM and POS as first and second columns.
#'
#' @return A merged data.frame containing all depth files joined by CHROM and POS.
#' @export
merge_depth_files <- function(depth_list) {
  files <- readLines(depth_list)
  files <- files[files != ""]

  merged <- read.table(files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  key_cols <- colnames(merged)[1:2]

  for (f in files[-1]) {
    df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    merged <- merge(merged, df, by = key_cols, all = TRUE)
  }
  merged
}

#' Clean depth file 2 first columns
#'
#' @param depths Read depth data frame (SNPs x samples) with CHROM and POS as first two columns.
#'
#' @return A cleaned read depth data frame with chrom_pos as row names and no CHROM/POS columns.
#' @export
clean_depth_file <- function(depths){
  rownames(depths) <- paste(depths$CHROM, depths$POS, sep = "_")
  depths <- depths[, -c(1:2)]
  return(depths)
}


#' Plot boxplots of sequencing depth per marker
#'
#' @param depths Read depth data frame (SNPs x samples) with CHROM and POS as first two columns.
#' @param hline Optional numeric value for a horizontal reference line (e.g. a depth threshold).
#' @param out_file Optional path to save the plot as a PNG. If NULL, the plot is printed.
#'
#' @return A ggplot object (printed or saved).
#' @export
plot_depth_per_marker <- function(depths, hline = NULL, out_file = NULL) {
  
  depths <- clean_depth_file(depths)
  depths <- depths[names(sort(rowMeans(depths))), , drop = FALSE]
  
  df_long <- reshape2::melt(as.matrix(depths))
  colnames(df_long) <- c("Marker", "Sample", "Depth")
  
  plot_boxplots(
    df = df_long,
    x_col = "Marker",
    y_col = "Depth",
    x_lab = "Markers",
    y_lab = "Depth",
    title = "Distribution of depth per marker",
    hline = hline,
    show_x_labels = FALSE,
    out_file = out_file
  )
}


#' Compute and mean depth per marker
#'
#' @param depths Read depth data frame (SNPs x samples) with CHROM and POS as first two columns.
#' @param out_file Optional Path to save the output as a TSV.
#'
#' @return A data frame of mean depths per marker (printed or saved).
#' @export
compute_mean_depth_per_marker <- function(depths, out_file = NULL){
  depths <- clean_depth_file(depths)
  mean_depths<-as.data.frame(cbind(row.names(depths),round(rowMeans(depths), 2)))
  colnames(mean_depths)<-c("CHROM_POS", "mean_depth")
  if (! is.null(out_file)) {
    write.table(mean_depths, out_file, row.names = FALSE, quote = F, sep ="\t")
    writeLines(c("", "Output file saved in:", out_file, ""))
  } else {
    return(mean_depths)
  }
}