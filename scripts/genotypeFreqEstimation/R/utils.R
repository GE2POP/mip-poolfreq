#' Check that all provided paths are non-NULL
#'
#' @param paths Named list of paths to check
#'
#' @return Stops if any value is NULL
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
check_file_exists <- function(name, path) {
  if (!file.exists(path)) {
    stop(sprintf("File or folder does not exist:\n- %s: %s", name, path), call. = FALSE)
  }
}

#' Check required and optional input file paths
#'
#' @param required_files Named list of required file/folder paths
#' @param optional_files Named list of optional file/folder paths (default: NULL)
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

#' Check that required columns are present in a data frame
#'
#' Verifies that all expected column names are present in the input data frame.
#'
#' @param df A data frame to check.
#' @param required_cols A character vector of required column names.
#'
#' @return Nothing.
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
melt_genotype_freqs<-function(genotype_freqs, freq_col_name){
  genotype_freqs_melt<-melt(genotype_freqs, id.vars = "Mixture")
  colnames(genotype_freqs_melt)[2:3]<-c("Component", freq_col_name)
  return(genotype_freqs_melt)
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
estimate_genotype_freqs <- function(genotyping_matrix, allele_freqs, snp_depths = NULL, expected_freqs_melt = NULL, est_freq_col_name = "EstFreq", out_file = NULL) {
  nb_mixtures <- ncol(allele_freqs)
  nb_components <- ncol(genotyping_matrix)
  genotype_freqs <- data.frame(matrix(nrow = nb_components, ncol = nb_mixtures))
  colnames(genotype_freqs) <- colnames(allele_freqs)
  rownames(genotype_freqs) <- colnames(genotyping_matrix)
  
  allele_freqs <- allele_freqs[rownames(genotyping_matrix),]
  snp_depths <- if (!is.null(snp_depths)) snp_depths[rownames(genotyping_matrix),] else NULL
  
  for (mixture_name in colnames(allele_freqs)) {
    weights <- if (!is.null(snp_depths)) snp_depths[, mixture_name] else NULL
    est_freqs <- coefficients(lm(allele_freqs[, mixture_name] ~ genotyping_matrix - 1, weights = weights))
    genotype_freqs[, mixture_name] <- est_freqs
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
#'
#' @return A ggplot object
plot_condition_effect_boxplots <- function(df_melt, x_col, y_col, fill_col, ref_col = NULL, y_lab = "", x_lab = "", fill_lab = "Group") {
  # Set fill as factor to ensure ordering
  df_melt[[fill_col]] <- factor(df_melt[[fill_col]], levels = sort(unique(df_melt[[fill_col]])))
  
  # Define palette
  fill_levels <- levels(df_melt[[fill_col]])
  color_palette <- scales::hue_pal()(length(fill_levels))
  names(color_palette) <- fill_levels
  
  p <- ggplot(df_melt, aes_string(x = x_col, y = y_col, fill = fill_col)) +
    geom_boxplot() +
    scale_fill_manual(values = color_palette) +
    xlab(x_lab) +
    ylab(y_lab) +
    labs(fill = fill_lab) +
    theme_minimal() +
    scale_x_discrete(expand = expansion(add = c(1, 0)))
  
  # Optional reference lines and labels
  if (!is.null(ref_col)) {
    p <- p +
      geom_hline(aes_string(yintercept = paste0("as.numeric(as.character(", ref_col, "))"),
                            color = fill_col),
                 linetype = "dashed") +
      geom_text(aes_string(x = 0,
                           y = paste0("as.numeric(as.character(", ref_col, "))"),
                           label = ref_col,
                           color = fill_col),
                hjust = 0, vjust = -0.5, size = 3) +
      scale_color_manual(values = color_palette, guide = "none")
  }
  
  return(p)
}

#' Plot boxplots of (1 or more sets of) estimated frequencies by expected frequencies
#'
#' @param freqs_df Dataframe with at least one estimated frequency column and a "ExpFreq" column (providing the expected frequencies)
#' @param variable_name Label for x-axis (below all conditions)
#' @param out_file Path to the output file to save the plot (default = NULL)
#' @return A ggplot object
compare_boxplots_est_freq <- function(freqs_df, variable_name, out_file = NULL) {
  df_melt <- melt(freqs_df[, -c(1:2)], id.vars = "ExpFreq")
  colnames(df_melt)[2:3] <- c("Condition", "EstFreq")
  
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
    ggsave(out_file, plot = p, width = 8, height = 6, bg = "white")
  }
  
  return(p)
}



#' Compute mean and SD per expected frequency
#'
#' @param freqs_df Dataframe with an "ExpFreq" colomn (providing the expected frequencies) and 1 or more estimated frequencies column(s)
#' @param out_dir Path the the output directory to write output files (default = NULL)
#' @param suffix Suffix to add in the output file names (default = "")
#' @return A list with two dataframes: expected frequencies mean and standard deviation
compute_stats_per_exp_freq <- function(freqs_df, out_dir = NULL, suffix = "") {
  output <- list()
  df_melt <- melt(freqs_df[, -c(1:2)], id.vars = "ExpFreq")
  mean_df <- aggregate(value ~ ExpFreq + variable, data = df_melt, mean)
  sd_df <- aggregate(value ~ ExpFreq + variable, data = df_melt, sd)
  output$mean <- dcast(mean_df, ExpFreq ~ variable)
  output$sd <- dcast(sd_df, ExpFreq ~ variable)
  
  if (! is.null(out_dir)){
    write.table(output$mean, glue("{out_dir}/est_geno_freqs_mean{suffix}.tsv"), row.names = FALSE, quote = F, sep ="\t")
    write.table(output$sd, glue("{out_dir}/est_geno_freqs_sd{suffix}.tsv"), row.names = FALSE, quote = F, sep ="\t")
  }
  
  return(output)
}

#' Plot estimated vs expected frequencies with regression
#'
#' @param freqs_df Dataframe with an "ExpFreq" colomn (providing the expected frequencies) and 1 or more estimated frequencies column(s)
#' @param extra_freqs_df Dataframe with additional points to plot (they won't be used to compute the regression) (default = NULL)
#' @param out_file Path to the output file to save the plot (default = NULL)
plot_correlation_with_exp_freq <- function(freqs_df, extra_freqs_df = NULL, out_file = NULL) {
  
  if (! is.null(out_file)) { png(out_file, width = 800, height = 600) }
  
  for (condition in colnames(freqs_df)[!colnames(freqs_df) %in% c("Component", "Mixture", "ExpFreq")]) {
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
      r2 <- cor(freqs_df[, condition], freqs_df$ExpFreq)^2
    } else {
      r2 <- NA
    }
    
    text(0.2, 0.8, paste0("Slope: ", round(coeffs[2], 3)), col = "darkred")
    text(0.2, 0.75, paste0("Intercept: ", round(coeffs[1], 3)), col = "darkred")
    text(0.2, 0.7, paste0("R^2: ", round(r2, 3)), col = "darkred")
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
compute_bias_from_errors <- function(freqs_df, error_col, condition) {
  stopifnot(startsWith(error_col, "error_"))
  
  formula_comp <- as.formula(glue("{error_col} ~ Component"))
  formula_exp  <- as.formula(glue("{error_col} ~ ExpFreq"))
  
  comp_mean <- aggregate(formula_comp, freqs_df, mean)
  expf_mean <- aggregate(formula_exp, freqs_df, mean)
  
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
plot_bias_boxplots <- function(
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
    fill_lab = "Component"
  )
  
  p_expfreq <- plot_condition_effect_boxplots(
    df_melt = errors_long_df,
    x_col = "Condition",
    y_col = "Error",
    fill_col = "ExpFreq",
    y_lab = "Error (estimated - expected)",
    x_lab = variable_name,
    fill_lab = "Expected frequency"
  )
  
  if (!is.null(out_dir)) {
    ggsave(glue("{out_dir}/bias_by_component{suffix}.png"), p_component, width = 8, height = 6, dpi = 300, bg = "white")
    ggsave(glue("{out_dir}/bias_by_exp_freq{suffix}.png"), p_expfreq, width = 8, height = 6, dpi = 300, bg = "white")
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
    ) %>%
    dplyr::mutate(Condition = factor(Condition, levels = est_freq_cols))
  
  plot_bias_boxplots(errors_long_df, variable_name, out_dir, suffix)
  
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
    out_file = glue("{out_dir}/correlation_plot{suffix}.png")
  )
  
  bias <- compute_bias(
    freqs_df = freqs_to_compare,
    variable_name = variable_name,
    out_dir = out_dir,
    suffix = suffix
  )
  
  return(list(stats = stats, bias = bias))
  
}


#' Compare the effect of weight vector on genotype frequency estimation
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
compare_weight_vector_effect<-function(genotyping_matrix, allele_freqs, snp_depths, expected_freqs_melt, out_dir){
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
#' @param expected_freqs_melt Melted dataframe of expected genotype frequencies (default = NULL)
#' @return Dataframe of estimated frequencies per condition
estimate_geno_freqs_snps_subsampling <- function(genotyping_matrix, allele_freqs, step_size, snp_depths = NULL, expected_freqs_melt = NULL) {
  total_nb_snps<-dim(genotyping_matrix)[1]
  nb_snps_to_sample<-unique(c(seq(step_size, total_nb_snps, by = step_size), total_nb_snps))
  
  subsampled_genotype_freqs <- data.frame(Mixture = colnames(allele_freqs))
  
  for (nb_snps in nb_snps_to_sample) {
    condition <- glue("{nb_snps}_SNPs")
    snps <- rownames(genotyping_matrix)[sample(nrow(genotyping_matrix), nb_snps)]
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
  
  return(subsampled_genotype_freqs)
}

#' Compare the effect of SNP subsampling on genotype frequency estimation
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
#' @param out_dir Path to the output directory to save plots and result files
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{stats}{A list of dataframes from \code{compute_stats_per_exp_freq()}: mean and standard deviation per expected frequency}
#'   \item{bias}{A list of dataframes from \code{compute_bias()}: bias per component and per expected frequency}
#' }
compare_subsampling_effect<-function(genotyping_matrix, allele_freqs, snp_depths, expected_freqs_melt, step_size, out_dir){
  writeLines("\n\nCompare the effect of random SNP subsampling on genotype frequency estimation:\n")
  
  freqs_to_compare_subsampling<-estimate_geno_freqs_snps_subsampling(
    genotyping_matrix = genotyping_matrix, 
    allele_freqs = allele_freqs, 
    snp_depths = snp_depths, 
    expected_freqs_melt = expected_freqs_melt, 
    step_size = step_size
  )
  
  stats<-compare_different_est_freqs(
    freqs_to_compare = freqs_to_compare_subsampling,
    variable_name = "Number of SNPs",
    out_dir = out_dir,
    suffix = "_subsampling_effect"
  )
  return(stats)
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
plot_MAF_hist<-function(vcf, x_lim_values, out_file = NULL){
  min_x_lim<-x_lim_values[1]
  max_x_lim<-x_lim_values[2]
  maf_values <- maf(vcf)
  
  if (!is.null(out_file)){ png(out_file, width = 800, height = 600) }
  hist(maf_values[,4],
       ylab = "Number of SNPs",
       xlab = "Minor Allele Frequency",
       main = "",
       breaks = seq(0, 1, by = 0.1),
       xlim = c(min_x_lim, max_x_lim),
       xaxt = "n")
  axis(1, at = seq(min_x_lim, max_x_lim, by = 0.1))
  if (!is.null(out_file)){ 
    dev.off()
    writeLines(c(
      "",
      "Output file saved in:",
      glue("{out_file}"),
      ""
    ))
    }
  
}