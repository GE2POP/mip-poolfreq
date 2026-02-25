#' Run the full frequency accuracy evaluation pipeline
#'
#' This function orchestrates the complete workflow to estimate genotype
#' frequencies, compare them to expected frequencies, compute biases,
#' and evaluate weight vector and SNP subsampling effects.
#'
#' @param genotyping_vcf_path Path to the genotyping VCF file
#' @param allele_freqs_mixtures_path Path to allele frequencies file for mixtures
#' @param snp_depths_mixtures_path Path to SNP depths file for mixtures
#' @param expected_freqs_mixtures_path Path to expected genotype frequencies for mixtures
#' @param allele_freqs_components_path Optional path to allele frequencies file for components
#' @param snp_depths_components_path Optional path to SNP depths file for components
#' @param expected_freqs_components_path Optional path to expected genotype frequencies for components
#' @param lib_names_corresp_path Optional path to library name correspondence file
#' @param min_depth Minimum read depth threshold
#' @param subsampling_step SNP subsampling step size
#' @param subsampling_reps Integer. Number of random subsampling replicates per SNP subset size (default = 5).
#' @param out_dir Output directory
#'
#' @return Invisibly returns a list with main output paths
#' @export
eval_pipeline <- function(
    genotyping_vcf_path,
    allele_freqs_mixtures_path,
    snp_depths_mixtures_path,
    expected_freqs_mixtures_path,
    allele_freqs_components_path = NULL,
    snp_depths_components_path = NULL,
    expected_freqs_components_path = NULL,
    lib_names_corresp_path = NULL,
    min_depth = 0,
    subsampling_step = 50,
    subsampling_reps = 5,
    out_dir
) {
    required_files <- list(
        vcf = genotyping_vcf_path,
        allele_freqs_mix = allele_freqs_mixtures_path,
        depths_mix = snp_depths_mixtures_path,
        exp_freqs_mix = expected_freqs_mixtures_path,
        out_dir = out_dir
    )

    optional_files <- list(
        libs = lib_names_corresp_path,
        allele_freqs_comp = allele_freqs_components_path,
        depths_comp = snp_depths_components_path,
        exp_freqs_comp = expected_freqs_components_path
    )

    check_missing_args(args = c(required_files, min_depth))

    check_input_files(
        required_files = required_files,
        optional_files = optional_files
    )

    inputs <- load_inputs(
        genotyping_vcf_path = genotyping_vcf_path,
        lib_names_corresp_path = lib_names_corresp_path,
        allele_freqs_mixtures_path = allele_freqs_mixtures_path,
        snp_depths_mixtures_path = snp_depths_mixtures_path,
        expected_freqs_mixtures_path = expected_freqs_mixtures_path,
        allele_freqs_components_path = allele_freqs_components_path,
        snp_depths_components_path = snp_depths_components_path,
        expected_freqs_components_path = expected_freqs_components_path,
        min_depth = min_depth
    )


    general_subdir<-glue("{out_dir}/expected_vs_estimated")
    weights_subdir<-glue("{out_dir}/weight_vector_effect")
    subsampling_subdir<-glue("{out_dir}/snp_subsampling_effect")

    for (subdir in c(general_subdir, weights_subdir, subsampling_subdir)){
        dir.create(subdir)
    }

    writeLines("\n\nEstimate genotype frequencies, compare to expected frequencies and compute biases:\n")
    genotype_frequencies_mixtures<-estimate_genotype_freqs(
        genotyping_matrix = inputs$genotyping_matrix,
        allele_freqs = inputs$allele_freqs_mixtures,
        snp_depths = inputs$snp_depths_mixtures,
        expected_freqs_melt = inputs$expected_freqs_mixtures_melt,
        out_file = glue("{general_subdir}/est_geno_freqs_mixtures.tsv")
    )

    compute_stats_per_exp_freq(
        freqs_df = genotype_frequencies_mixtures,
        out_dir = general_subdir
    )

    if (! is.null(inputs$allele_freqs_components)){
    genotype_frequencies_components<-estimate_genotype_freqs(
        genotyping_matrix = inputs$genotyping_matrix,
        allele_freqs = inputs$allele_freqs_components,
        snp_depths = inputs$snp_depths_components,
        expected_freqs_melt = inputs$expected_freqs_components_melt,
        out_file = glue("{general_subdir}/est_geno_freqs_components.tsv")
    )
    } else {
        genotype_frequencies_components<-NULL
    }

    plot_correlation_with_exp_freq(
        freqs_df = genotype_frequencies_mixtures,
        extra_freqs_df = genotype_frequencies_components,
        out_file = glue("{general_subdir}/correlation_plot.pdf")
    )

    compute_bias(
        freqs_df = genotype_frequencies_mixtures,
        out_dir = general_subdir
    )

    evaluate_weight_vector_effect(
        genotyping_matrix = inputs$genotyping_matrix,
        allele_freqs = inputs$allele_freqs_mixtures,
        snp_depths = inputs$snp_depths_mixtures,
        expected_freqs_melt = inputs$expected_freqs_mixtures_melt,
        out_dir = weights_subdir
    )

    evaluate_subsampling_effect(
        genotyping_matrix = inputs$genotyping_matrix,
        allele_freqs = inputs$allele_freqs_mixtures,
        snp_depths = inputs$snp_depths_mixtures,
        expected_freqs_melt = inputs$expected_freqs_mixtures_melt,
        step_size = subsampling_step,
        nb_reps = subsampling_reps,
        out_dir = subsampling_subdir
    )

    writeLines(c(
        "",
        "Output files saved in:",
        glue("- {general_subdir}"),
        glue("- {weights_subdir}"),
        glue("- {subsampling_subdir}"),
        ""
    ))

    invisible(list(
        expected_vs_estimated = general_subdir,
        weight_vector_effect = weights_subdir,
        snp_subsampling_effect = subsampling_subdir
    ))
}
