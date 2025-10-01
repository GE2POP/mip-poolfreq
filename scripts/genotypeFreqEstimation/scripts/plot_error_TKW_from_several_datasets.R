## NB
# - tous les datasets ont le même fichier de TKW !!
# - il faudra sélectionner les bons témoins dans expectedGenoFreqs, ref_allelic et total_depths de R2022 et R2023

files_dir="~/MIPs/PUBLI/data_for_bias_TKW_corr"
genotyping_vcf_path=paste0(files_dir, "/GENOTYPING_MATRIX/04__Genotype_Locus1_Sample_Locus2_Filtered.vcf")
lib_names_corresp_path=paste0(files_dir, "/info_files/corresp_comp_genotypes_libnames.tsv")
allele_freqs_paths=c(
  test3=paste0(files_dir, "/test3/ref_allelic_freqs_test3_246SNPs.tsv"),
  R2022=paste0(files_dir, "/2022/ref_allelic_freqs_R2022_243SNPs.tsv"),
  R2023=paste0(files_dir, "/2023/ref_allelic_freqs_R2023_244SNPs.tsv")
)
snp_depths_paths=c(
  test3=paste0(files_dir, "/test3/total_depths_test3_246SNPs.tsv"),
  R2022=paste0(files_dir, "/2022/total_depths_R2022_243SNPs.tsv"),
  R2023=paste0(files_dir, "/2023/total_depths_R2023_244SNPs.tsv")
)
expected_freqs_paths=c(
  test3=paste0(files_dir, "/test3/expectedGenoFreqs_3témoins.tsv"),
  R2022=paste0(files_dir, "/2022/expectedGenoFreqs_5temoins.tsv"),
  R2023=paste0(files_dir, "/2023/expectedGenoFreqs_6temoins.tsv")
)
# tkw_paths=c(
#   test3=paste0(files_dir, "/test3/components_PMG.txt"),
#   R2022=paste0(files_dir, "/2022/components_PMG.txt"),
#   R2023=paste0(files_dir, "/2023/components_PMG.txt")
# )
tkw_paths=c(
  test3=paste0(files_dir, "/test3/components_TKW.tsv"),
  R2022=paste0(files_dir, "/2022/components_TKW.tsv"),
  R2023=paste0(files_dir, "/2023/components_TKW.tsv")
)

## Import dependencies
user_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(user_lib)

pkgs<-c(
  "optparse",
  "devtools",
  "this.path",
  "glue",
  "vcfR",
  "reshape2",
  "scales",
  "ggplot2",
  "multcompView",
  "tidyverse",
  "dplyr",
  "lme4"
)

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing missing package: ", pkg)
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

script_dir <- dirname(this.path::this.path())
project_root <- normalizePath(file.path(script_dir, ".."))
setwd(project_root)

load_all(".")


required_files=c(genotyping_vcf_path=genotyping_vcf_path, lib_names_corresp_path=lib_names_corresp_path, 
                 allele_freqs_paths, snp_depths_paths, expected_freqs_paths, tkw_paths)

check_input_files(
  required_files = required_files
)


## Import input files
genotyping_matrix<-vcf_to_numeric_matrix(
  vcf_path = genotyping_vcf_path,
  lib_names_corresp_path = lib_names_corresp_path
)

allele_freqs_dfs<-c()
snp_depths_dfs<-c()
expected_freqs_melt_dfs<-c()
tkw_dfs<-c()
for (dataset in c("test3", "R2022", "R2023")){
  allele_freqs<-read.table(allele_freqs_paths[dataset], header=T)
  snp_depths<-read.table(snp_depths_paths[dataset], header=T)
  expected_freqs<-read.table(expected_freqs_paths[dataset], header=T)
  tkw<-read.table(tkw_paths[dataset], header=T)
  
  allele_freqs<-set_rownames_from_chrom_pos(allele_freqs)
  snp_depths<-set_rownames_from_chrom_pos(snp_depths)
  expected_freqs_melt<-melt_genotype_freqs(expected_freqs, "ExpFreq")
  
  allele_freqs_dfs[[dataset]]<-allele_freqs
  snp_depths_dfs[[dataset]]<-snp_depths
  expected_freqs_melt_dfs[[dataset]]<-expected_freqs_melt
  tkw_dfs[[dataset]]<-tkw
}

# Estimate genotype frequencies and compute biases
errors_dfs<-c()
for (dataset in c("test3", "R2022", "R2023")){
  genotype_frequencies_mixtures<-estimate_genotype_freqs(
    genotyping_matrix = genotyping_matrix, 
    allele_freqs = allele_freqs_dfs[[dataset]], 
    snp_depths = snp_depths_dfs[[dataset]],
    expected_freqs_melt = expected_freqs_melt_dfs[[dataset]]
  )
  
  errors_output<-compute_errors(
    freqs_df = genotype_frequencies_mixtures
  )
  errors_dfs[[dataset]]<-errors_output
}

# Plot bias vs TKW
bias_TKW <- map2_dfr(
  errors_dfs,
  tkw_dfs,
  ~ inner_join(.x, .y, by = c("Mixture", "Component")),
  .id = "dataset"
)



ggplot(bias_TKW, aes(x = TKW, y = error_EstFreq, color = Component, shape = dataset)) +
  geom_point(size = 3) +
  labs(
    x = "TKW",
    y = "Error",
    title = "Error vs TKW by Component and Dataset",
    color = "Component",
    shape = "Dataset"
  ) +
  theme_minimal()

summary(lm(error_EstFreq~TKW+Component+Origin, data = bias_TKW))
summary(lm(error_EstFreq~TKW+Component, data = bias_TKW))
summary(lm(error_EstFreq~TKW, data = bias_TKW))
summary(lm(error_EstFreq~Component, data = bias_TKW))
boxplot(error_EstFreq~Origin, data = bias_TKW) #trop de bacs (modèle surparamétré et sans doute overfitté)

anova(lmer(error_EstFreq ~ TKW + Component + (1|Origin), data = bias_TKW)) #Origin passée en effet aléatoire

