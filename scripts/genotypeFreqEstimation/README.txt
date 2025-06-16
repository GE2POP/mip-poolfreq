cd /mnt/c/Users/girodolle/Documents/MIPs/PUBLI
genotypeFreqEstimation_path="/mnt/c/Users/girodolle/Documents/MIPs/PUBLI/03_scripts/genotypeFreqEstimation"
scripts_dir=${genotypeFreqEstimation_path}/scripts

# run estimate_genotype_frequencies.R
vcf="${genotypeFreqEstimation_path}/tests/data/GENOTYPING_MATRIX/04__Genotype_Locus1_Sample_Locus2_Filtered.vcf"
allele_freqs="${genotypeFreqEstimation_path}/tests/data/ALL_FREQ/27MIXTURES/ref_allelic_freqs_27mixtures_246SNPs.tsv"
depths="${genotypeFreqEstimation_path}/tests/data/ALL_FREQ/27MIXTURES/total_depths_27mixtures_246SNPs.tsv"
libs="${genotypeFreqEstimation_path}/tests/data/infos_files/corresp_comp_genotypes_libnames.tsv"
output_file="${genotypeFreqEstimation_path}/tests/results/estimate_genotype_frequencies/genotype_frequencies.tsv"

Rscript ${scripts_dir}/estimate_genotype_frequencies.R \
  -v $vcf \
  -a $allele_freqs \
  -d $depths \
  -l $libs \
  -o $output_file

# run compare_with_expected_frequencies.R
allele_freqs_mix="${genotypeFreqEstimation_path}/tests/data/ALL_FREQ/27MIXTURES/ref_allelic_freqs_27mixtures_246SNPs.tsv"
depths_mix="${genotypeFreqEstimation_path}/tests/data/ALL_FREQ/27MIXTURES/total_depths_27mixtures_246SNPs.tsv"
exp_freqs_mix="${genotypeFreqEstimation_path}/tests/data/infos_files/expectedGenoFreqs_mixtures.tsv"
allele_freqs_comp="${genotypeFreqEstimation_path}/tests/data/ALL_FREQ/4COMPONENTS/ref_allelic_freqs_4components_246SNPs.tsv"
depths_comp="${genotypeFreqEstimation_path}/tests/data/ALL_FREQ/4COMPONENTS/total_depths_4components_246SNPs.tsv"
exp_freqs_comp="${genotypeFreqEstimation_path}/tests/data/infos_files/expectedGenoFreqs_components.tsv"
lib_names_corresp_path="${genotypeFreqEstimation_path}/tests/data/infos_files/corresp_comp_genotypes_libnames.tsv"
out_dir="${genotypeFreqEstimation_path}/tests/results/compare_with_expected_frequencies"

Rscript ${scripts_dir}/compare_with_expected_frequencies.R \
  -v $vcf \
  --allele_freqs_mix $allele_freqs_mix \
  --depths_mix $depths_mix \
  --exp_freqs_mix $exp_freqs_mix \
  --allele_freqs_comp $allele_freqs_comp \
  --depths_comp $depths_comp \
  --exp_freqs_comp $exp_freqs_comp \
  -l $libs \
  -o $out_dir
  
# run plot_MAF_hist.R
vcf="${genotypeFreqEstimation_path}/tests/data/GENOTYPING_MATRIX/04__Genotype_Locus1_Sample_Locus2_Filtered.vcf"
output_file="${genotypeFreqEstimation_path}/tests/results/plot_MAF_hist/MAF_hist.png"

Rscript ${scripts_dir}/plot_MAF_hist.R \
  -v $vcf \
  -o $output_file
