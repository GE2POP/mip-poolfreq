# genotypeFreqEstimation

*R toolkit for estimating genotype frequencies in mixtures of known components*  

The contribution of each component to every mixture is inferred from a VCF file describing the component genotypes at a set of SNPs, combined with the allele frequencies observed at the same SNPs in the mixtures, and their corresponding sequencing depths to weight the estimation. Estimated frequencies can optionally be compared with theoretical expectations.

---

## Requirements

- **R ≥ 4.1**
- The scripts automatically install missing dependencies:

`optparse, devtools, this.path, vcfR, reshape2,
ggplot2, UpSetR, glue, scales, multcompView`

---

## Scripts and usage

### 1. `estimate_genotype_frequencies.R`

Estimate genotype frequencies per SNP based on allele frequencies and depth information.

#### Arguments
```
-v, --vcf              Components VCF file
-a, --allele_freqs     Mixtures allele frequencies TSV file (including columns CHROM, POS)
-d, --depths           Mixtures depth per SNP TSV file (including columns CHROM, POS)
-l, --libs             Library name correspondence TSV file *(optional)*
-o, --output_file      Genotype frequencies output TSV file
```

#### Example
```
Rscript scripts/estimate_genotype_frequencies.R \
  -v data/genotypes.vcf \
  -a data/allele_freqs.tsv \
  -d data/depths.tsv \
  -l data/libnames.tsv \
  -o results/genotype_frequencies.tsv
```

#### Output: `genotype_frequencies.tsv`
Contains the **estimated genotype frequencies** of each component within each mixture.
Each row corresponds to one **component–mixture** pair, with the estimated frequency of that component in the mixture.
```text
Component	Mixture	EstFreq
EL4X_35	Tm1310	0.31356190061712
GQ4X_83	Tm1310	0.239210978119877
EL4X_199	Tm1310	0.23199517936899
EL4X_482	Tm1310	0.23539587920789
EL4X_35	Tm1313	0.0637226242816768
GQ4X_83	Tm1313	0.313386942157224
EL4X_199	Tm1313	0.30785570431057
EL4X_482	Tm1313	0.348394182207116
...
```


### 2. `compare_with_expected_frequencies.R`

Compare estimated genotype frequencies to expected values and compute bias, variance, and correlation metrics.

#### Arguments
```
-v, --vcf                      Components VCF file
--allele_freqs_mix             Mixtures allele frequencies TSV file
--depths_mix                   Mixtures depth per SNP TSV file (including columns CHROM, POS)
--exp_freqs_mix                Mixtures expected genotype frequencies TSV file
--allele_freqs_comp            Components allele frequencies TSV file *(optional)*
--depths_comp                  Components depth per SNP TSV file (including columns CHROM, POS) *(optional)*
--exp_freqs_comp               Components expected genotype frequencies TSV file *(optional)*
-l, --libs                     Library name correspondence TSV file *(optional)*
-o, --out_dir                  Output directory
```

#### Example
```
Rscript scripts/compare_with_expected_frequencies.R \
  -v data/genotypes.vcf \
  --allele_freqs_mix data/ALL_FREQ/27MIXTURES/ref_allelic_freqs.tsv \
  --depths_mix data/ALL_FREQ/27MIXTURES/total_depths.tsv \
  --exp_freqs_mix data/infos/expectedGenoFreqs_mixtures.tsv \
  --allele_freqs_comp data/ALL_FREQ/4COMPONENTS/ref_allelic_freqs.tsv \
  --depths_comp data/ALL_FREQ/4COMPONENTS/total_depths.tsv \
  --exp_freqs_comp data/infos/expectedGenoFreqs_components.tsv \
  -l data/infos/corresp_comp_genotypes_libnames.tsv \
  -o results/compare_with_expected_frequencies
```

#### Outputs
•	general/ — bias and correlation summary tables and plots (correlations.tsv, bias_summary.tsv, freq_comparison_scatter.png, genotype_bias_boxplot.png).  
•	weight_vector_effect/ — analysis of weight effects (weight_vs_bias.png, stats.tsv).  
•	snp_subsampling_effect/ — analysis of SNP subsampling impact (subsampling_correlation.png, subsampling_variance.tsv).  


### 3. `plot_MAF_hist.R`
Compute Minor Allele Frequency (MAF) values from a VCF file and plot their distribution.

#### Arguments
```
-v, --vcf      Components VCF file
-l, --libs     Library name correspondence TSV file *(optional)*
-o, --out_dir  Output directory
```

#### Example
```
Rscript scripts/plot_MAF_hist.R \
  -v data/genotypes.vcf \
  -l data/infos/corresp_comp_genotypes_libnames.tsv \
  -o results/plot_MAF_hist
```

#### Outputs
•	MAF_hist.png — histogram of MAF values.  
•	MAFs.tsv — table of MAF values.  
•	marker_set_upsetplot.png — UpSet plot showing the number of SNPs are shared or unique among the components.


### 4. `plot_depth_per_marker.R`
Visualize boxplots of read depth distributions per SNP.

#### Arguments
```
-d, --depth_files_list   Text file listing all depth TSV files
-l, --hline              Read depth reference value to show as a horizontal line on the plot *(optional)*
-o, --out_dir            Output directory
```

#### Example
```
Rscript scripts/plot_depth_per_marker.R \
  -d data/ALL_FREQ/depth_files.list \
  -l 50 \
  -o results/plot_depth_per_marker
```

#### Outputs
•	mean_depth_per_marker.tsv — mean read depth per SNP across samples.
•	depth_per_marker_boxplots.png — boxplot of depth distributions, with optional horizontal reference line.


---

## Input files

### • VCF file (`--vcf`)

Standard multi-sample **VCF** (biallelic SNPs only) providing the genotypes of the component samples.  
Sample names correspond to libraries listed in the `--libs` file (if provided)


### • Allele frequency file (`--allele_freqs`)

TSV file giving the **reference allele frequencies at each SNP in the mixtures**.  
Used to estimate genotype frequencies from observed allele proportions.

**Expected columns**

| Column | Description |
|---------|-------------|
| `CHROM` | Chromosome or contig name |
| `POS` | SNP position |
| One column per mixture | Allele frequency of the REF allele in each mixture |


### • Read depth file (`--depths`)

TSV file giving the **sequencing depth per SNP in the mixtures**.  
Used to weight allele frequencies during genotype frequency estimation.

**Expected columns**

| Column | Description |
|---------|-------------|
| `CHROM` | Chromosome or contig name |
| `POS` | SNP position |
| One column per mixture | Read depth at this SNP in each mixture |

### • Expected genotype frequency file (`--exp_freqs`)

TSV file giving the **expected genotype proportions of each component within each mixture**.  
Each row corresponds to one mixture, and the values represent the expected fraction (0–1) of each component genotype contributing to it.  
Used as the theoretical genotype frequencies that estimated values will be compared to.

**Expected columns**

| Column | Description |
|---------|-------------|
| `Mixture` | Identifier of the mixture |
| One column per component | Expected proportion of each component in the mixture |

### • Library mapping file (`--libs`)

Optional TSV file used to rename VCF library IDs with genotype names in the outputs and figures.  
It does not affect computations, only the display of sample names.

**Expected columns**

| Column | Description |
|---------|-------------|
| `Sample_ID` | Library name as found in the VCF |
| `Library_Name` | Corresponding genotype name to display in outputs |








