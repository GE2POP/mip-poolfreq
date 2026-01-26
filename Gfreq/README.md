# Gfreq

*An R package with a command-line interface for estimating genotype frequencies in mixtures of known components*

Gfreq infers the contribution of each component to a set of mixtures using a VCF file providing component genotypes at a set of SNPs, combined with allele frequencies observed at the same SNPs in the mixtures and their corresponding sequencing depths to weight the estimation. For artificial mixtures of known proportions, estimated frequencies can be compared with theoretical expectations.

---

## Table of Contents
- [Installation](#installation)
- [Command-line usage](#command-line-usage)
  - [1. estimate_genotype_frequencies](#1-estimate_genotype_frequencies)
  - [2. evaluate_frequency_accuracy](#2-evaluate_frequency_accuracy)
  - [3. plot_MAF_hist](#3-plot_maf_hist)
  - [4. plot_depth_per_marker](#4-plot_depth_per_marker)
- [Input files](#input-files)

---

## Installation

### Requirements

- **R ≥ 4.1**

### Install from source

```bash
git clone git@github.com:GE2POP/mip-poolfreq.git
R CMD INSTALL Gfreq
```

### Expose the command-line interface

After installation, add the Gfreq command-line tools to your PATH:
```bash
export PATH="$(Rscript -e 'cat(system.file("bin", package="Gfreq"))'):$PATH"
```

To make this permanent, add the line above to your ~/.bashrc.

### Check installation

```bash
Gfreq --help
```

---

## Command-line usage

*Example input and output files are provided in the example_data/ directory to illustrate the expected formats and typical results.*

### 1. `estimate_genotype_frequencies`

Estimate genotype frequencies per SNP based on allele frequencies and depth information.

#### Arguments
```
-v, --vcf              Components VCF file
-a, --allele_freqs     Mixture allele frequencies TSV file (including columns CHROM, POS)
-d, --depths           Mixture depths per SNP TSV file (including columns CHROM, POS)
-t, --min_depth        Minimum read depth threshold *(default = 0)*
-l, --libs             Library name correspondence TSV file *(optional)*
-o, --output_file      Genotype frequencies output TSV file
```

#### Example
```
cd Gfreq
Gfreq estimate_genotype_frequencies \
  -v example_data/input_files/comp_genotypes.vcf \
  -a example_data/input_files/mix_ref_all_freqs.tsv \
  -d example_data/input_files/mix_read_depths.tsv \
  -t 40 \
  -l example_data/input_files/comp_libnames_corresp.tsv \
  -o genotype_frequencies.tsv
```

#### SNP filtering by minimum depth
The `--min_depth` option removes all SNPs for which at least one mixture has a sequencing depth below the specified threshold.
This can help exclude low-confidence loci with unreliable allele frequency estimates due to poor sequencing coverage.

#### Output
• `genotype_frequencies.tsv`
**Estimated genotype frequencies** of each component within each mixture.
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


### 2. `evaluate_frequency_accuracy`

*Compare estimated genotype frequencies to expected values and compute bias, variance, and correlation metrics.*

This script is intended for test datasets composed of **artificial mixtures with known proportions**, where expected genotype frequencies are available for direct comparison with the estimated values.

Three complementary analyses will be performed:

1. **Direct comparison between expected and estimated frequencies**  
Estimated genotype frequencies in mixtures are compared to their expected values.  
In addition to estimating genotype frequencies in real mixtures, the same procedure can optionally be applied to the component libraries themselves, treated as if they were mixtures (thereafter referred to as “one-component mixtures”). To include them, the user must provide the components input files through the optional arguments --allele_freqs_comp, --depths_comp, and --exp_freqs_comp. This can serve as a control to verify that the estimation correctly returns a frequency of 1 for the corresponding component and 0 for the others.

3. **Effect of weighting in the regression model**  
The regression model used to estimate genotype frequencies can include a weight vector to assign more importance to certain observations (i. e. allele frequencies) than others.  
In our implementation, SNP read depths are used as the weight vector when estimating genotype frequencies in a given mixture. This allows to give greater influence to SNPs with higher sequencing depth, as these are expected to provide more reliable allele frequency estimates.  
This analysis compares results obtained with and without read depth weighting.  
*One-component mixtures are not included in this analysis.*

4. **Effect of reducing the number of SNPs used for estimation**  
To evaluate how decreasing the number of SNPs impacts estimation accuracy, SNPs are gradually subsampled, genotype frequencies are estimated across multiple replicates, and estimation error distributions are visualized using boxplots.  
*One-component mixtures are not included in this analysis.*

#### Arguments
```
-v, --vcf                      Components VCF file
--allele_freqs_mix             Mixture allele frequencies TSV file
--depths_mix                   Mixture depths per SNP TSV file (including columns CHROM, POS)
--exp_freqs_mix                Mixture expected genotype frequencies TSV file
-t, --min_depth                Minimum read depth threshold *(default = 0)*
--allele_freqs_comp            Component allele frequencies TSV file *(optional)*
--depths_comp                  Component depths per SNP TSV file (including columns CHROM, POS) *(optional)*
--exp_freqs_comp               Component expected genotype frequencies TSV file *(optional)*
-l, --libs                     Library name correspondence TSV file *(optional)*
-s, --subsampling_step         Number of SNPs added between subsampling iterations *(default = 50)*
-r, --subsampling_reps         Number of random subsampling replicates per SNP subset size *(default = 5)*
-o, --out_dir                  Output directory
```

#### Example
```
Gfreq evaluate_frequency_accuracy \
  -v example_data/input_files/comp_genotypes.vcf \
  --allele_freqs_mix example_data/input_files/mix_ref_all_freqs.tsv \
  --depths_mix example_data/input_files/mix_read_depths.tsv \
  --exp_freqs_mix example_data/input_files/mix_exp_geno_freqs.tsv \
  -t 40 \
  --allele_freqs_comp example_data/input_files/comp_ref_all_freqs.tsv \
  --depths_comp example_data/input_files/comp_read_depths.tsv \
  --exp_freqs_comp example_data/input_files/comp_exp_geno_freqs.tsv \
  -l example_data/input_files/comp_libnames_corresp.tsv \
  -s 5 \
  -r 3 \
  -o .
```

#### SNP filtering by minimum depth
As [previously described](#snp-filtering-by-minimum-depth), SNPs whose minimum sequencing depth across mixtures is below the provided threshold are excluded from the analysis.
Filtering is based exclusively on the mixture depth file (`--depths_mix`), but the same SNPs are also removed from the corresponding allele frequency and component tables to keep all inputs consistent before estimating genotype frequencies.

#### Outputs

- **`expected_vs_estimated` folder:**
*Comparison analysis of estimated vs expected genotype frequencies*
  - `est_geno_freqs_mixtures.tsv`: estimated genotype frequencies in mixtures
  - `est_geno_freqs_components.tsv`: estimated genotype frequencies in components
  - tables reporting mean and standard deviation of estimated frequencies (in mixtures) per expected value
  - tables of estimated biases (mean errors) per component and per expected frequency, with corresponding error boxplots
  - scatter plot of estimated vs expected genotype frequencies; each dot represents the frequency of one component in a mixture. If --allele_freqs_comp, --depths_comp and --exp_freqs_comp were provided, frequencies estimated in one-component mixtures (i.e. component libraries treated as mixtures for validation) are shown in grey and are excluded from the regression fitting.
<p align="center">
  <img width="532" height="388" alt="image" src="https://github.com/user-attachments/assets/5073bff6-db9e-4e2f-ade4-144c40e3982a" />
</p>

- **`weight_vector_effect` folder:**
*Analysis of the effect of using read depth as a weight in the estimation model*
  - tables reporting mean and standard deviation of estimated frequencies (in mixtures) per expected value, computed with or without using read depth as weight
  - tables of estimated biases (mean errors) per component and per expected frequency, with corresponding error boxplots, with or without using read depth as weight
  - scatter plot of estimated vs expected genotype frequencies, with or without using read depth as weight
  - `est_geno_freqs_boxplot_weight_effect.png`: boxplots of estimated values per expected frequency, with or without using read depth as weight
<p align="center">
  <img width="616" height="463" alt="image" src="https://github.com/user-attachments/assets/f7575565-95ef-49f9-ae8d-f430c699c440" />
</p>

- **`snp_subsampling_effect` folder:**
*Analysis of the effect of gradually reducing the number of SNPs through random subsampling on genotype frequency estimation*
  - `subsampling_freq_errors_boxplot.png`: boxplots of estimation errors (estimated minus expected genotype frequencies) across random SNP subsampling replicates, for increasing SNP subset sizes
<p align="center">

</p>


### 3. `plot_MAF_hist`
Compute Minor Allele Frequency (MAF) values from a VCF file and plot their distribution.

#### Arguments
```
-v, --vcf      Components VCF file
-l, --libs     Library name correspondence TSV file *(optional)*
-o, --out_dir  Output directory
```

#### Example
```
Gfreq plot_MAF_hist \
  -v example_data/input_files/comp_genotypes.vcf \
  -l example_data/input_files/comp_libnames_corresp.tsv \
  -o .
```

#### Outputs
- `MAF_hist.png`: histogram of MAF values.
- `MAFs.tsv`: table of MAF values.
- `marker_set_upsetplot.png`: UpSet plot showing the number of shared or unique SNPs among the components. **Only loci that are homozygous in every component are included.** For example, in the plot below, the first vertical bar indicates that GQ4X-83 carries a different homozygous allele than the other three components at 33 SNPs.
<p align="center">
  <img width="794" height="497" alt="image" src="https://github.com/user-attachments/assets/10195c83-bd6f-4837-96f8-75b0a5a2176f" />
</p>



### 4. `plot_depth_per_marker`
Plot read depth distribution (boxplots) per SNP.

#### Arguments
```
-d, --depth_files_list   Text file listing all depth TSV files
-l, --hline              Read depth reference value to show as a horizontal line on the plot *(optional)*
-o, --out_dir            Output directory
```

#### Example
```
Gfreq plot_depth_per_marker \
  -d example_data/input_files/depth_files.list \
  -l 50 \
  -o .
```

#### Outputs
- `mean_depth_per_marker.tsv`: mean read depth per SNP across samples.
- `depth_per_marker_boxplots.png`: boxplot of depth distributions, with optional horizontal reference line.
<p align="center">
  <img width="871" height="482" alt="image" src="https://github.com/user-attachments/assets/7c331c4b-c744-4a04-b347-ffb7961872f1" />
</p>

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
| One column per mixture | Allele frequency of the REF allele at this SNP in each mixture |


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
