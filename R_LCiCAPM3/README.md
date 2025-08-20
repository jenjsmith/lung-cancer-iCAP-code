---
title: "LC-iCAP Model M3 Generation"
---


# Summary: This repository contains analysis code and r markdown document for generating model M3 described in: A Multivariate Cell-Based Liquid Biopsy for Lung Nodule Risk Stratification: Analytical Validation and Early Clinical Evaluation

## Quick Start

If you just want to load and inspect the exact Model M3 from the manuscript:
```r

# Load original results list (unwrapped, legacy format)
res_3 <- readRDS("results/res_3.rds")

# Load wrapped results object (with results + meta)
res_3_wrapped <- readRDS("results/res_3_wrapped.rds")

# Extract the chosen model (Model 3)
model_3 <- res_3_wrapped$results[[ res_3_wrapped$meta$chosen_index ]]$final_model

# Inspect metadata for Model 3
# Note: in this analysis, seed number = index number (seed not stored in legacy output)
res_3_wrapped$meta
```

If you want to run blind set predictions exactly as in the paper:

```r
# Load input data
ex_pos <- readRDS("data/for_modeling/LCII.035.056.hk3-batch_corrected-clean-pos-3.0.20231016.rds")

# Run prediction function
make_predictions("results/res_3_wrapped.rds", ex_pos, "model_3")
```

>Note: Re-running model generation from scratch may produce slightly different models unless you use the chosen_seed from res_3_wrapped$meta in rf_nestedcv_stable().
>Note: Model M3’s chosen_seed may be NA in legacy runs; if you rerun from scratch, set seed in rf_nestedcv_stable() to the documented seed (or the run index) to reproduce the selected model.
---
## Table of Contents
- [Installation](#installation)
- [Project Structure](#project-structure)
- [Usage](#usage)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Training Setup](#Training Setup)
- [Blind Validation](#Blind Validation)
- [Reproducibility](#reproducibility)
- [Dependencies](#Dependencies)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

## Installation
Clone the repository and install required R packages:

```bash
git clone (https://github.com/jenjsmith/lung-cancer-iCAP-code.git)
```
### Option 1: Manual install (quick start, may use newer packages)
```r
# Install required packages yourself
install.packages(c("randomForest", "caret", "tidyverse","doParallel", "SummarizedExperiment"))
```
### Option 2: Reproducible environment with renv (recommended)
This project includes an renv.lock file that pins the exact package versions used in the manuscript
(R 4.4.3, renv 1.1.5, Bioconductor 3.20).

To restore that environment:
```r
install.packages("renv")
renv::restore()
```
This will install the precise package versions from the lockfile.

>Note: You don’t have to use renv. If you prefer, follow Option 1 with your own package versions.
>Results may differ slightly if APIs or defaults have changed in newer packages.

## Project Structure
```r
/R-LCiCAPMS2
├── data/                     # Raw and processed data
├── models/                   # model M3 in manuscript
├── R/                        # source files needed to run .rmd files
├── results/                  # Output files
├── LCiCAP_RFM3generation.rmd # r markdown file to reproduce model generation and selection of model M3
├── README.md                 # Project documentation
└── R_LCiCAPMS2.Rproj
```

## Usage
Open the R Project file (R_LCiCAPMS2.Rproj) in R Studio to run interactively.

### Script description:

### 1. LCiCAP_RFM3generation.rmd
This script includes: 
1) Model training by nested cross validation (with 25 seeds and 5 outer folds) 
2) Selection of M3 based on performance (AUC q75)
3) Generation of M3 predictions on the training set (OOB) and predictions on the blind test set.
Analysis of performance on the blind test set (using the predictions and true class labels) was done outside of this script.
Readout was used to produces results shown in Data File 4E and Figure 7 in the manuscript.

To exactly reproduce M3 as reported in the paper, run this file end-to-end, loading the pre-computed results from results/res_3.rds


### Rendering the analysis report
  - **From RStudio:** open `LCiCAP_RFM3generation.Rmd` and click *Knit*.
  - **From the terminal:**  
    Make sure Pandoc is available. You can either:
    - Install Pandoc system-wide (`brew install pandoc` on macOS), or  
    - Point R to RStudio’s bundled Pandoc:
  
  
  ```bash
  export RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools"
  ```
  Then render with:

  ```bash
  Rscript -e 'renv::load("."); rmarkdown::render("LCiCAP_RFM3generation.Rmd", output_format="html_document")'
  ```
  The output LCiCAP_RFM3generation.html will be written next to the Rmd file.


### 2. rf_nested_cv_mc.R
Source file containing the rf_nestedcv_stable() function used in M3 generation.

## Inputs
1. **data/for_modeling/LCII.035.056.hk3-batch_corrected-clean-pos-3.0.20231016.rds**
    - a summarizedExperiment containing
    - 2 assays: NanoString PlexSet gene expression data from all patient serum samples used for modeling (96 training samples + 80 blind test samples)
       - log2 normalized counts
       - normalized counts
       - columns are patient data and rows are expression data + feature metadata
    - patient/sample metadata
       - experimental metada (e.g. icap batch)
       - sample metadata (e.g. hemolysis score)
       - patient metadat (e.g. nodule size)
       - patient status is in column called "unblind" (note that only the 96 training samples are benign or malignant; the 80 blind test samples are labelled "blind")
       - sample partition into train or external blind set is in column called "partition_new".  labels are either "train" or "external_validation"
      

2. **data/stable_genes.csv**
    - List of “stable” genes retained after filtering to remove features showing batch effects between training and blind test batches.
    - Filtering was based on mis-calibration identified from DMOG chemical control replicates run in both batches.

3. **data/blind_testset_unblind_cohort4.xlsx** (not required for running code)
    - Contains unblinded class labels for cohort 4 (samples of the blind tst set) — needed only for evaluating blind set performance.
  

4. **data/nanostring/LCII.035.056.hk3-batch_corrected-clean-pos-3.0.csv**
    - processed NanoString Plexset data in a cvs file. This is not used for the script.

## Outputs
1. **results/preds/231016_model_3_blind_validation_predictions.csv**
    - Predictions from Model M3 on the blind test set.
    - One sample in the blind set is duplicated (two serum samples from same patient, same draw).
    - This sample corresponds to precyte_uids 20230317_09_1 and 20230324_81_1 and patient barcode 620624
    - Performance metrics for Models 3 and 4 are calculated with the two prediction probabilities averaged for this duplicate.  
    

2. **results/preds/231016_model_3_rf_features.csv**
    - Random Forest features for Model M3 with mean decrease accuracy and mean decrease gini.

3. **results/res_3_wrapped_meta.txt**
    - Metadata for M3 (seed, index, selection method).
    - note that seed is NA for M3, but seed number = index number

4. **results/res_3_wrapped.rds**
    - Full results object for M3:
    - Structure:
       - results — list of runs (length = 25 seeds), each containing:
       - final_model — trained randomForest object.
       - auc — mean outer-fold AUC.
       - Additional CV metrics and objects from rf_nestedcv_stable().
    - meta — metadata:
       - chosen_index — run index whose AUC is closest (min abs diff) to the 75th percentile across seeds.
       - chosen_seed — random seed for the chosen run, which is recorded as but chosen seed = chosen index).
       - method_note — description of the selection procedure.  
      
      
5. **results/preds_train/model_3.csv**
    -out of bag predictions of Model M3 on the training set

6. **results/res_3.rds**
    - Original results list (unwrapped, legacy format)

7. **LCiCAP_RFM3generation.html**
    -HTML file from running LCiCAP_RFM3generation.rmd describing model development and shows resulting tables and graphs
      
## Training setup
  - Algorithm: Random Forest (randomForest R package, ntree = 500)
  - Validation: Balanced group 5-fold nested cross-validation
  - Feature selection: none used other than selection of 'stable genes' (see stable_genes.csv for description)
  - Interaction terms: smoking_status allowed only as interaction terms with gene expression features
  - Random seeds: 1–25 for outer runs; internal seeds fixed within each run.

## Model M3 selection
  The model with AUC is closest to the 75th percentile across seeds of nested cross validation
  
## Blind validation
  Model M3 was used to make prediction probabilities of the blind test set. 
  
## Reproducibility
```r
res_3 <- readRDS("results/res_3.rds")
res_3_wrapped <- readRDS("results/res_3_wrapped.rds")
model <- res_3_wrapped$results[[ res_3_wrapped$meta$chosen_index ]]$final_model
```
For exact reproducibility, use the same dataset, preprocessing, and chosen_seed indicated above.
Performance metrics in the manuscript are from this run.

## Dependencies and Environment

This project uses [renv](https://rstudio.github.io/renv/) for reproducibility.  
All package versions are recorded in [`renv.lock`](renv.lock).

- **R version**: 4.4.3 (2025-02-28)
- **Platform**: macOS Ventura 13.0 (Apple Silicon)
- **Bioconductor release**: 3.20
- **Key packages**: `randomForest`, `caret`, `glmnet`, `SummarizedExperiment`, `tidyverse`, `pROC`, `ggplot2`, and others pinned in `renv.lock`.

To recreate the exact environment:

```r
install.packages("renv")
renv::restore()
```
to see all versions used in this project:
```r
renv::dependencies()   # lists packages referenced in code
renv::status()         # shows if library matches lockfile
sessionInfo()          # prints full details for current run
```
>If you prefer not to use renv, you can install packages manually (see Installation), but results may differ slightly with newer package versions.

- **other details:
  -Seeds used: 1:n_seeds (set inside the modeling function)
  -Function used: rf_nestedcv_stable()

## Citation
If you use this model or code, please cite:

Jason D. Berndt, Fergal J. Duff, Mark D. D’Ascenzo, Leslie R. Miller, Yijun Qi, G. Adam 
Whitney, Samuel A. Danziger, Anil Vachani, Pierre P. Massiom, Stephen A. Deppen, 
Robert J. Lipshutz, John D. Aitchison, Jennifer J. Smith, 
**A Multivariate Cell-Based Liquid Biopsy for Lung Nodule Risk Stratification: 
Analytical Validation and Early Clinical Evaluation.**
Journal of Liquid Biopsy, 2025. <a href="https://doi.org/10.1016/j.jlb.2025.100313" target="_blank">DOI</a>


## License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

## Contact
Full pipeline repository: [GitHub repository](https://github.com/jenjsmith/lung-cancer-iCAP-code.git)

For questions, contact [Jennifer Smith](mailto:jennifer@precyte.net)

> This work was supported by the National Institutes of Health under award numbers  
> [R43CA203455], [R44CA203455], [R44AG051282], [P30CA068485], [P30ES013508], [U01CA152662].  
> The content is solely the responsibility of the authors and does not necessarily represent the official views of the NIH.

