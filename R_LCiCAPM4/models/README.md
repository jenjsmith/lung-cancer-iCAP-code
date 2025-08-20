# `res_3.rds` and `res_4.rds` – Final Model Ensembles from Stability-Selected Nested CV 
*(corresponding to M3 and M4 in the manuscript, respectively)*

These two files (`res_3.rds` and `res_4.rds`) contains trained classification models developed using stability-based feature selection and nested cross-validation on NanoString Plexset gene expression data. The models were generated using the `rf_nestedcv_stable()` function in R and are intended for use in lung nodule risk stratification or related diagnostic modeling.
---
## Contents

The file is an R list object containing results from multiple random seeds:

Each entry `res_3[[i]]` or `res_4[[i]]`includes:

- `accuracy`: Average classification accuracy on outer test folds
- `auc`: Area under the ROC curve
- `stable_features`: Features selected as stable (i.e., high-frequency across subsamples)
- `final_model`: The final trained `randomForest` model using those stable features
- `outer_res`: Performance metrics and models from each outer fold
- `preprocess_values`: Centering/scaling information (if used)
- `error_code`: Indicates success (0) or failure

The final models (M3 and M4) referenced in the manuscript were selected as the models with AUC values closest to the 75th percentile of all AUCs within res_3 and res_4, respectively.
---
## Reproducibility


To load and extract the final model from `res_3.rds`:

```r
res_3 <- readRDS("res_3.rds")

# Extract AUCs
aucs <- sapply(res_3, function(x) x$auc)

# Select the model closest to 75th percentile AUC
q75 <- quantile(aucs, 0.75)
model_3 <- res_3[[which.min(abs(aucs - q75))]]$final_model
```
Repeat similarly for res_4.rds.
## Environment

- R version: 4.4.3 (2025-02-28)
- Platform: macOS Ventura 13.0 (Apple Silicon)
- Key packages:
  - `randomForest` 4.7-1.2
  - `glmnet` 4.1-8
  - `caret` 7.0-1
  - `SummarizedExperiment` 1.36.0
  - `tidyverse` 2.0.0
- other details:
  -Seeds used: 1:n_seeds (set inside the modeling function)
  -Function used: rf_nestedcv_stable()

## Related Files
- rf_nestedcv_stable.R: Function used to generate the models
-stable_genes.csv: List of input genes used for modeling
-SummarizedExperiment object: Preprocessed gene expression data

# Citation
If you use this model or code, please cite:

Jason D. Berndt, Fergal J. Duff, Mark D. D’Ascenzo, Leslie R. Miller, Yijun Qi, G. Adam 
Whitney, Samuel A. Danziger, Anil Vachani, Pierre P. Massiom, Stephen A. Deppen, 
Robert J. Lipshutz, John D. Aitchison, Jennifer J. Smith, 
*A Multivariate Cell-Based Liquid Biopsy for Lung Nodule Risk Stratification: 
Analytical Validation and Early Clinical Evaluation.
Journal of Liquid Biopsy, 20225. DOI: DOI: 10.1016/j.jlb.2025.100313 *

## Contact
Full pipeline repository: [GitHub repository](https://github.com/jenjsmith/lung-cancer-iCAP-code.git)

For questions, contact [Jennifer Smith](mailto:jennifer@precyte.net)

> This work was supported by the National Institutes of Health under award numbers  
> [R43CA203455], [R44CA203455], [R44AG051282], [P30CA068485], [P30ES013508], [U01CA152662].  
> The content is solely the responsibility of the authors and does not necessarily represent the official views of the NIH.


