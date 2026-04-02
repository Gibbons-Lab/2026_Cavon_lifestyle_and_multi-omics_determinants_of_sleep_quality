# Project Structure and Organization

This repository is organized into subdirectories that contain the cohort dataset subsetting (generating the master working_df), analyses, intermediate outputs, and figure generation for this study. Each directory corresponds to a specific component of the project as described below.

## Directory Overview
---
### Subsetted Cohort Working Dataset
The master analysis dataset is constructed in two stages:
  1. `fitbit_data_processing`: cleaning and aggregation of Fitbit-derived metrics.
  2. `working_df`: integration of all additional data features onto the processed Fitbit dataset.

#### `fitbit_data_processing`
Contains the workflow used to generate the aggregated Fitbit-derived sleep and activity features (i.e., monthly means and variability metrics).

#### `working_df`
Contains the notebook used to merge questionnaire data, microbiome features, metabolomics, proteomics, clinical chemistries, and other covariates onto the processed Fitbit-derived dataset generated in `fitbit_data_processing`. Together, these steps produce the final merged dataframe used for all downstream analyses.

---

### `cohort_statistics_and_covariate_correlations`
Summary statistics of the subsetted cohort and correlation analyses between covariates used across models.

**Associated manuscript figures/data:**
- Figure 2b
- Figure S1

---

### `activity_sleep_happiness_questionnaire_alignment_with_fitbit`
Analyses and visualizations comparing Fitbit-derived activity and sleep metrics with questionnaire-based measures, including Oxford Happiness Questionnaire (OHQ) alignment.

**Associated manuscript figures/data:**
- Figure 1

---

### `random_forest_partitioning_variance_analysis`
Random forest analyses used to quantify variance explained by different feature sets.

**Associated manuscript figures/data:**
- Figure 2a
- Figure S2
- Table S1

---

### `lifestyle_and_multi-omics_on_sleep_regression_analysis`
Primary regression analyses linking lifestyle factors and multi-omic features to sleep outcomes, along with associated visualizations.

**Associated manuscript figures/data:**
- Figures 3-7
- Figure S3

---

### `microbe_activity_interaction_effects_regression_analysis`
Regression models assessing interaction effects between gut microbial presence/absence features and physical activity on sleep outcomes.

**Associated manuscript figures/data:**
- Figure 8
- Figure S4

---

### `regression_and_plotting_functions`
Reusable functions for regression modeling used across notebooks.

---
