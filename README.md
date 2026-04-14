# Lifestyle and Multi-Omic Determinants of Sleep Quality

## Project Structure and Organization

This repository is organized into subdirectories that contain the subsetted cohort dataset, analyses, intermediate outputs, and figure generation for this study. Each directory corresponds to a specific component of the project as described below.

## Directory Overview
---
### Subsetted Cohort Working Dataset
The master analysis dataset was constructed in two stages:
  1. `fitbit_data_processing`: cleaning and aggregation of Fitbit-derived metrics.
  2. `working_df`: integration of all additional data features onto the processed Fitbit dataset.

#### `fitbit_data_processing`
Contains the workflow used to generate the aggregated Fitbit-derived sleep and activity features (i.e., monthly means and variability metrics).

#### `working_df`
Contains the notebook used to merge questionnaire data, microbiome features, metabolomics, proteomics, clinical chemistries, and other covariates onto the processed Fitbit-derived dataset generated in `fitbit_data_processing`. Together, these steps produce the final merged dataframe used for all downstream analyses.

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

### `cohort_statistics_and_covariate_correlations`
Summary statistics of the subsetted cohort and correlation analyses between covariates used across models.

**Associated manuscript figures/data:**
- Figure 2b
- Figure S1

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

## Computational Environment
All analyses were performed on a Linux-based high-performance computing environment (x86_64 architecture) with the following specifications:
- CPU: 40 cores
- RAM: ~1 TB

The analyses and data visualizations were executed using a Jupyter Notebook environment and are expected to be compatible with standard Linux and macOS systems.

## Software Requirements

---

The following software environments and packages were used for the analyses:
- Python (v3.11.5)
  #### Fitting individual-level microbial community-scale metabolic models
  - MICOM (v0.37.0)
  #### Random forest analysis `random_forest_partitioning_variance_analysis`
  - sklearn (v1.3.1)
  #### Main effects and interaction effects regression model specification and fitting
  - statsmodels (v0.14.0)
  - patsy (v0.5.3)
  #### All main-text figures except for Figure 7
  - altair (v5.0.1)
  - matplotlib (v3.8.0)
  #### Figure 7
  - plotly (v5.17.0)
  #### Data wrangling
  - pandas (v2.1.1)
  - numpy (v1.26.0)

---

- R (v4.4.2)
  #### Rarefying microbiome count table and computing microbiome alpha diversities
  - vegan (v2.7.2)
  #### Intermediate data visualizations for generating microbiome rarefied count table
  - ggplot2 (v4.0.1)
  #### Data wrangling
  - tidyr (v1.3.2)
  - dplyr (v1.1.4)

---

## Installation Guide

### Instructions

To reproduce the analyses, clone this repository and install the required Python and R dependencies listed above.

#### Python environment
We recommend using a virtual environment:

```bash
git clone https://github.com/Gibbons-Lab/2026_Cavon_lifestyle_and_multi-omic_determinants_of_sleep_quality.git
cd 2026_Cavon_lifestyle_and_multi-omic_determinants_of_sleep_quality
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

#### R environment
Install required R packages:

```r
install.packages(c("vegan", "ggplot2", "tidyr", "dplyr"))
```

### Typical Install Time

Installation of all dependencies typically takes approximately **15 minutes** on a standard desktop computer, depending on internet speed and environment configuration.

---

## Demo

### Instructions to Run Example Code

Due to data availability restrictions, the dataset used for analyses can not be distributed with this repository. However, intermediate output data files are provided.

To run a data visualization of a representative analysis:

1. Navigate to:
   `microbe_activity_interaction_effects_regression_analysis/`

2. Open and execute:
   `microbe_activity_interaction_effects_regression_analysis.ipynb`

As outlined above, this notebook contains the code to conduct the interaction effects analysis (Figure 8), and has all the necesarry output files to generate the data visualization.

### Expected Output

Running the notebook will generate:
- Regression model outputs for activity × microbe interaction effects
- Corresponding visualizations (Figure 8)
- Intermediate result tables saved within the analysis directory

### Expected Run Time

On a standard desktop computer (4 cores, 16 GB RAM), this demo is expected to run in approximately **10 minutes**.

---

## Instructions for Use

### Running the Full Analysis Pipeline

The full analysis pipeline consists of the following steps:

1. **Process Fitbit data**
   - Run the notebook in `fitbit_data_processing/` to clean and aggregate raw Fitbit-derived sleep and activity data.

2. **Construct working dataset**
   - Run the notebook in `working_df/` to integrate questionnaire data, microbiome features, metabolomics, proteomics, clinical chemistries, and covariates.

3. **Perform analyses**
   - Activity–sleep alignment and questionnaire comparisons:
     `activity_sleep_happiness_questionnaire_alignment_with_fitbit/`
   - Random forest variance partitioning:
     `random_forest_partitioning_variance_analysis/`
   - Cohort statistics and covariate correlations:
     `cohort_statistics_and_covariate_correlations/`
   - Main regression analyses:
     `lifestyle_and_multi-omics_on_sleep_regression_analysis/`
   - Interaction effect models:
     `microbe_activity_interaction_effects_regression_analysis/`

4. **Generate figures**
   - Each analysis directory contains notebooks that reproduce the corresponding manuscript figures and supplementary materials as described above.

---

## Data Availability

The datasets used in this study are not publicly available due to restrictions related to participant privacy and data use agreements.

Access to the underlying data may be requested through the appropriate data access channels associated with the original study (e.g., the Arivale cohort), subject to approval (see `Data Availability` statement in Manuscript for further details).

Where permitted, intermediate output data files required to reproduce figures and analyses are provided within this repository.
