import math
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.formula.api import ols, logit
from statsmodels.miscmodels.ordinal_model import OrderedModel
from statsmodels.stats.multitest import multipletests
from patsy import dmatrices
import pandas as pd
import numpy as np
import os
from sklearn.metrics import r2_score

# quick plot histogram function
def plot_histograms(df, features_to_plot_list, num_bins=100):
    """Plot histograms in a fixed 3-column layout with dynamic rows."""
    n = len(features_to_plot_list)
    cols = 3
    rows = math.ceil(n / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4.5 * rows))
    axes = axes.flatten()

    for i, col in enumerate(features_to_plot_list):
        axis = axes[i]
        x = df[col].dropna()
        mean, std = x.mean(), x.std()
        axis.hist(x, bins=num_bins)
        axis.set_title(f'{col}\nMean = {round(mean, 2)}\nStd = {round(std, 2)}', fontsize=12)
        axis.tick_params(axis='both', labelsize=10)

    # Remove any unused axes
    for j in range(n, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.show()


# Called in the regression function below for the DV cognate variable covariate functionality
# Added 08/28/2025
# This function was ultimately NOT USED in any of the analyses
def cognate_from(dep: str) -> str | None:
    """
    Return the cognate sleep variable name:
    - if dep contains '_std_' before 'log1p_resid', return the avg version (remove 'std_')
    - otherwise, insert 'std_' before 'log1p_resid' to get the SD version
    Only handles names ending with 'log1p_resid'; returns None otherwise.
    """
    if not isinstance(dep, str) or not dep.endswith("log1p_resid"):
        return None
    if "_std_" in dep:
        return dep.replace("_std_", "_")
    return dep.replace("log1p_resid", "std_log1p_resid")


# Regression function that iterates through sets of variables
def linear_regression_results(dependent_features, independent_features, drop_subset_list, 
                              formula_string, df, test_df, metabolite_anns=False, categorical=False,
                              categorical_poly=False, n_per_predictor=10,
                              model_type='linear', train_test=False, auto_cognate=False, check_iv_levels=True):
    
    """
    Fit regression models across all dependent–independent feature pairs and
    return a tidy results table plus an error log.

    For each dependent/independent feature combination, this function:
    - builds the model formula (optionally adding a cognate covariate),
    - drops rows with missing values in required columns,
    - checks minimum variation, prevalence/level-count thresholds, and
      minimum sample size per predictor,
    - fits a linear, logistic, or ordinal regression model,
    - extracts key model statistics (beta, test statistic, p-value, N, R2/pseudo-R2),
    - optionally evaluates holdout test performance,
    - records failed or skipped models in an error log.

    Supports continuous, categorical, and polynomial-coded categorical
    independent variables, and can optionally merge metabolite annotations
    onto the output table.

    Returns
    -------
    tests : pandas.DataFrame
        Tidy table of successfully fit model results.
    error_log : pandas.DataFrame
        Table describing models that were skipped or failed, with reasons.
    """
    
    # Initialize an empty list to store the results
    all_results = []

    # Initialize error log to capture failed models
    error_log = []

    # Iterate through all combinations of dependent and independent features
    for dependent_feature_input in dependent_features:
        for independent_feature_input in independent_features:

            # --- Added this part 08/28/2025 as a way to include DVs cognate variable as covariate ---
            # build base formula and add cognate to formula automatically
            cognate = cognate_from(dependent_feature_input) if auto_cognate else None
            formula = formula_string.format(
                dependent_feature=dependent_feature_input,
                independent_feature=independent_feature_input
            )

            # ensure we drop NA for any column used in the model
            drop_subset_list_final = [dependent_feature_input, independent_feature_input] + drop_subset_list
            if cognate:
                drop_subset_list_final.append(cognate)
                formula += f" + {cognate}"
            # --- END cognate variable covariate part ---

            # drop the nan from the train and test dfs:
            df_dropped_nan = df.dropna(subset=drop_subset_list_final)
            test_df_dropped_nan = test_df.dropna(subset=drop_subset_list_final)

            # Check that there are at least 2 unique values in the independent
            if df_dropped_nan[independent_feature_input].nunique() < 2:
                error_log.append({
                    "dependent_feature": dependent_feature_input,
                    "independent_feature": independent_feature_input,
                    "formula": formula,
                    "error": "Independent variable has only one unique value",
                    "model_type": model_type
                })
                continue

            # Check that there are at least 2 unique values in the dependent
            if df_dropped_nan[dependent_feature_input].nunique() < 2:
                error_log.append({
                    "dependent_feature": dependent_feature_input,
                    "independent_feature": independent_feature_input,
                    "formula": formula,
                    "error": "Dependent variable has only one unique value",
                    "model_type": model_type
                })
                continue
            
            # Check that there are at least more than 70 non nan
            if len(df_dropped_nan) < 70:
                # skip to next feature in the loop bc N of the final df is too small, will trigger error if
                # allowed to proceed. Just picked 70 = 10*num_predictors, which is number of base predictors
                # ie sex, age, bmi, pc1, pc2, pc3, independent feature
                continue

            if (
                # Check if the independent_feature_input is categorical, binary or multiple categories
                # Then check that each level of independent_feature has count of at least 10
                # If not, then skip to next iteration in loop
                # For binary variable additionally check if proportion of 1s is at least 10%
                df_dropped_nan[independent_feature_input].isin([0, 1]).all()
                or df_dropped_nan[independent_feature_input].dtype in ['object', 'category']
            ):
                level_counts = df_dropped_nan[independent_feature_input].value_counts()
            
                # If binary, additionally check for imbalance of 1s
                if df_dropped_nan[independent_feature_input].isin([0, 1]).all():
                    prop_ones = df_dropped_nan[independent_feature_input].mean()
                    if prop_ones < 0.1 or prop_ones > 0.90:
                        error_log.append({
                            "dependent_feature": dependent_feature_input,
                            "independent_feature": independent_feature_input,
                            "formula": formula,
                            "error": f"1s in INDEPENDENT variable make up {prop_ones:.2%} of the total N (outside 10-90% prevalence)",
                            "model_type": model_type
                        })
                        continue

                if ((level_counts < 10).any()) & (check_iv_levels == True):
                    error_log.append({
                        "dependent_feature": dependent_feature_input,
                        "independent_feature": independent_feature_input,
                        "formula": formula,
                        "error": f"Low level counts in INDEPENDENT variable {level_counts[level_counts < 10].index.tolist()}",
                        "model_type": model_type
                    })
                    continue
            
            # Initialize the model
            if model_type == 'linear':
                model = ols(formula, data=df_dropped_nan)
            
            elif model_type == 'logistic':
                # Check that the 1s make up at least 10% of the total N
                prop_ones = df_dropped_nan[dependent_feature_input].mean()
                if prop_ones < 0.1:
                    error_log.append({
                        "dependent_feature": dependent_feature_input,
                        "independent_feature": independent_feature_input,
                        "formula": formula,
                        "error": f"1s in DEPENDENT variable make up only {prop_ones:.2%} of the total N",
                        "model_type": model_type
                    })
                    continue

                # Make sure that this minimum 10% is at least 10 values
                level_counts = df_dropped_nan[dependent_feature_input].value_counts()
                if (level_counts < 10).any():
                    error_log.append({
                        "dependent_feature": dependent_feature_input,
                        "independent_feature": independent_feature_input,
                        "formula": formula,
                        "error": f"Low level counts in DEPENDENT variable {level_counts[level_counts > 1000].index.tolist()}",
                        "model_type": model_type
                    })
                    continue

                # Fit the model if checks are passed
                model = logit(formula, data=df_dropped_nan)

            elif model_type == 'ordinal':
                # Ultimately did NOT use this part of the function in ANY of the analyses
                # CHECK if there are sparse levels (less than 10 observations) in the dependent ordinal variable
                level_counts = df_dropped_nan[dependent_feature_input].value_counts()
                if (level_counts < 10).any():
                    error_log.append({
                        "dependent_feature": dependent_feature_input,
                        "independent_feature": independent_feature_input,
                        "formula": formula,
                        "error": f"Low level counts in DEPENDENT variable {level_counts[level_counts > 1000].index.tolist()}",
                        "model_type": model_type
                    })
                    continue

                # Fit the model if above check is passed
                # Ensure that the dependent variable is encoded as pandas category and ordered
                df_dropped_nan[dependent_feature_input] = df_dropped_nan[dependent_feature_input].astype('category')
                df_dropped_nan[dependent_feature_input] = df_dropped_nan[dependent_feature_input].cat.as_ordered()

                test_df_dropped_nan[dependent_feature_input] = test_df_dropped_nan[dependent_feature_input].astype('category')
                test_df_dropped_nan[dependent_feature_input] = test_df_dropped_nan[dependent_feature_input].cat.as_ordered()
                
                model = OrderedModel.from_formula(formula, df_dropped_nan, distr='logit')

            else:
                print("No valid model_type parameter passed, skipping to next regression pair")
                continue

            # Find the number of predictors of the initialized model
            num_predictors = model.exog.shape[1]

            # add an additional rule to above bc different formula have different number of predictors
            if (len(df_dropped_nan) >= num_predictors*n_per_predictor):
                # fit the model
                if model_type == 'logistic':
                    try:
                        fitted = model.fit(maxiter=100, disp=0)
                        # If model didn't converge, skip
                        if not fitted.mle_retvals['converged']:
                            error_log.append({
                                "dependent_feature": dependent_feature_input,
                                "independent_feature": independent_feature_input,
                                "formula": formula,
                                "error": "Failed to converge",
                                "model_type": model_type
                            })
                            continue

                    except Exception as e:
                        error_log.append({
                            "dependent_feature": dependent_feature_input,
                            "independent_feature": independent_feature_input,
                            "formula": formula,
                            "error": str(e),
                            "model_type": model_type
                        })
                        continue
                
                elif model_type == 'ordinal':
                    try:
                        fitted = model.fit(maxiter=5000, disp=0)
                        # If model didn't converge skip
                        if not fitted.mle_retvals['converged']:
                            error_log.append({
                                "dependent_feature": dependent_feature_input,
                                "independent_feature": independent_feature_input,
                                "formula": formula,
                                "error": "Failed to converge",
                                "model_type": model_type
                            })
                            continue

                    except Exception as e:
                        error_log.append({
                            "dependent_feature": dependent_feature_input,
                            "independent_feature": independent_feature_input,
                            "formula": formula,
                            "error": str(e),
                            "model_type": model_type
                        })
                        continue
                    
                    # Check that the fitted coefficients for the dependent thresholds are monotonically increasing
                    # Extract threshold coefficients
                    # This is assuming NO OTHER VARIABLE NAMES in the model contain a slash
                    thresholds = fitted.params.filter(like='/')
                    
                    # Convert to array and check monotonicity by subtracting each threshold coefficient from the preceeding one
                    # The coefficients should already be in order
                    threshold_values = thresholds.values
                    is_monotonic = np.all(np.diff(threshold_values) > 0)

                    if not is_monotonic:
                        error_log.append({
                            "dependent_feature": dependent_feature_input,
                            "independent_feature": independent_feature_input,
                            "formula": formula,
                            "error": "Predicted thresholds for ordinal regression do NOT monotonically increase",
                            "model_type": model_type
                        })
                        continue
                
                elif model_type == 'linear':
                    try:
                        fitted = model.fit()

                    # If there's an error fitting the model, skip
                    except Exception as e:
                        error_log.append({
                            "dependent_feature": dependent_feature_input,
                            "independent_feature": independent_feature_input,
                            "formula": formula,
                            "error": str(e),
                            "model_type": model_type
                        })
                        continue

                # NOTE: Still no functionality added here for OrderedModel
                # but doesn't matter for current analyses bc OrderedModel (aka ordinal regression)
                # was not used in ANY of the analyses
                if (train_test == True) and (model_type == 'linear'):
                    # Compute R2 for the fitted model on the holdout test df
                    # Get the residuals and calculate their sum of squares
                    y_test = test_df_dropped_nan[dependent_feature_input]
                    predictions = fitted.predict(test_df_dropped_nan)
                    r2_test = r2_score(y_test, predictions)

                elif (train_test == True) and (model_type != 'linear'):
                    r2_test = fitted.prsquared
                    
                elif train_test == False:
                    y_test = test_df_dropped_nan[dependent_feature_input]
                    r2_test = np.nan
                
                # Write some code to handle a categorical independent variable of interest
                if (categorical == True) & (categorical_poly == False):
                    # Create a Series with the results
                    result_series = pd.Series({
                    "dependent_feature": dependent_feature_input,
                    "independent_feature": independent_feature_input,
                    "n_train": fitted.nobs,
                    "n_test": len(y_test),
                    "r2_train": fitted.rsquared,
                    "r2_test": r2_test,
                    "formula": formula
                    })
                    
                    # Get relevant betas and p-values for categorical variables
                    relevant_params = fitted.params.filter(like=f"C({independent_feature_input}").index
                    relevant_pvalues = fitted.pvalues.filter(like=f"C({independent_feature_input}").index
        
                    # Extract betas
                    model_betas = pd.Series({
                        f"beta_{param.split('[T.')[1][:-1]}": fitted.params[param]
                        for param in relevant_params
                    })
        
                    # Extract p-values and add to test_dict
                    model_pvals = pd.Series({
                        f"pval_{param.split('[T.')[1][:-1]}": fitted.pvalues[param]
                        for param in relevant_pvalues
                    })
    
                    # Extract the model's F-statistic p_value
                    model_f_pval = pd.Series({
                        "f_stat_pval": fitted.f_pvalue
                    })
        
                    # Concatenate the betas and pvals to the result_series
                    result_series = pd.concat([result_series, model_betas, model_pvals, model_f_pval])
    
                    all_results.append(result_series)
                            
                # This code handles categorical_poly independent variables that use the polynomial contrast coding
                elif (categorical == False) & (categorical_poly == True):
                    # Create a Series with the results
                    result_series = pd.Series({
                    "dependent_feature": dependent_feature_input,
                    "independent_feature": independent_feature_input,
                    "beta": fitted.params.filter(regex=f'{independent_feature_input}.*.Linear')[0],
                    "t_statistic": fitted.tvalues.filter(regex=f'{independent_feature_input}.*.Linear')[0],
                    "p": fitted.pvalues.filter(regex=f'{independent_feature_input}.*.Linear')[0],
                    "n_train": fitted.nobs,
                    "n_test": len(y_test),
                    "r2_train": fitted.rsquared if model_type == "linear" else fitted.prsquared,
                    "r2_test": r2_test,
                    "formula": formula
                    })
                    
                    all_results.append(result_series)
                
                elif (categorical == False) & (categorical_poly == False):
                    # Create a Series with the results
                    result_series = pd.Series({
                    "dependent_feature": dependent_feature_input,
                    "independent_feature": independent_feature_input,
                    "beta": fitted.params[independent_feature_input],
                    "t_statistic": fitted.tvalues[independent_feature_input],
                    "p": fitted.pvalues[independent_feature_input],
                    "n_train": fitted.nobs,
                    "n_test": len(y_test),
                    "r2_train": fitted.rsquared if model_type == "linear" else fitted.prsquared,
                    "r2_test": r2_test,
                    "formula": formula
                    })
                    
                    all_results.append(result_series)

            else:
                # skip to next feature in the loop bc N too small
                continue

    # Turn the all_results list and error log into pandas dfs
    tests = pd.DataFrame(all_results)
    error_log = pd.DataFrame(error_log)

    if metabolite_anns == True:
        # add the metabolite metadata
        anns = get_snapshot("metabolomics_metadata").iloc[:, 0:8]
        anns["independent_feature"] = "metabolite_" + anns.CHEMICAL_ID.astype(str)
        tests = pd.merge(tests, anns, on="independent_feature")
        print("This is metabolite regression")

    else:
        print("This is microbe regression")
    
    return tests, error_log
