# Binary_SELECT_PLCO

A reprository for the code used for the binary outcome causal inference analysis.

## Step 1: Data cleaning
SELECT_clean_2020_function.R
SELECT_clean_2020_manually.R (for checking)
PLCO_clean_2020.R

Table_1.Rmd (characteristics Table 1)

## Step 2: Model development on PLCO
Prepare_ReadIn.R (pre-processing cleaned data, run it before running Final_Models.R)
Final_Models.R 

## Step 3: Validate developed model on SELECT 
bootstrap_report.Rmd
calculation of weighted/unweighted calibration-in-the-large (CIL) and area under the receiver operating characteristics curve (AUC), and the respective 95% confidence intervals. The 95% confidence in serval for unweighted CIL is calculated as the paper, while for the weighted CIL and weighted or unweighted AUC are from the 95% percentiles confidence interval from bootstrapping.
