# Binary_SELECT_PLCO

A repository for the code used for the paper of "A causal inference approach to the external validation of prediction".

## Step 1: Data cleaning
SELECT_clean_2020_function.R

SELECT_clean_2020_manually.R (for checking)

PLCO_clean_2020.R


## Step 2: Model development on PLCO
Prepare_ReadIn.R (pre-processing cleaned data, run it before running the following code)

Table_1.Rmd (characteristics Table 1)

Final_Model.R 

## Step 3: Validate developed model on SELECT 
bootstrap_report.Rmd

Calculation of weighted/unweighted calibration-in-the-large (CIL) and area under the receiver operating characteristics curve (AUC), and the respective 95% confidence intervals. The 95% confidence in serval for unweighted CIL is calculated as shown in the paper, while for the weighted CIL and weighted or unweighted AUC are from the 95% percentiles confidence interval from bootstrapping.
