EEG Signal Analysis and Regression Modeling

This project focuses on analyzing Electroencephalography (EEG) signals, building regression models to predict an output signal from input signals, and estimating key parameters using Approximate Bayesian Computation (ABC).

Project Goals

•
Explore EEG signals: Understand the characteristics of four input signals (x1-x4) and one output signal (y).

•
Build regression models: Create and compare multiple linear regression models to predict the output signal (y) from the input signals.

•
Select the best model: Use statistical criteria (AIC, BIC, residual analysis) to identify the most effective regression model.

•
Estimate key parameters: Apply Approximate Bayesian Computation (ABC) to determine the posterior distributions of the most influential parameters in the chosen model.

Data

The project uses three simulated data files, all included in this repository:

•
x.csv: Contains the four input EEG signals.

•
y.csv: Contains the single output EEG signal.

•
time.csv: Contains the time vector for the signals.

These files simulate data sampled at 0.002-second intervals over 0.4 seconds, matching the paper's specifications for full reproducibility.

Summary of Analysis and Findings

1. Preliminary Data Analysis

•
Visualization: Signals show volatility and cyclic patterns but no long-term trend. The output signal (y) has sudden spikes, suggesting potential outliers.

•
Univariate Analysis:

•
Input signals (x1-x4) are roughly symmetric and close to normal distributions.

•
The output signal (y) distribution is skewed with heavy tails, indicating potential outliers.



•
Bivariate Analysis:

•
y has a strong positive correlation with x1 and a moderate one with x4.

•
The correlation between y and x2 is weak.



2. Regression Modeling

Five linear regression models were proposed, including polynomial and interaction terms.

Model 3 (y = θ₁x₂ + θ₂x₁³ + θ₃x₄⁴ + θ_bias + ϵ) was selected as the best model. It had the lowest AIC and BIC values, and its residuals were normally distributed with no significant serial correlation.

3. Approximate Bayesian Computation (ABC)

ABC was used to estimate the posterior distributions for θ_bias (intercept) and θ₁ (coefficient for x2) from Model 3.

The expected values of the marginal distributions for each parameter were very close to the original Ordinary Least Squares (OLS) estimates, confirming the initial regression results.

Notes and Corrections

•
File Paths: The analysis.R script now uses relative file paths.

•
Skewness Value: A typo in the paper (page 8) regarding x4's skewness has been noted; the correct value is -0.244741.

•
Random Number Generation: Replaced a custom random number function with R's standard runif for better efficiency.

How to Run the Analysis

1.
Install R Packages: Ensure you have the following R packages installed: tidyverse, GGally, corrplot, ggcorrplot, lmtest, mrfDepth, ggExtra, and MASS. You can install them using install.packages("package_name").

2.
Place Data Files: Download x.csv, y.csv, and time.csv and place them in the same directory as the analysis.R script.

3.
Execute Script: Open analysis.R in R or RStudio and run the entire script. It will display statistical results and generate all plots.

