Statistical Analysis of Simulated EEG Time-Series Data
Project Overview
This project reproduces the statistical analysis presented in the paper. The analysis investigates a simulated EEG dataset comprising four input time-series signals (x1, x2, x3, x4) and one output time-series signal (y).

The core objectives of the analysis are:

To perform a preliminary exploratory data analysis (univariate and bivariate) to understand the characteristics of the EEG signals.

To build and compare several multiple linear regression models to predict the output signal y based on the input signals.

To select the best regression model using statistical criteria like AIC, BIC, and residual analysis.

To use Approximate Bayesian Computation (ABC) to estimate the posterior distributions for the most influential parameters of the chosen model.

This repository contains the R script (analysis.R) to replicate the entire analysis and generate all the plots and results shown in the paper, along with the necessary data files.

Data Source
The analysis requires three data files, which are included in this repository:

x.csv: A CSV file containing the 4 input EEG signals.

y.csv: A CSV file containing the single output EEG signal.

time.csv: A CSV file containing the time vector for the signals.

The paper states this is simulated data, sampled at equally spaced intervals of 0.002 seconds over a total duration of 0.4 seconds. The data files provided here have been simulated to match the statistical properties described in the paper, allowing for a full reproduction of the analysis.

Summary of Analysis and Findings
1. Preliminary Data Analysis (Task 1)
Time Series Visualization: The signals are plotted over time. They exhibit volatility and cyclic patterns but no clear long-term trend. The output signal shows some sudden spikes, hinting at potential outliers.

Univariate Analysis:

Distributions of the input signals (x1 to x4) are roughly symmetric and close to a normal distribution.

The output signal distribution is skewed and has heavy tails, indicating potential outliers.

The time variable follows a continuous uniform distribution, as expected.

Bivariate Analysis:

A correlation matrix reveals the linear relationships between all variables.

The output y has a strong positive linear correlation with input x1 and a moderate one with x4.

The correlation between y and x2 is weak.

2. Regression Modeling (Task 2)
Five different linear regression models were proposed to predict y. These models include polynomial and interaction terms of the input variables.

Model Selection:
The models were evaluated based on several criteria:

Residual Sum of Squares (RSS): Measures the total prediction error.

Log-Likelihood: Measures the goodness of fit.

Akaike Information Criterion (AIC) & Bayesian Information Criterion (BIC): These are penalized-likelihood criteria that balance model fit with model complexity to prevent overfitting. Lower values are better.

Residual Analysis: The residuals of the models were checked for normality (using Q-Q plots) and serial correlation (using the Breusch-Godfrey test).

Conclusion:
Model 3 (y = θ1*x2 + θ2*x1^3 + θ3*x3^4 + θ_bias + ε) was chosen as the best model. It had the lowest AIC and BIC values, and its residuals closely followed a normal distribution with no significant serial correlation, satisfying the key modeling assumptions.

3. Approximate Bayesian Computation (ABC) (Task 3)
The final task focused on the two coefficients from Model 3 with the largest absolute values: θ_bias (the intercept) and θ1 (the coefficient for x2).

Method: A rejection ABC algorithm was used to approximate the posterior distribution of these two parameters.

Prior: Uniform priors were defined for each parameter, centered around their Ordinary Least Squares (OLS) estimates from Task 2.

Result: The analysis generated a joint posterior distribution, showing the most plausible combinations of the two parameters given the data. The marginal distributions for each parameter were found to have expected values very close to the original OLS estimates, confirming the initial regression results.

Notes and Corrections
File Paths: The original R code used hardcoded absolute file paths. The provided analysis.R script has been modified to use relative paths, assuming the data files are in the same directory.

Skewness Value: There is a discrepancy in the paper. On page 8, the skewness for input x4 is reported as 0.2062838. However, the R code's output on page 31 correctly calculates it as -0.244741. The text in the paper likely contains a typo.

Random Number Generation: The original code implemented a custom function myrunif for generating uniformly distributed random numbers. This has been replaced with R's standard and more efficient runif function in the corrected script.

How to Run the Analysis
Install R Packages: Make sure you have the following R packages installed. You can install them using the command install.packages("package_name").

tidyverse (for data manipulation and ggplot2)

GGally (for the correlation matrix plot)

corrplot and ggcorrplot (for correlation plots)

lmtest (for the Breusch-Godfrey test)

mrfDepth (for the medcouple skewness statistic)

ggExtra (for marginal plots)

MASS (used internally for data simulation)

Place Data Files: Download the x.csv, y.csv, and time.csv files and place them in the same folder as the analysis.R script. (They are provided here).

Execute the Script: Open analysis.R in R or RStudio and run the entire script. It will print the statistical results to the console and generate all the plots from the paper.
