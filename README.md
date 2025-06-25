Statistical Analysis of Simulated EEG Time-Series Data
This project analyzes simulated EEG time-series data to reproduce the statistical findings of a research paper. It involves three main stages:

Project Goals

Explore EEG signals: Understand the characteristics of four input signals (x1-x4) and one output signal (y).
Build regression models: Create and compare multiple linear regression models to predict the output signal (y) from the input signals.
Select the best model: Use statistical criteria (AIC, BIC, residual analysis) to identify the most effective regression model.
Estimate key parameters: Apply Approximate Bayesian Computation (ABC) to determine the posterior distributions of the most influential parameters in the chosen model.


Data
The project uses three simulated data files, all included in this repository:
x.csv: Contains the four input EEG signals.
y.csv: Contains the single output EEG signal.
time.csv: Contains the time vector for the signals.

These files simulate data sampled at 0.002-second intervals over 0.4 seconds, matching the paper's specifications for full reproducibility.


How to Run the Analysis
Install R Packages: Ensure you have the following R packages installed: tidyverse, GGally, corrplot, ggcorrplot, lmtest, mrfDepth, ggExtra, and MASS. You can install them using install.packages("package_name").
Place Data Files: Download x.csv, y.csv, and time.csv and place them in the same directory as the analysis.R script.
Execute Script: Open analysis.R in R or RStudio and run the entire script. It will display statistical results and generate all plots.
