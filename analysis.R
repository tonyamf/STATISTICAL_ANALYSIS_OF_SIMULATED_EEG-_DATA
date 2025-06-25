# -----------------------------------------------------------------------------
# STATISTICAL ANALYSIS OF SIMULATED EEG DATA
# -----------------------------------------------------------------------------
# This script reproduces the analysis from the paper "7089CEM: Statistical
# Methods Coursework". It covers data loading, preliminary analysis,
# regression modeling, and Approximate Bayesian Computation.
#
# HOW TO USE:
# 1. Install all required packages (listed below).
# 2. Place 'x.csv', 'y.csv', and 'time.csv' in the same directory as this script.
# 3. Run the script.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 1. LOAD REQUIRED PACKAGES
# -----------------------------------------------------------------------------
# Note: Run install.packages("package_name") for any missing packages.

# For data manipulation and plotting (ggplot2)
library(tidyverse)
# For advanced plots like correlation matrices
library(GGally)
# For correlation plotting
library(corrplot)
library(ggcorrplot)
# For statistical tests like Breusch-Godfrey
library(lmtest)
# For the medcouple skewness statistic used in boxplots
library(mrfDepth)
# For plotting marginal histograms/densities on scatter plots
library(ggExtra)


# -----------------------------------------------------------------------------
# 2. DATA LOADING AND PREPARATION
# -----------------------------------------------------------------------------
# The original code used absolute paths. This has been changed to relative paths.
# Please ensure the .csv files are in your R working directory.

# --- Load Time Data ---
# Using try-catch to handle potential file not found errors gracefully.
tryCatch({
  t_raw <- read.csv("time.csv", header = FALSE, stringsAsFactors = FALSE)
  time <- as.numeric(t_raw[,1])
}, error = function(e) {
  stop("Error loading time.csv. Make sure the file is in the working directory.")
})

# --- Load Input Data (x) ---
tryCatch({
  x_raw <- read.csv("x.csv", header = FALSE, stringsAsFactors = FALSE)
  # Convert to a numeric matrix
  input <- data.matrix(x_raw)
}, error = function(e) {
  stop("Error loading x.csv. Make sure the file is in the working directory.")
})


# --- Load Output Data (y) ---
tryCatch({
  y_raw <- read.csv("y.csv", header = FALSE, stringsAsFactors = FALSE)
  output <- data.matrix(y_raw[,1])
}, error = function(e) {
  stop("Error loading y.csv. Make sure the file is in the working directory.")
})


# --- Combine into a single data frame for analysis ---
# This is the main data frame used for most plots and analyses.
df <- data.frame(
  time = time,
  x1 = input[,1],
  x2 = input[,2],
  x3 = input[,3],
  x4 = input[,4],
  y = output[,1]
)

# -----------------------------------------------------------------------------
# 3. TASK 1: PRELIMINARY DATA ANALYSIS
# -----------------------------------------------------------------------------

# --- 3.1 Time Series Plot ---
# Reshape data into a long format for faceting with ggplot
df_long <- df %>%
  select(-time) %>%
  pivot_longer(everything(), names_to = "signal_type", values_to = "value") %>%
  mutate(time = rep(df$time, ncol(df) - 1)) %>%
  # Rename for consistency with the paper's plots
  mutate(signal_type = recode(signal_type,
                              "x1" = "Input x1", "x2" = "Input x2",
                              "x3" = "Input x3", "x4" = "Input x4",
                              "y" = "Output"))

time_series_plot <- ggplot(df_long, aes(x = time, y = value, colour = signal_type)) +
  geom_line() +
  facet_grid(signal_type ~ ., scales = "free_y") +
  theme_bw() +
  guides(color = "none") +
  labs(title = "Time Series of Input and Output EEG Signals", x = "Time (s)", y = "Signal Value")

print(time_series_plot)


# --- 3.2 Univariate Analysis: Box Plots ---
# Using the same long data frame
box_plot <- ggplot(df_long, aes(y = value, x = signal_type, colour = signal_type)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 16) +
  facet_wrap(~ signal_type, scales = "free") +
  theme_bw() +
  guides(color = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(title = "Distribution of Each EEG Signal", x = "Signal", y = "Signal Value")

print(box_plot)


# --- 3.3 Univariate Analysis: Histograms and Density Plots ---
# Calculate means for vertical lines on the plot
means <- df_long %>%
  group_by(signal_type) %>%
  summarise(mean_val = mean(value))

hist_plot <- ggplot(df_long, aes(x = value)) +
  geom_histogram(aes(y = ..density..), binwidth = 1, colour = "black", fill = "white") +
  geom_density(aes(color = signal_type), alpha = 0.2, fill = "#FF6666") +
  geom_vline(data = means, aes(xintercept = mean_val, color = signal_type), linetype = "dashed", size = 1) +
  facet_wrap(~ signal_type, scales = "free") +
  theme_bw() +
  guides(color = "none") +
  labs(title = "Histogram and Density Plot for Each Signal",
       subtitle = "Dashed line indicates the mean",
       x = "Signal Value", y = "Density")

print(hist_plot)


# --- 3.4 Descriptive Statistics ---
# Function to calculate skewness
skewness_fn <- function (x){
  (sqrt(length(x)) * sum( (x - mean(x))^3 ) ) / (sqrt( sum( (x - mean(x))^2 ) ))^3
}

# Function to calculate kurtosis
kurtosis_fn <- function(data){
  (sum((data - mean(data))^4) / length(data)) / ((sum((data - mean(data))^2) / length(data))^2)
}

# Print summary statistics
cat("\n--- Summary Statistics ---\n")
print(summary(df))

cat("\n--- Standard Deviations ---\n")
print(sapply(df, sd))

cat("\n--- Skewness ---\n")
print(sapply(df, skewness_fn))

cat("\n--- Kurtosis ---\n")
print(sapply(df, kurtosis_fn))


# --- 3.5 Bivariate Analysis: Correlation Plot ---
# Custom function to add both lm and loess smoothers to scatter plots
my_scatter_fn <- function(data, mapping, ...){
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = loess, formula = 'y ~ x', color = "red", se = FALSE, ...) +
    geom_smooth(method = lm, formula = 'y ~ x', color = "blue", se = TRUE, ...)
}

cat("\n--- Generating Correlation Matrix Plot (this may take a moment) ---\n")
correlation_plot <- ggpairs(df,
                            lower = list(continuous = my_scatter_fn),
                            diag = list(continuous = "blankDiag")) +
  labs(title = "Correlation Matrix of EEG Signals",
       subtitle = "Blue line: Linear Model, Red line: LOESS Smoother") +
  theme_bw()

print(correlation_plot)


# -----------------------------------------------------------------------------
# 4. TASK 2: REGRESSION ANALYSIS
# -----------------------------------------------------------------------------

# --- 4.1 Define the 5 Models ---
ones <- matrix(1, nrow(df), 1)

# Model 1: 𝑦 = 𝜃1*𝑥4 + 𝜃2*𝑥1^2 + 𝜃3*𝑥1^3 + 𝜃4*𝑥3^4 + 𝜃_bias + 𝜀
model1Input <- cbind(ones, df$x4, df$x1^2, df$x1^3, df$x3^4)

# Model 2: 𝑦 = 𝜃1*𝑥3^3 + 𝜃2*𝑥3^4 + 𝜃_bias + 𝜀
model2Input <- cbind(ones, df$x3^3, df$x3^4)

# Model 3: 𝑦 = 𝜃1*𝑥2 + 𝜃2*𝑥1^3 + 𝜃3*𝑥3^4 + 𝜃_bias + 𝜀  (Chosen as best model)
model3Input <- cbind(ones, df$x2, df$x1^3, df$x3^4)

# Model 4: 𝑦 = 𝜃1*𝑥4 + 𝜃2*𝑥1^3 + 𝜃3*𝑥3^4 + 𝜃_bias + 𝜀
model4Input <- cbind(ones, df$x4, df$x1^3, df$x3^4)

# Model 5: 𝑦 = 𝜃1*𝑥4 + 𝜃2*𝑥1^2 + 𝜃3*𝑥1^3 + 𝜃4*𝑥3^4 + 𝜃5*𝑥1^4 + 𝜃_bias + 𝜀
model5Input <- cbind(ones, df$x4, df$x1^2, df$x1^3, df$x3^4, df$x1^4)


# --- 4.2 Master Function for Regression Analysis ---
# This function calculates all the required metrics for a given model.
analyze_regression_model <- function(X_train, y_train, model_name) {
  # Estimate parameters using Ordinary Least Squares (OLS)
  # Using try-catch for singularity issues
  Theta_hat <- tryCatch({
    solve(t(X_train) %*% X_train) %*% t(X_train) %*% y_train
  }, error = function(e) {
    stop(paste("Matrix is singular for", model_name, ". Cannot compute inverse."))
  })

  # Predict values
  y_hat <- X_train %*% Theta_hat
  residuals <- y_train - y_hat
  n <- nrow(X_train)
  k <- ncol(X_train) # Number of parameters

  # Calculate metrics
  rss <- sum(residuals^2)
  sigma_2_hat <- rss / (n - k)
  log_likelihood <- -n/2 * log(2 * pi) - n/2 * log(sigma_2_hat) - rss / (2 * sigma_2_hat)
  aic <- 2 * k - 2 * log_likelihood
  bic <- k * log(n) - 2 * log_likelihood

  # Create a Q-Q plot of residuals
  qqplot_data <- data.frame(residuals = as.vector(residuals))
  qqplot <- ggplot(qqplot_data, aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line() +
    theme_bw() +
    ggtitle(paste("Normal Q-Q Plot for", model_name, "Residuals"))

  # Return results
  list(
    model_name = model_name,
    coefficients = Theta_hat,
    rss = rss,
    log_likelihood = log_likelihood,
    aic = aic,
    bic = bic,
    residuals = residuals,
    qqplot = qqplot
  )
}

# --- 4.3 Analyze All Models ---
cat("\n--- Analyzing Regression Models ---\n")
model1_results <- analyze_regression_model(model1Input, df$y, "Model 1")
model2_results <- analyze_regression_model(model2Input, df$y, "Model 2")
model3_results <- analyze_regression_model(model3Input, df$y, "Model 3")
model4_results <- analyze_regression_model(model4Input, df$y, "Model 4")
model5_results <- analyze_regression_model(model5Input, df$y, "Model 5")

results_list <- list(model1_results, model2_results, model3_results, model4_results, model5_results)

# --- 4.4 Display Model Comparison Table ---
comparison_df <- data.frame(
  Model = sapply(results_list, function(x) x$model_name),
  RSS = sapply(results_list, function(x) x$rss),
  LogLikelihood = sapply(results_list, function(x) x$log_likelihood),
  AIC = sapply(results_list, function(x) x$aic),
  BIC = sapply(results_list, function(x) x$bic)
)

cat("\n--- Model Comparison Table ---\n")
print(comparison_df)

cat("\nModel 3 has the lowest AIC and BIC, suggesting it's the best model.\n")

# --- 4.5 Display Residual Plots and Tests for Best Model (Model 3) ---
cat("\n--- Residual Analysis for Best Model (Model 3) ---\n")
print(model3_results$qqplot)

# Breusch-Godfrey test for serial correlation
# We create a temporary data frame for the lm() function
temp_df_model3 <- as.data.frame(model3Input)
bg_test_model3 <- bgtest(lm(df$y ~ . -1, data = temp_df_model3), order = 1)
cat("\nBreusch-Godfrey Test for Model 3 (p > 0.05 indicates no serial correlation):\n")
print(bg_test_model3)


# -----------------------------------------------------------------------------
# 5. TASK 2.7: MACHINE LEARNING PREDICTION (TRAIN/TEST SPLIT)
# -----------------------------------------------------------------------------
cat("\n--- Training Model 3 on 70% of data and testing on 30% ---\n")

set.seed(42) # for reproducibility
train_size <- floor(0.7 * nrow(df))
train_indices <- 1:train_size

# Split data
X_train <- model3Input[train_indices, ]
y_train <- df$y[train_indices]
X_test <- model3Input[-train_indices, ]
y_test <- df$y[-train_indices]
x_values_test <- df[-train_indices, ] # For plotting

# Retrain model on training data
trained_model <- analyze_regression_model(X_train, y_train, "Trained Model 3")

# Predict on test data
y_pred <- X_test %*% trained_model$coefficients
test_residuals <- y_test - y_pred
test_rss <- sum(test_residuals^2)
r_squared <- 1 - (test_rss / sum((y_test - mean(y_test))^2))

cat(paste("\nTest R-squared for Model 3:", round(r_squared, 4), "\n"))

# --- Confidence Intervals for Predictions ---
# Function to calculate confidence intervals for predictions
calculate_ci <- function(X, sigma_2_hat, confidence_level = 0.95) {
  z_score <- qnorm(1 - (1 - confidence_level) / 2)
  cov_theta_hat <- sigma_2_hat * solve(t(X) %*% X)
  var_y_hat <- diag(X %*% cov_theta_hat %*% t(X))
  ci_margin <- z_score * sqrt(var_y_hat)
  return(ci_margin)
}

# Calculate CIs on the test set predictions
n_train <- nrow(X_train)
k_train <- ncol(X_train)
sigma_2_hat_train <- sum(trained_model$residuals^2) / (n_train - k_train)
ci_margin <- calculate_ci(X_test, sigma_2_hat_train)

# Plot predictions with confidence intervals
pred_df <- data.frame(
  x2 = x_values_test$x2,
  y_true = y_test,
  y_pred = as.vector(y_pred),
  lower_ci = as.vector(y_pred) - ci_margin,
  upper_ci = as.vector(y_pred) + ci_margin
)

prediction_plot <- ggplot(pred_df, aes(x = x2)) +
  geom_point(aes(y = y_true), color = "blue", alpha = 0.6, shape = 1) +
  geom_point(aes(y = y_pred), color = "red") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), color = "red", width = 0.1, alpha = 0.5) +
  theme_bw() +
  labs(
    title = "Model 3 Predictions on Test Data with 95% Confidence Intervals",
    subtitle = "Blue Points: Actual Values, Red Points/Bars: Predicted Values ± CI",
    x = "Input x2 Value (from Test Set)",
    y = "Output Signal y"
  )
print(prediction_plot)


# -----------------------------------------------------------------------------
# 6. TASK 3: APPROXIMATE BAYESIAN COMPUTATION (ABC)
# -----------------------------------------------------------------------------
cat("\n--- Performing Approximate Bayesian Computation for Model 3's largest coefficients ---\n")

# Parameters of interest from full model (Model 3)
# Coeffs are: bias, theta1 (for x2), theta2 (for x1^3), theta3 (for x3^4)
# We focus on bias (coeff 1) and theta1 (coeff 2)
theta_hat_full <- model3_results$coefficients
params_of_interest <- theta_hat_full[1:2]

# Calculate standard errors for the priors
n_full <- nrow(model3Input)
k_full <- ncol(model3Input)
sigma_2_hat_full <- model3_results$rss / (n_full - k_full)
cov_theta_hat_full <- sigma_2_hat_full * solve(t(model3Input) %*% model3Input)
se_theta <- sqrt(diag(cov_theta_hat_full))
se_of_interest <- se_theta[1:2]

# --- Setup ABC ---
n_sim <- 100000
# Using R's built-in runif is more standard and efficient
# Priors are uniform, centered around the OLS estimate +/- 3 standard errors
prior_samples_bias <- runif(n_sim, params_of_interest[1] - 3 * se_of_interest[1], params_of_interest[1] + 3 * se_of_interest[1])
prior_samples_theta1 <- runif(n_sim, params_of_interest[2] - 3 * se_of_interest[2], params_of_interest[2] + 3 * se_of_interest[2])

# Set a tolerance epsilon
# A reasonable tolerance is a small fraction of the original RSS
tolerance <- 0.5 # As distance metric, based on paper's logic of RSS diff

# Store accepted samples
accepted_samples <- data.frame(bias = numeric(), theta1 = numeric())

# --- Run ABC Rejection Algorithm ---
cat("Running ABC loop... this may take some time.\n")
for (i in 1:n_sim) {
  # Current parameter proposal from priors
  current_theta <- c(prior_samples_bias[i], prior_samples_theta1[i])

  # Create full parameter vector for prediction
  theta_vector_sim <- c(current_theta, theta_hat_full[3], theta_hat_full[4])

  # Simulate data (i.e., calculate predicted y)
  y_sim <- model3Input %*% theta_vector_sim

  # Calculate summary statistic (RSS in this case)
  rss_sim <- sum((df$y - y_sim)^2)

  # Compare with observed data's summary statistic
  distance <- abs(rss_sim - model3_results$rss)

  # Accept or reject
  if (distance < tolerance) {
    accepted_samples <- rbind(accepted_samples, data.frame(bias = current_theta[1], theta1 = current_theta[2]))
  }
}

cat(paste("Finished ABC. Accepted", nrow(accepted_samples), "samples out of", n_sim, "\n"))

# --- Plot Posterior Distributions ---
# Joint distribution plot
joint_plot <- ggplot(accepted_samples, aes(x = bias, y = theta1)) +
  geom_density_2d_filled(alpha = 0.7) +
  theme_bw() +
  guides(fill = "none") +
  labs(
    title = "ABC Joint Posterior Distribution",
    x = "Posterior for Theta_Bias",
    y = "Posterior for Theta_1 (for x2)"
  )

# Add marginal densities
final_abc_plot <- ggMarginal(joint_plot, type = "densigram", fill = "lightblue")

print(final_abc_plot)

# --- Posterior Summary ---
cat("\n--- ABC Posterior Summary Statistics ---\n")
summary(accepted_samples)

cat("\n--- OLS Estimates for Comparison ---\n")
print(params_of_interest)

