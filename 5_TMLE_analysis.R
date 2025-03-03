# Load necessary libraries
library(tmle)           # For targeted maximum likelihood estimation
library(SuperLearner)   # For the Super Learner algorithm
library(dplyr)          # For data manipulation
library(ggplot2)        # For visualization

# ------------------------------
# Step 0: Load the Cleaned and Imputed Data
# ------------------------------
# Replace the file path as needed
df <- read.csv("data/imputed_case_study_data.csv")

# Check the first few rows of the dataset
head(df)

# ------------------------------
# Step 1: Define Variables for TMLE
# ------------------------------
# Outcome variable: 'event' (binary indicator for kidney injury)
Y <- df$event

# Treatment variable: 'treatment' (binary: 1 = treated, 0 = untreated)
A <- df$treatment

# Covariates (assumed to be fully imputed)
W <- df %>% subset(., select=c(age, baseline_gfr, bmi, sex, diabetes, hypertension))

# ------------------------------
# Step 2: Specify the Super Learner Library
# ------------------------------
# Define a Super Learner library with a mix of parametric and machine learning algorithms.
# Adjust this library as needed based on data characteristics.
SL.library <- c("SL.glm",      # Standard logistic regression
                "SL.glmnet",   # LASSO/elastic-net regression
                "SL.ranger",   # Random forest (via ranger)
                "SL.xgboost")  # Gradient boosting

# You can add more learners (e.g., "SL.nnet") if desired.
# The same library will be used for estimating both the outcome regression (Q) and the propensity score (g).

# ------------------------------
# Step 3: Implement TMLE
# ------------------------------
# Use the tmle() function to estimate the causal effect.
# The family is set to "binomial" since the outcome is binary.
tmle_fit <- tmle(Y = Y, A = A, W = W,
                 family = "binomial",
                 Q.SL.library = SL.library,
                 g.SL.library = SL.library)

# ------------------------------
# Step 4: Review and Diagnose TMLE Output
# ------------------------------
# Print a summary of the TMLE results
summary(tmle_fit)

# Key output includes:
#   - The estimated risk difference (or other parameter as specified)
#   - The estimated mean outcome under treatment and no treatment (Qbar1 and Qbar0)
#   - Influence curve-based standard errors and confidence intervals

# You may visualize the initial and updated outcome predictions if desired.
# For example, plot the estimated probabilities:
df$Qbar1 <- tmle_fit$Qstar[,2]  # predicted outcome if treated
df$Qbar0 <- tmle_fit$Qstar[,1]  # predicted outcome if untreated

ggplot(df, aes(x = Qbar1)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
  labs(title = "Estimated Outcome Probabilities under Treatment", x = "Predicted Probability", y = "Frequency") +
  theme_minimal()

ggplot(df, aes(x = Qbar0)) +
  geom_histogram(binwidth = 0.05, fill = "tomato", color = "black") +
  labs(title = "Estimated Outcome Probabilities under No Treatment", x = "Predicted Probability", y = "Frequency") +
  theme_minimal()

# ------------------------------
# Step 5: Sensitivity and Additional Diagnostics (Optional)
# ------------------------------
# You can examine the Super Learner weights to see which algorithms contributed most:

#NOTE: Need to fix and learn how to get the weights

print(tmle_fit$g$SL.weights)  # for the treatment model
print(tmle_fit$Q$SL.weights)  # for the outcome model

# Optionally, compare TMLE estimates with those from other methods (e.g., IPW, PSM) to assess robustness.

# End of TMLE Analysis Script
