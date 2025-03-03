

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(naniar)       # For missing data visualization
library(tableone)     # For summary tables
library(missForest)   # For random-forest based imputation

df <- read.csv("data/simulated_case_study_data.csv")


# ------------------------------
# Step 1: Exploratory Data Tabulation & Visualization
# ------------------------------

# 1.1 Visualize missingness patterns across variables
vis_miss(df)

# 1.2 Create a summary table of key covariates (including missingness information)
covariate_vars <- c("age", "baseline_gfr", "bmi", "sex", "diabetes", "hypertension")
table1 <- CreateTableOne(vars = covariate_vars, data = df, test = FALSE)
print(table1, missing = TRUE)

# 1.3 Plot histograms for continuous covariates
# Histogram for Age
ggplot(df, aes(x = age)) +
  geom_histogram(binwidth = 2, fill = "skyblue", color = "black", na.rm = TRUE) +
  labs(title = "Histogram of Age (non-missing)", x = "Age", y = "Frequency") +
  theme_minimal()

# Histogram for Baseline GFR
ggplot(df, aes(x = baseline_gfr)) +
  geom_histogram(binwidth = 2, fill = "lightgreen", color = "black", na.rm = TRUE) +
  labs(title = "Histogram of Baseline GFR (non-missing)", x = "Baseline GFR", y = "Frequency") +
  theme_minimal()

# Histogram for BMI
ggplot(df, aes(x = bmi)) +
  geom_histogram(binwidth = 1, fill = "lightcoral", color = "black", na.rm = TRUE) +
  labs(title = "Histogram of BMI (non-missing)", x = "BMI", y = "Frequency") +
  theme_minimal()

# 1.4 Scatter plot example: Age vs. Baseline GFR, colored by diabetes status
ggplot(df, aes(x = age, y = baseline_gfr, color = factor(diabetes))) +
  geom_point(alpha = 0.6, na.rm = TRUE) +
  labs(title = "Age vs. Baseline GFR", x = "Age", y = "Baseline GFR", color = "Diabetes") +
  theme_minimal()

# ------------------------------
# Step 2: Imputation of Missing Covariates Using Random Forest (missForest)
# ------------------------------

# We will impute only the covariates. For this example, assume we impute: age, baseline_gfr, and bmi.
# (Outcome 'event' and other variables will remain unchanged.)

# Select the covariates to impute
covariates_to_impute <- c("age", "baseline_gfr", "bmi")

# Create a subset for imputation
df_cov <- df[, covariates_to_impute]

# Apply missForest for imputation (using default settings, which uses random forest)
impute_result <- missForest(df_cov, verbose = TRUE)

# Extract the imputed data
df_cov_imputed <- impute_result$ximp

# Replace original covariates in the full dataset with the imputed values
df_imputed <- df
df_imputed[, covariates_to_impute] <- df_cov_imputed

# ------------------------------
# Step 3: Diagnostics After Imputation
# ------------------------------

# Check that missingness in the imputed covariates is removed
sapply(df_imputed[, covariates_to_impute], function(x) sum(is.na(x)))

# Confirm that the outcome variable remains with its original missingness pattern
cat("Missing values in outcome (event):", sum(is.na(df_imputed$event)), "\n")

# Optionally, visualize the distributions after imputation to check consistency:
ggplot(df_imputed, aes(x = age)) +
  geom_histogram(binwidth = 2, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Imputed Age", x = "Age", y = "Frequency") +
  theme_minimal()

# Save the imputed dataset if needed
write.csv(df_imputed, "data/imputed_case_study_data.csv", row.names = FALSE)
