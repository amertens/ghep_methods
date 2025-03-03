# Load necessary libraries
library(dplyr)
library(ggplot2)
library(MatchIt)
library(cobalt)
library(tableone)

# -----------------------------
# Step 0: Load the cleaned and imputed data
# -----------------------------
# Replace the file path with your actual path if needed.
df <- read.csv("data/imputed_case_study_data.csv")

# Assume the dataset has the following columns:
# - Covariates: age, baseline_gfr, bmi, sex, diabetes, hypertension
# - Treatment: treatment (binary: 0/1)
# - Outcome: event (binary outcome for kidney injury)
# - Other variables: time_to_event, observed_time, censored

# -----------------------------
# Step 1: Exploratory Analysis & Pre-Matching Diagnostics
# -----------------------------

# 1.1: Create a summary table of covariates by treatment group
covariates <- c("age", "baseline_gfr", "bmi", "sex", "diabetes", "hypertension")
table1 <- CreateTableOne(vars = covariates, strata = "treatment", data = df, test = FALSE)
print(table1, smd = TRUE)

# 1.2: Estimate the propensity score using logistic regression
# Here we model treatment as a function of key covariates.
ps_model <- glm(treatment ~ age + baseline_gfr + bmi + sex + diabetes + hypertension,
                data = df, family = binomial)
df$pscore <- predict(ps_model, type = "response")

# 1.3: Visualize the distribution of propensity scores by treatment group (Pre-Matching)
ggplot(df, aes(x = pscore, fill = factor(treatment))) +
  geom_density(alpha = 0.5) +
  labs(title = "Propensity Score Distribution (Pre-Matching)",
       x = "Propensity Score", fill = "Treatment") +
  theme_minimal()

# -----------------------------
# Step 2: Perform Propensity Score Matching
# -----------------------------
# We use nearest-neighbor matching with a caliper to ensure good matches.
# Adjust caliper width if needed.
matchit_out <- matchit(treatment ~ age + baseline_gfr + bmi + sex + diabetes + hypertension,
                       data = df, method = "nearest", caliper = 0.2)

# View a summary of the matching
summary(matchit_out)

# -----------------------------
# Step 3: Assess Covariate Balance After Matching
# -----------------------------
# Use cobalt to assess balance
bal_tab <- bal.tab(matchit_out, un = TRUE, m.threshold = 0.1)
print(bal_tab)

# Create a love plot to visually assess balance
love.plot(bal_tab, threshold = 0.1, var.order = "unadjusted", 
          abs = TRUE, line = TRUE, colors = c("red", "blue"),
          title = "Covariate Balance Before and After Matching")

# -----------------------------
# Step 4: Extract the Matched Data and Estimate the Treatment Effect
# -----------------------------
# Extract matched data from the matchit object
matched_df <- match.data(matchit_out)

# Check the size of the matched sample
cat("Matched sample size:", nrow(matched_df), "\n")

# 4.1: Estimate the treatment effect on the outcome (e.g., risk difference)
# Option A: Directly compare proportions in the matched data
risk_treated <- mean(matched_df$event[matched_df$treatment == 1], na.rm = TRUE)
risk_untreated <- mean(matched_df$event[matched_df$treatment == 0], na.rm = TRUE)
risk_difference <- risk_treated - risk_untreated
cat("Risk (Treated):", risk_treated, "\n")
cat("Risk (Untreated):", risk_untreated, "\n")
cat("Risk Difference:", risk_difference, "\n")

# Option B: Fit a logistic regression model on the matched data
outcome_model <- glm(event ~ treatment, data = matched_df, family = binomial)
summary(outcome_model)

# -----------------------------
# Step 5: Visualize Propensity Score Distribution Post-Matching
# -----------------------------
# Plot the distribution of propensity scores in the matched sample
ggplot(matched_df, aes(x = pscore, fill = factor(treatment))) +
  geom_density(alpha = 0.5) +
  labs(title = "Propensity Score Distribution (Post-Matching)",
       x = "Propensity Score", fill = "Treatment") +
  theme_minimal()

# -----------------------------
# Step 6: Additional Diagnostics & Sensitivity Analyses
# -----------------------------
# For sensitivity analyses, you might try:
# - Changing the caliper value (e.g., 0.1 or 0.3) and re-assessing balance.
# - Using alternative matching methods (e.g., optimal matching).
# - Checking robustness of treatment effect estimates with different methods.

# Example: Re-run matching with a different caliper
matchit_out_alt <- matchit(treatment ~ age + baseline_gfr + bmi + sex + diabetes + hypertension,
                           data = df, method = "nearest", caliper = 0.1)
summary(matchit_out_alt)
bal_tab_alt <- bal.tab(matchit_out_alt, un = TRUE, m.threshold = 0.1)
love.plot(bal_tab_alt, threshold = 0.1, var.order = "unadjusted", 
          abs = TRUE, line = TRUE, colors = c("red", "blue"),
          title = "Covariate Balance with Caliper = 0.1")

