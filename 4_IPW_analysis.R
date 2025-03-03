# Load necessary libraries
library(dplyr)
library(ggplot2)
library(ipw)
library(survey)
library(cobalt)
library(tableone)

# -----------------------------
# Step 0: Load the cleaned and imputed dataset
# -----------------------------
df <- read.csv("data/imputed_case_study_data.csv")

# -----------------------------
# Step 1: Estimate the Propensity Score
# -----------------------------
# Model treatment as a function of key covariates using logistic regression
ps_model <- glm(treatment ~ age + baseline_gfr + bmi + sex + diabetes + hypertension,
                data = df, family = binomial)

# Obtain the predicted propensity scores (P(A=1|X))
df$pscore <- predict(ps_model, type = "response")

# -----------------------------
# Step 2: Calculate Stabilized Weights
# -----------------------------
# Compute the marginal probability of treatment assignment
pA1 <- mean(df$treatment == 1)
pA0 <- 1 - pA1

# Calculate stabilized weights:
# For treated: weight = P(A=1) / pscore
# For controls: weight = P(A=0) / (1 - pscore)
df <- df %>%
  mutate(stab_weight = ifelse(treatment == 1, pA1 / pscore, pA0 / (1 - pscore)))

# Examine summary of weights
summary(df$stab_weight)

# -----------------------------
# Step 3: Diagnostics for Weighting
# -----------------------------
# 3.1: Visualize the distribution of stabilized weights
ggplot(df, aes(x = stab_weight)) +
  geom_histogram(binwidth = 0.5, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Stabilized Weights", x = "Stabilized Weight", y = "Frequency") +
  theme_minimal()

# 3.2: Assess covariate balance after weighting
# Create a survey design object using the stabilized weights
design <- svydesign(ids = ~1, data = df, weights = ~stab_weight)

# Compare covariate balance before weighting using tableone:
covariates <- c("age", "baseline_gfr", "bmi", "sex", "diabetes", "hypertension")
table1_before <- CreateTableOne(vars = covariates, strata = "treatment", data = df, test = FALSE)
print(table1_before, smd = TRUE)

# Assess balance after weighting using cobalt
bal_tab <- bal.tab(treatment ~ age + baseline_gfr + bmi + sex + diabetes + hypertension, 
                   data = df, weights = df$stab_weight, estimand = "ATE")
print(bal_tab)

# Visualize balance with a love plot
love.plot(bal_tab, threshold = 0.1, var.order = "unadjusted", 
          abs = TRUE, line = TRUE, colors = c("red", "blue"),
          title = "Covariate Balance Before and After Weighting")

# -----------------------------
# Step 4: Fit the Weighted Outcome Model
# -----------------------------
# Use the survey design to fit a weighted logistic regression model for the outcome
weighted_model <- svyglm(event ~ treatment, design = design, family = binomial)
summary(weighted_model)

# Optional: Calculate predicted probabilities to obtain risk estimates.
# For example, compute marginal risk for treated and untreated:
predicted <- predict(weighted_model, type = "response")
# (Alternatively, one can use the design-based approach to estimate group means.)

# -----------------------------
# Step 5: Sensitivity Analyses (Optional)
# -----------------------------
# For example, truncate extreme weights and re-run the analysis:
df <- df %>%
  mutate(trunc_weight = ifelse(stab_weight > 20, 20, stab_weight))
design_trunc <- svydesign(ids = ~1, data = df, weights = ~trunc_weight)
weighted_model_trunc <- svyglm(event ~ treatment, design = design_trunc, family = binomial)
summary(weighted_model_trunc)

# -----------------------------
# End of IPW Analysis Script
# -----------------------------
