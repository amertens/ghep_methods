
# Load necessary libraries
library(MASS)
library(tidyverse)

# Set seed for reproducibility
set.seed(42)

# Define sample size
n <- 5000

# Simulate baseline covariates
age <- rnorm(n, mean = 60, sd = 10)  # Age centered at 60
sex <- rbinom(n, 1, 0.5)  # Binary sex variable (0 or 1)
diabetes <- rbinom(n, 1, 0.3)  # 30% prevalence of diabetes
hypertension <- rbinom(n, 1, 0.4)  # 40% prevalence of hypertension
baseline_gfr <- rnorm(n, mean = 90, sd = 15)  # Baseline kidney function (eGFR)
bmi <- rnorm(n, mean = 28, sd = 5)  # BMI distribution

# Introduce non-linearity in the relationship between age and treatment
age_effect <- exp(-0.05 * (age - 60)^2)  # Exponential effect centered at 60

# Interaction terms affecting treatment assignment
interaction_term <- (diabetes * hypertension) + 0.5 * (sex * bmi)

# Simulate treatment assignment (propensity depends on confounders non-linearly)
treatment_prob <- plogis(-1 + 0.02 * age + 0.4 * diabetes + 0.3 * hypertension +
                           0.2 * sex + 0.01 * baseline_gfr + 0.02 * bmi +
                           interaction_term - age_effect)

treatment <- rbinom(n, 1, treatment_prob)

# Introduce non-linearity in the relationship between age and outcome
age_outcome_effect <- sin(age / 10)

# Simulate time to kidney injury outcome (survival outcome)
baseline_hazard <- 0.02  # Baseline hazard rate
outcome_prob <- plogis(-2 + 0.03 * age + 0.5 * diabetes + 0.4 * hypertension +
                         0.3 * sex + 0.015 * baseline_gfr + 0.025 * bmi +
                         interaction_term - 0.5 * treatment - age_outcome_effect)

event <- rbinom(n, 1, outcome_prob)  # Binary outcome for kidney injury

# Simulate time-to-event (higher probabilities lead to shorter time-to-event)
time_to_event <- rexp(n, rate = baseline_hazard * (1 + outcome_prob))

# Censor some observations randomly
censoring_time <- rexp(n, rate = 1/10)  # Administrative censoring at random times
observed_time <- pmin(time_to_event, censoring_time)
censored <- as.integer(time_to_event > censoring_time)

# Create a data frame
df <- data.frame(
  age, sex, diabetes, hypertension, baseline_gfr, bmi,
  treatment, event, time_to_event, observed_time, censored
)


# Load necessary packages
library(dplyr)

# Set seed for reproducibility of missingness
set.seed(123)

# Assume 'df' is your simulated data frame with the following columns:
# age, sex, diabetes, hypertension, baseline_gfr, bmi, treatment, event, time_to_event, observed_time, censored

# --------------------------
# Introduce missingness in covariates:
# --------------------------

# 1. Age: Missingness probability depends on diabetes status and deviation from age 60.
df <- df %>%
  mutate(
    p_age_missing = plogis(-2 + 0.05 * diabetes + 0.01 * (age - 60)),
    age = ifelse(runif(n()) < p_age_missing, NA, age)
  )

# 2. Baseline kidney function (baseline_gfr): Missingness depends on hypertension and the current gfr value.
df <- df %>%
  mutate(
    p_gfr_missing = plogis(-1 + 0.03 * hypertension - 0.01 * baseline_gfr),
    baseline_gfr = ifelse(runif(n()) < p_gfr_missing, NA, baseline_gfr)
  )

# 3. BMI: Missingness influenced by sex and age.
df <- df %>%
  mutate(
    p_bmi_missing = plogis(-1.5 + 0.02 * sex + 0.01 * (age - 60)),
    bmi = ifelse(runif(n()) < p_bmi_missing, NA, bmi)
  )

# --------------------------
# Introduce missingness in the outcome:
# --------------------------

# 4. Outcome (event): Now, missingness probability depends on age, diabetes, and treatment status.
# Here, we assume that being treated (treatment == 1) reduces the probability of missing outcome data.
df <- df %>%
  mutate(
    p_event_missing = plogis(-1 + 0.03 * age - 0.5 * diabetes - 0.5 * treatment),
    event = ifelse(runif(n()) < p_event_missing, NA, event)
  )

# Optionally, remove the temporary probability columns if not needed
df <- df %>% subset(., select=-c(p_age_missing, p_gfr_missing, p_bmi_missing, p_event_missing))

# Check missingness summary for key variables
summary(df[, c("age", "baseline_gfr", "bmi", "event")])


# Save to CSV
write.csv(df, "data/simulated_case_study_data.csv", row.names = FALSE)

# View first few rows
head(df)


#-------------------------------------------------------------------------------
# Calculate truth
#-------------------------------------------------------------------------------

# Define the function for probability of kidney injury
calc_risk <- function(age, diabetes, hypertension, sex, baseline_gfr, bmi, interaction_term, treatment) {
  lin_pred <- -2 + 0.03 * age + 0.5 * diabetes + 0.4 * hypertension +
    0.3 * sex + 0.015 * baseline_gfr + 0.025 * bmi + interaction_term - 0.5 * treatment
  return(1 / (1 + exp(-lin_pred)))
}

# Compute expected risks under treatment and no treatment
df <- df %>%
  mutate(
    risk_treated = calc_risk(age, diabetes, hypertension, sex, baseline_gfr, bmi, (diabetes * hypertension + 0.5 * sex * bmi), 1),
    risk_untreated = calc_risk(age, diabetes, hypertension, sex, baseline_gfr, bmi, (diabetes * hypertension + 0.5 * sex * bmi), 0)
  )

# Compute true causal estimand
true_risk_treated <- mean(df$risk_treated)
true_risk_untreated <- mean(df$risk_untreated)

true_risk_difference <- true_risk_treated - true_risk_untreated
true_risk_ratio <- true_risk_treated / true_risk_untreated

# Print results
cat("True Risk (Treated):", true_risk_treated, "\n")
cat("True Risk (Untreated):", true_risk_untreated, "\n")
cat("True Risk Difference (RD):", true_risk_difference, "\n")
cat("True Risk Ratio (RR):", true_risk_ratio, "\n")
