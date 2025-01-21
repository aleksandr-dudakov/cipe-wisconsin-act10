# ===============================================================
# Causal Inference and Policy Evaluation
# Project: Effects of Wisconsin Act 10 on Unions, Wages, and Benefits
# Author: Aleksandr Dudakov
# ===============================================================
# -------------------------
# Load Required Libraries
# -------------------------
# To install the package for loading CPS ORG data from the EPI microdata service, 
# uncomment the line below and run it:
#install.packages("epiextractr", repos = c("https://economic.r-universe.dev", "https://cloud.r-project.org"))
library(tidyverse)    # Data manipulation and visualization
library(epiextractr)  # CPS data import
library(estimatr)     # Robust regression analysis
library(VIM)          # Hot deck imputation
library(fixest)       # Fixed effects models
library(haven)        # Reading .dta files
library(HonestDiD)    # Honest Difference-in-Differences
library(Synth)        # Synthetic Control Method
library(SCtools)      # Tools for Synthetic Control
library(SimDesign)    # Quiet output

# -------------------------
# Set Working Directory
# -------------------------
# setwd("/Users/adudakov/Documents/data science for economics/causal inference and policy evaluation/CIPE_project_Dudakov")

# ===============================================================
# 1. Unions and wages
# ===============================================================
# ---------------------------------------------------------------
# 1.1 Data Loading and Preprocessing
# ---------------------------------------------------------------
# Reference for CPS ORG EPI microdata service: https://microdata.epi.org/
# Load CPS ORG data
cps_data <- load_org(
  1994:2023,
  year, month, age, female, citizen, wbhao, married, metstat,
  statecensus, gradeatn, cow1, emp, mind16, mocc10, 
  unmem, union, wage, weekpay, wageotc, orgwgt,
  .extracts_dir = "data" # here write the name of the folder with data files
) %>%
  as_factor() %>%
  # Handle age: convert "80+" to 80 and others to integer
  mutate(
    age = as.integer(ifelse(age == "80+", 80, as.character(age)))
  ) %>%
  # Filter for employed, wage/salary workers of working age
  filter(
    emp == "Employed",
    cow1 != "Without pay",
    age >= 18 & age <= 64,
    # union friendly states and WI
    statecensus %in% c("SD", "NM", "MO", "FL", "KS", "NE", "DE", "MT", 
                       "OH", "MI", "IL", "VT", "NH", "MA", "OR", "HI",
                       "MN", "WA", "ME", "CA", "RI", "PA", "NJ", "NY", "WI")
  ) %>%
  # Rename variables for clarity
  rename(
    occupation = mocc10,
    sex = female,
    race = wbhao,
    metropolitan = metstat,
    state = statecensus,
    education_level = gradeatn,
    class = cow1,
    industry = mind16,
    weight = orgwgt
  ) %>%
  # Create treatment, post, and time-to-treatment indicators; recode union, union_member, and weekpay
  mutate(
    treatment = as.integer(state == "WI"),
    post = case_when(
      #year == 2011 & month > 6 ~ 1,
      year > 2011               ~ 1,
      TRUE                      ~ 0
    ),
    time_to_treatment = factor(year - 2011, levels = sort(unique(year - 2011))),
    union = as.integer(union == "Union represented"),
    union_member = as.integer(unmem == "Union member"),
    public = as.integer(class %in% c("Government - State", "Government - Local") & 
                          occupation != "Protective service"),
    weekpay = case_when(
      weekpay == 0 & is.na(wage) ~ NA, # Set weekpay to NA if it's zero and hourly wage is not null
      weekpay == 0               ~ 0.1, # Set weekpay to 0.1 for log transformation to work
      TRUE                       ~ weekpay
    ),
    year_factor = as.character(factor(year, levels = sort(unique(year))))
  ) %>%
  # Select only the relevant variables
  select(
    year, year_factor, month, age, sex, citizen, race, married, metropolitan, state, education_level, industry,
    wage, wageotc, weekpay, union, union_member, public, weight, treatment, post, time_to_treatment
  )

# Quick structure check
str(cps_data)

# ---------------------------------------------------------------
# 1.2 Missing Value Imputation
# ---------------------------------------------------------------
# Summarize missing values before imputation
missing_summary_before <- cps_data %>%
  summarise(across(everything(), ~ mean(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "missing_percentage")

cat("Missing Data Summary Before Imputation:\n")
print(missing_summary_before, n = Inf)

# Perform hot deck imputation
set.seed(123)
cps_data <- cps_data %>%
  hotdeck(
    variable = c("metropolitan", "wage", "wageotc", "weekpay", "union", "union_member"),
    domain_var = c("year", "state", "public")
  ) %>%
  as_tibble() %>%
  mutate(
    is_imputed = rowSums(select(., ends_with("_imp"))) > 0
  ) %>%
  select(-ends_with("_imp"))

# Summarize missing values after imputation (before dropping NAs)
missing_summary_after <- cps_data %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "missing_count")

cat("\nMissing Data Summary After Imputation (Count of NAs):\n")
print(missing_summary_after, n = Inf)
# No missing values

# Drop rows with any remaining NA values
cps_data <- cps_data %>% drop_na() 

# At this point, 'cps_data' is cleaned and imputed, ready for analysis.
# ---------------------------------------------------------------
# 1.3 Visualizing Raw Data Trends
# ---------------------------------------------------------------
# Precompute log-transformed variables for analysis
cps_data <- cps_data %>%
  mutate(
    log_wage = log(wage),
    log_wageotc = log(wageotc),
    log_weekpay = log(weekpay)
  )

# Define the outcome variables to analyze with their labels
variables_to_analyze <- list(
  "Union Membership Rate" = "union_member",
  "Union Representation Rate" = "union",
  "Log of Hourly Earnings" = "log_wage",
  "Log of Hourly Earnings (OTC)" = "log_wageotc",
  "Log of Weekly Earnings (OTC)" = "log_weekpay"
)

academic_theme <- theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

# Function to plot raw trends for a given outcome with better colors for academic papers
plot_raw_trends <- function(data, outcome_var, outcome_label) {
  raw_summary <- data %>%
    group_by(year, treatment) %>%
    summarise(
      avg_value = mean(!!sym(outcome_var), na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(raw_summary, aes(x = year, y = avg_value, color = factor(treatment))) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    geom_vline(xintercept = 2011, linetype = "dotted", color = "darkgray", linewidth = 0.8) +
    annotate(
      "text", x = 2011.5, y = max(raw_summary$avg_value, na.rm = TRUE) * 0.95,
      label = "Year of Policy Implementation", hjust = 0, vjust = 0.5, color = "darkgray", size = 4
    ) +
    labs(
      title = paste("Temporal Trends in", outcome_label),
      x = "Year",
      y = paste("Average", outcome_label),
      color = "Group"
    ) +
    scale_color_manual(
      name = "Group",
      values = c("0" = "#4E79A7", "1" = "#F28E2B"), # Muted navy blue for control, soft orange for treatment
      labels = c("Control Group", "WI")
    ) +
    scale_x_continuous(
      breaks = seq(min(raw_summary$year) + 2, max(raw_summary$year), by = 3) # Year labels every 5 years
    ) +
    academic_theme
}

# Plot raw trends for each outcome variable
for (var_name in names(variables_to_analyze)) {
  var <- variables_to_analyze[[var_name]]
  print(plot_raw_trends(
    cps_data %>% filter(public == 1, state %in% c("WI", "IL", "KS", "MI", "MN", "SD", "MA", "PA", "WA")), 
    var, var_name))
}

print(plot_raw_trends(
  cps_data %>% filter(state %in% c("WI", "IL", "KS", "MI", "MN", "SD", "MA", "PA", "WA")), 
  "public", "Share of Public Sector Employees"))

# ---------------------------------------------------------------
# 1.4 Event Study Analysis for DD (with HonestDiD)
# ---------------------------------------------------------------
# General function to run event study and perform HonestDiD analysis
run_event_study_honestdid <- function(data, outcome_var, outcome_label, numPrePeriods, numPostPeriods, include_honest = TRUE) {
  # Construct the formula for the event study
  formula_str <- paste0(
    outcome_var, " ~ i(year_factor, treatment, ref='2010') + year_factor + treatment +",
    "age + sex + married + race + metropolitan + industry + education_level | state"
  )
  
  # Fit the event study model using fixed effects and cluster by state
  event_study_model <- feols(
    as.formula(formula_str),
    data = data,
    weights = ~weight,
    cluster = ~state
  )
  
  cat("\n\nEvent Study Results for:", outcome_label, "\n")
  
  # Plot the event study estimates
  iplot(
    event_study_model,
    main = paste0("Event Study: Impact of Treatment on ", outcome_label),
    xlab = "Year",
    ylab = paste0("Effect on ", outcome_label),
    col.point = "blue",
    col.line = "blue",
    x.cross = TRUE,
    zero.at = 0,
    grid = TRUE,
    order = as.character(sort(as.numeric(unique(cps_data$year_factor))))
  )
  
  if (include_honest) {
    
    # Extract coefficients and covariance matrix for event study periods
    all_coefficients <- coef(event_study_model)
    all_covariance <- vcov(event_study_model)
    
    # Identify event study coefficients (those prefixed with time_to_treatment)
    event_study_indices <- grep("^year_factor::", names(all_coefficients))
    
    # Extract betahat and sigma for HonestDiD
    betahat <- all_coefficients[event_study_indices]
    sigma <- all_covariance[event_study_indices, event_study_indices, drop = FALSE]
    
    # Run HonestDiD sensitivity analysis
    delta_rm_results <- HonestDiD::createSensitivityResults_relativeMagnitudes(
      betahat = betahat,
      sigma = sigma,
      numPrePeriods = numPrePeriods,
      numPostPeriods = numPostPeriods,
      Mbarvec = seq(0.5, 2, by = 0.5), 
      parallel = TRUE
    )
    
    cat("\nHonestDiD Sensitivity Results (Relative Magnitudes) for:", outcome_label, "\n")
    print(delta_rm_results)
    
    # Construct original confidence set results
    originalResults <- HonestDiD::constructOriginalCS(
      betahat = betahat,
      sigma = sigma,
      numPrePeriods = numPrePeriods,
      numPostPeriods = numPostPeriods
    )
    
    # Plot sensitivity results
    HonestDiD::createSensitivityPlot_relativeMagnitudes(delta_rm_results, originalResults) %>% plot()
  }
  # Return the model invisibly
  invisible(event_study_model)
}

# Run the event study + HonestDiD for all main outcomes
for (var_name in names(variables_to_analyze)) {
  var <- variables_to_analyze[[var_name]]
  run_event_study_honestdid(
    data = cps_data %>%
      filter(public == 1, state %in% c("WI", "IL", "KS", "MI", "MN", "SD", "MA", "PA", "WA")),
    outcome_var = var,
    outcome_label = var_name,
    numPrePeriods = 17,
    numPostPeriods = 12,
    include_honest = T
  )
}

# Run separately for share of public employees
run_event_study_honestdid(
  data = cps_data %>% 
    filter(state %in% c("WI", "IL", "KS", "MI", "MN", "SD", "MA", "PA", "WA")),
  outcome_var = "public",
  outcome_label = "Public Worker Status",
  numPrePeriods = 17,
  numPostPeriods = 12
)

# ---------------------------------------------------------------
# 1.5 Difference-in-Differences (DiD) Analysis
# ---------------------------------------------------------------
# Specifications for the main analysis
did_specifications <- list(
  "Without Imputed Data and Controls" = list(
    filter = function(data) data %>% filter(state %in% c("WI", "IL", "KS", "MI", "MN", "SD", "MA", "PA", "WA"), is_imputed == 0),
    controls_formula = "",
    fixed_effects = NULL
  ),
  "With Imputed Data and Without Controls" = list(
    filter = function(data) data %>% filter(state %in% c("WI", "IL", "KS", "MI", "MN", "SD", "MA", "PA", "WA")),
    controls_formula = "",
    fixed_effects = NULL
  ),
  "With Imputed Data and With Controls" = list(
    filter = function(data) data %>% filter(state %in% c("WI", "IL", "KS", "MI", "MN", "SD", "MA", "PA", "WA")),
    controls_formula = "age + sex + married + race + metropolitan + education_level",
    fixed_effects = NULL
  ),
  "With Imputed Data and With Controls and State FE" = list(
    filter = function(data) data %>% filter(state %in% c("WI", "IL", "KS", "MI", "MN", "SD", "MA", "PA", "WA")),
    controls_formula = "age + sex + married + race + metropolitan + education_level",
    fixed_effects = "state"
  ),
  "With Imputed Data and With Controls and State and Year FE" = list(
    filter = function(data) data %>% filter(state %in% c("WI", "IL", "KS", "MI", "MN", "SD", "MA", "PA", "WA")),
    controls_formula = "age + sex + married + race + metropolitan + education_level + industry",
    fixed_effects = "state+time_to_treatment"
  ),
  "With Imputed Data and With Controls and State and Year FE and MA-PA-WA" = list(
    filter = function(data) data %>% filter(state %in% c("WI", "MA", "PA", "WA")),
    controls_formula = "age + sex + married + race + metropolitan + education_level + industry",
    fixed_effects = "state+time_to_treatment"
  ),
  "With Imputed Data and With Controls and State and Year FE and Year <= 2019" = list(
    filter = function(data) data %>% filter(state %in% c("WI", "IL", "KS", "MI", "MN", "SD", "MA", "PA", "WA"), year <= 2019),
    controls_formula = "age + sex + married + race + metropolitan + education_level + industry",
    fixed_effects = "state+time_to_treatment"
  ),
  "With Imputed Data and With Controls and State and Year FE and Year <= 2015" = list(
    filter = function(data) data %>% filter(state %in% c("WI", "IL", "KS", "MI", "MN", "SD", "MA", "PA", "WA"), year <= 2015),
    controls_formula = "age + sex + married + race + metropolitan + education_level + industry",
    fixed_effects = "state+time_to_treatment"
  )
)

run_did_model_and_print <- function(data, dep_var, spec_filter) {
  model_data <- spec_filter$filter(data)
  controls_formula <- spec_filter$controls_formula
  fixed_effects_str <- spec_filter$fixed_effects
  
  # Construct formula for the DiD model
  formula_str <- paste0(dep_var, " ~ treatment * post", if (controls_formula != "") paste0(" + ", controls_formula) else "")
  
  # Run the DiD model based on the fixed_effects_str
  if (is.null(fixed_effects_str)) {
    # No fixed effects
    model <- lm_robust(
      as.formula(formula_str),
      data = model_data,
      weights = weight,
      clusters = state,
      se_type = "stata"
    )
  } else if (fixed_effects_str == "state") {
    # Fixed effects for state
    model <- lm_robust(
      as.formula(formula_str),
      data = model_data,
      weights = weight,
      clusters = state,
      fixed_effects = ~state,
      se_type = "stata"
    )
  } else if (fixed_effects_str == "time_to_treatment") {
    # Fixed effects for time_to_treatment
    model <- lm_robust(
      as.formula(formula_str),
      data = model_data,
      weights = weight,
      clusters = state,
      fixed_effects = ~time_to_treatment,
      se_type = "stata"
    )
  } else if (fixed_effects_str == "state+time_to_treatment") {
    # Fixed effects for state and time_to_treatment
    model <- lm_robust(
      as.formula(formula_str),
      data = model_data,
      weights = weight,
      clusters = state,
      fixed_effects = ~state + time_to_treatment,
      se_type = "stata"
    )
  } else {
    stop("Invalid fixed_effects specification provided in spec_filter.")
  }
  
  # Extract and print the treatment:post coefficient
  filtered_term <- broom::tidy(model) %>% filter(term == "treatment:post")
  print(filtered_term)
  
  invisible(model) # Return the model invisibly if needed
}

# Loop through variables and specifications (excluding `public`)
did_results <- lapply(names(variables_to_analyze), function(var_name) {
  dep_var <- variables_to_analyze[[var_name]]
  
  spec_results <- lapply(names(did_specifications), function(spec_name) {
    cat("\n\nDependent Variable:", var_name, "\n")
    cat("Specification:", spec_name, "\n")
    
    # Run the DiD model for this variable and specification
    model <- run_did_model_and_print(
      data = cps_data |> filter(public == 1), 
      dep_var = dep_var, 
      spec_filter = did_specifications[[spec_name]]
    )
    list(spec_name = spec_name, model = model)
  })
  
  list(dep_var = var_name, results = spec_results)
})

# Separate DiD Analysis for `public`
cat("\n\nDependent Variable: Public Worker Status\n")
did_results_public <- lapply(names(did_specifications), function(spec_name) {
  cat("\nSpecification:", spec_name, "\n")
  
  # Run the DiD model for `public`
  model <- run_did_model_and_print(
    data = cps_data,
    dep_var = "public",
    spec_filter = did_specifications[[spec_name]]
  )
  list(spec_name = spec_name, model = model)
})

# ---------------------------------------------------------------
# 1.6 DDD Event Study
# ---------------------------------------------------------------
# A helper function to run and plot a DDD event study. 
run_ddd_event_study <- function(data, outcome_var, outcome_label) {
  # Construct formula for DDD event study model
  # This model includes interaction terms with treatment, public, and time_to_treatment
  formula_str <- paste0(
    outcome_var, " ~ i(year_factor, treatment*public, ref='2010') + ",
    "year_factor*public + year_factor*treatment + ",
    "age + sex + education_level + married + race + metropolitan + industry | state"
  )
  
  # Fit the DDD event study model using fixest
  ddd_event_model <- feols(
    as.formula(formula_str),
    data = data,
    weights = ~weight,
    cluster = ~state
  )
  
  cat("\n\nDDD Event Study for:", outcome_label, "\n")
  print(summary(ddd_event_model))
  
  # Plot the DDD event study results using iplot from fixest
  iplot(
    ddd_event_model,
    main = paste0("DDD Event Study: Impact on ", outcome_label),
    xlab = "Year",
    ylab = paste0("Effect on ", outcome_label),
    col.point = "blue",
    col.line = "blue",
    x.cross = TRUE,
    zero.at = 0,
    grid = TRUE,
    order = as.character(sort(as.numeric(unique(cps_data$year_factor))))
  )
  
  invisible(ddd_event_model)
}


for (var_name in names(variables_to_analyze)) {
  var <- variables_to_analyze[[var_name]]
  run_ddd_event_study(data = cps_data %>%
                        filter(state %in% c("WI", "IL", "KS", "MI", "MN", "SD", "MA", "PA", "WA")), 
                      outcome_var = var, outcome_label = var_name)
}

# ---------------------------------------------------------------
# 1.7 DDD Regression
# ---------------------------------------------------------------
# Create a helper function to run the DDD model and print the treatment:post:public coefficient
run_ddd_model_and_print <- function(data, dep_var, spec_filter) {
  model_data <- spec_filter$filter(data)
  controls_formula <- spec_filter$controls_formula
  fixed_effects_str <- spec_filter$fixed_effects
  
  # Construct formula for the DDD model
  formula_str <- paste0(dep_var, " ~ treatment * post * public", 
                        if (controls_formula != "") paste0(" + ", controls_formula) else "")
  
  # Run the DDD model based on the fixed_effects_str
  if (is.null(fixed_effects_str)) {
    # No fixed effects
    model <- lm_robust(
      as.formula(formula_str),
      data = model_data,
      weights = weight,
      clusters = state,
      se_type = "stata"
    )
  } else if (fixed_effects_str == "state") {
    # Fixed effects for state
    model <- lm_robust(
      as.formula(formula_str),
      data = model_data,
      weights = weight,
      clusters = state,
      fixed_effects = ~state,
      se_type = "stata"
    )
  } else if (fixed_effects_str == "state+time_to_treatment") {
    # Fixed effects for state and time_to_treatment
    model <- lm_robust(
      as.formula(formula_str),
      data = model_data,
      weights = weight,
      clusters = state,
      fixed_effects = ~state + time_to_treatment,
      se_type = "stata"
    )
  } else {
    stop("Invalid fixed_effects specification provided in spec_filter.")
  }
  
  # Extract and print the treatment:post:public coefficient
  filtered_term <- broom::tidy(model) %>% filter(term == "treatment:post:public")
  print(filtered_term)
  
  invisible(model) # Return the model invisibly 
}

# Run the DDD regressions for each variable and specification
ddd_results <- lapply(names(variables_to_analyze), function(var_name) {
  dep_var <- variables_to_analyze[[var_name]]
  
  spec_results <- lapply(names(did_specifications), function(spec_name) {
    cat("\n\nDependent Variable:", var_name, "\n")
    cat("Specification:", spec_name, "\n")
    
    # Run the DDD model for this variable and specification
    model <- run_ddd_model_and_print(
      data = cps_data %>% filter(state %in% c("WI", "IL", "KS", "MI", "MN", "SD", "MA", "PA", "WA")),
      dep_var = dep_var,
      spec_filter = did_specifications[[spec_name]]
    )
    list(spec_name = spec_name, model = model)
  })
  
  list(dep_var = var_name, results = spec_results)
})

# ===============================================================
# 2. Benefits
# ===============================================================
# ---------------------------------------------------------------
# 2.1 Load and Preprocess Benefits Data
# ---------------------------------------------------------------
benefits_data <- read_dta("data/benefits_data.dta") %>%
  filter(asecflag == 1) %>%
  as_factor() %>%
  mutate(
    # Convert statecensus to two-letter abbreviations
    state = recode(
      statecensus,
      "maine" = "ME", "new hampshire" = "NH", "vermont" = "VT", "massachusetts" = "MA",
      "rhode island" = "RI", "connecticut" = "CT", "new york" = "NY", "new jersey" = "NJ",
      "pennsylvania" = "PA", "ohio" = "OH", "indiana" = "IN", "illinois" = "IL",
      "michigan" = "MI", "wisconsin" = "WI", "minnesota" = "MN", "iowa" = "IA",
      "missouri" = "MO", "north dakota" = "ND", "south dakota" = "SD", "nebraska" = "NE",
      "kansas" = "KS", "delaware" = "DE", "maryland" = "MD", "district of columbia" = "DC",
      "virginia" = "VA", "west virginia" = "WV", "north carolina" = "NC", 
      "south carolina" = "SC", "georgia" = "GA", "florida" = "FL", "kentucky" = "KY",
      "tennessee" = "TN", "alabama" = "AL", "mississippi" = "MS", "arkansas" = "AR",
      "louisiana" = "LA", "oklahoma" = "OK", "texas" = "TX", "montana" = "MT",
      "idaho" = "ID", "wyoming" = "WY", "colorado" = "CO", "new mexico" = "NM",
      "arizona" = "AZ", "utah" = "UT", "nevada" = "NV", "washington" = "WA",
      "oregon" = "OR", "california" = "CA", "alaska" = "AK", "hawaii" = "HI"
    ),
    
    age = as.integer(as.character(age)),
    married = ifelse(marst == "married, spouse present", 1, 0),
    citizen = ifelse(citizen %in% c("born in u.s", "born in u.s. outlying", 
                                    "born abroad of american parents", "naturalized citizen"), 1, 0),
    health_insurance_included = ifelse(inclugh == "yes", 1, 0),
    health_paid_by_employer = ifelse(paidgh %in% c("yes, paid for part", "yes, paid for all"), 1, 0),
    health_paid_by_employer_full = ifelse(paidgh == "yes, paid for all", 1, 0),
    pension_plan_at_work = case_when(
      pension %in% c("included in pension plan at work", "pension plan at work, but not included") ~ 1,
      pension == "no pension plan at work" ~ 0,
      TRUE ~ 0
    ),
    race = case_when(
      race == "white" ~ "white",
      race == "black" ~ "black",
      race == "asian only" ~ "asian",
      race %in% c("american indian/aleut/eskimo", "white-american indian") ~ "american indian",
      TRUE ~ "other"
    ),
    metropolitan = ifelse(metro == "in central/principal city", 1, 0),
    treatment = ifelse(state == "WI", 1, 0),
    post = ifelse(year > 2011, 1, 0),
    public = ifelse(
      classwkr %in% c("state government employee", "local government employee") &
        ind1990 != "justice, public order, and safety", 1, 0
    ),
    time_to_treatment = factor(year - 2011),
    work_expenses = ifelse(wkxpns == 0, 0.1, wkxpns),
    employer_insurance_contribution = ifelse(emcontrb == 0, 0.1, emcontrb),
    year_factor = as.character(factor(year, levels = sort(unique(year)))),
    is_imputed = 0
  ) %>%
  filter(
    year >= 1994 & year <= 2018,
    age >= 18 & age < 65,
    labforce == "yes, in the labor force",
    classwkr != "niu"
  ) %>%
  rename(
    education_level = educ,
    industry = ind1990,
    weight = asecwt
  ) %>%
  select(
    year, year_factor, month, age, sex, citizen, race, married, metropolitan, state, education_level, industry,
    pension_plan_at_work, employer_insurance_contribution, health_insurance_included, 
    health_paid_by_employer, health_paid_by_employer_full,
    treatment, post, time_to_treatment, public, weight, is_imputed
  )

# Check structure
str(benefits_data)

# Summarize missing values
missing_summary_benefits <- benefits_data %>%
  summarise(across(everything(), ~ mean(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "missing_percentage")
print(missing_summary_benefits, n = Inf)

# ---------------------------------------------------------------
# 2.2 Visualizing Raw Data Trends for Benefits
# ---------------------------------------------------------------
# Precompute log-transformed variables for analysis
benefits_data <- benefits_data %>%
  mutate(
    log_employer_insurance_contribution = log(employer_insurance_contribution)
  )

# Define the outcome variables to analyze with their labels
benefits_variables_to_analyze <- list(
  "Participation in Workplace Pension" = "pension_plan_at_work",
  "Log of Employer Health Insurance Contribution" = "log_employer_insurance_contribution",
  "Group Health Insurance Coverage" = "health_insurance_included"
)

# Plot raw trends for each outcome variable in benefits_data using the previously defined plot_raw_trends function
for (var_name in names(benefits_variables_to_analyze)) {
  var <- benefits_variables_to_analyze[[var_name]]
  print(plot_raw_trends(
    benefits_data %>% filter(public == 1, state %in% c("WI", "IL", "KS", "MI", "MN", "SD", "MA", "PA", "WA")), 
    var, var_name))
}

# ---------------------------------------------------------------
# 2.3 Event Study Analysis for DD (Benefits)
# ---------------------------------------------------------------
# Using the previously defined run_event_study_honestdid function
for (var_name in names(benefits_variables_to_analyze)) {
  var <- benefits_variables_to_analyze[[var_name]]
  run_event_study_honestdid(
    data = benefits_data %>%
      filter(public == 1, state %in% c("WI", "IL", "KS", "MI", "MN", "SD", "MA", "PA", "WA")),
    outcome_var = var,
    outcome_label = var_name,
    numPrePeriods = 17,
    numPostPeriods = 7
  )
}

# ---------------------------------------------------------------
# 2.4 Difference-in-Differences (DiD) Analysis (Benefits)
# ---------------------------------------------------------------
benefits_did_results <- lapply(names(benefits_variables_to_analyze), function(var_name) {
  dep_var <- benefits_variables_to_analyze[[var_name]]
  
  spec_results <- lapply(names(did_specifications), function(spec_name) {
    cat("\n\nDependent Variable:", var_name, "\n")
    cat("Specification:", spec_name, "\n")
    
    # Run the DiD model for this variable and specification
    model <- run_did_model_and_print(
      data = benefits_data |> filter(public == 1), 
      dep_var = dep_var, 
      spec_filter = did_specifications[[spec_name]]
    )
    list(spec_name = spec_name, model = model)
  })
  
  list(dep_var = var_name, results = spec_results)
})

# ---------------------------------------------------------------
# 2.5 DDD Event Study (Benefits)
# ---------------------------------------------------------------
for (var_name in names(benefits_variables_to_analyze)) {
  var <- benefits_variables_to_analyze[[var_name]]
  run_ddd_event_study(
    data = benefits_data %>%
      filter(state %in% c("WI", "IL", "KS", "MI", "MN", "SD", "MA", "PA", "WA")),
    outcome_var = var,
    outcome_label = var_name
  )
} 

for (var_name in names(benefits_variables_to_analyze)) {
  var <- benefits_variables_to_analyze[[var_name]]
  run_ddd_event_study(
    data = benefits_data %>%
      filter(state %in% c("WI", "MA", "PA", "WA")),
    outcome_var = var,
    outcome_label = var_name
  )
} 

# ---------------------------------------------------------------
# 2.6 DDD Regression (Benefits)
# ---------------------------------------------------------------
ddd_results <- lapply(names(benefits_variables_to_analyze), function(var_name) {
  dep_var <- benefits_variables_to_analyze[[var_name]]
  
  spec_results <- lapply(names(did_specifications), function(spec_name) {
    cat("\n\nDependent Variable:", var_name, "\n")
    cat("Specification:", spec_name, "\n")
    
    # Run the DDD model for this variable and specification
    model <- run_ddd_model_and_print(
      data = benefits_data %>% filter(state %in% c("WI", "IL", "KS", "MI", "MN", "SD", "MA", "PA", "WA")),
      dep_var = dep_var,
      spec_filter = did_specifications[[spec_name]]
    )
    list(spec_name = spec_name, model = model)
  })
  
  list(dep_var = var_name, results = spec_results)
})

# ===============================================================
# 3. Synthetic Control Analysis
# ===============================================================
# ---------------------------------------------------------------
# 3.1 Data Preparation for Synthetic Control (CPS)
# ---------------------------------------------------------------
# Aggregate CPS data at the state-year level for public employees
scm_data <- cps_data %>%
  filter(public == 1) %>%
  group_by(state, year) %>%
  summarise(
    avg_union_member = weighted.mean(union_member, weight),
    avg_union = weighted.mean(union, weight),
    avg_log_wage = weighted.mean(log(wage), weight),
    avg_log_wageotc = weighted.mean(log(wageotc), weight),
    avg_log_weekpay = weighted.mean(log(weekpay), weight),
    avg_age = weighted.mean(age, weight),
    pct_female = weighted.mean(sex == "Female", weight),
    pct_white = weighted.mean(race == "White", weight),
    pct_college_educated = weighted.mean(education_level %in% c("Bachelor's degree", "Master's degree", "Doctorate"), weight),
    pct_married = weighted.mean(married == "Married", weight),
    pct_metropolitan = weighted.mean(metropolitan == "Metropolitan", weight),
    .groups = "drop"
  ) %>%
  mutate(
    state_id = as.numeric(factor(state)),
    state = as.character(state)
  )

# Define pre/post periods and treatment unit
pre_years <- 1994:2011
treatment_year <- 2012
post_years <- 2013:2023

treated_state <- "WI"
treated_unit <- scm_data %>%
  distinct(state_id, state) %>%
  filter(state == treated_state) %>%
  pull(state_id)

# Define donor pool (control states)
union_friendly_states <- c("SD", "NM", "MO", "KS", "NE", "DE", "MT", 
                           "OH", "MI", "IL", "VT", "NH", "MA", "OR", "HI",
                           "MN", "WA", "ME", "CA", "RI", "PA", "NJ", "NY")

donor_units <- scm_data %>%
  distinct(state_id, state) %>%
  filter(state %in% union_friendly_states) %>%
  pull(state_id)

# ---------------------------------------------------------------
# 3.2 Synthetic Control Helper Function
# ---------------------------------------------------------------
run_scm <- function(outcome_var, pred_vars, special.predictors = NULL, treated_unit, donor_units, pre_years, treatment_year, post_years, data, y_label) {
  
  # Filter data to include only treated unit and donor units
  filtered_data <- data[data$state_id %in% c(treated_unit, donor_units), ]
  
  # Prepare data for SC
  dataprep.out <- dataprep(
    foo = filtered_data %>% as.data.frame(),
    predictors = pred_vars,
    predictors.op = "mean",
    dependent = outcome_var,
    unit.variable = "state_id",
    time.variable = "year",
    special.predictors = special.predictors,
    treatment.identifier = treated_unit,
    controls.identifier = donor_units,
    time.predictors.prior = pre_years,
    time.optimize.ssr = pre_years,
    time.plot = c(pre_years, treatment_year, post_years),
    unit.names.variable = "state"
  )
  
  # Compute custom weights
  custom_weights <- (1 / sapply(pred_vars, function(var) 
    sd(filtered_data[[var]][filtered_data$year %in% pre_years], na.rm = TRUE))) / 
    sum(1 / sapply(pred_vars, function(var) 
      sd(filtered_data[[var]][filtered_data$year %in% pre_years], na.rm = TRUE)))
  
  # Run Synth quietly (suppresses console output)
  synth.out <- quiet(synth(
    dataprep.out))
  
  # Extract tables
  synth.tables <- synth.tab(dataprep.res = dataprep.out, synth.res = synth.out)
  
  # Print predictor balance and unit weights
  cat("\n------------------------------------------\n")
  cat("Synthetic Control Results for:", outcome_var, "\n")
  cat("Unit Weights:\n")
  print(synth.tables$tab.w)
  cat("\nPredictor Balance:\n")
  print(synth.tables$tab.pred)
  cat("------------------------------------------\n\n")
  
  # Plot Path
  path.plot(
    dataprep.res = dataprep.out,
    synth.res = synth.out,
    tr.intake = 2011,
    Ylab = y_label,
    Xlab = "Year",
    Legend = c("Wisconsin","Synthetic Wisconsin"),
    Legend.position = "bottomleft"
  ) 
  
  # Plot Gaps
  gaps.plot(
    dataprep.res = dataprep.out,
    synth.res = synth.out,
    tr.intake = 2011,
    Ylab = paste("Gap (", treated_state, " - Synthetic)", sep = ""),
    Xlab = "Year"
  ) 
  
  # Compute treatment effect for the last post-treatment period
  observed_outcomes <- dataprep.out$Y1plot
  synthetic_outcomes <- dataprep.out$Y0plot %*% synth.out$solution.w
  
  # Identify the index for the last post-treatment period
  last_period_index <- which(dataprep.out$tag$time.plot == max(dataprep.out$tag$time.plot[dataprep.out$tag$time.plot > treatment_year]))
  
  # Calculate the treatment effect for the last period
  TE <- observed_outcomes[last_period_index] - synthetic_outcomes[last_period_index]
  
  cat("Treatment Effect in the last period:", round(TE, 4), "\n\n")
  
  # Placebo Tests
  placebos <- quiet(generate.placebos(dataprep.out, synth.out, strategy = "multisession"))
  
  # Plot placebo distribution
  plot_placebos(placebos) %>% plot()
  
  # Plot MSPE distribution (discard extremes for readability)
  mspe.plot(placebos) %>% plot()
  #mspe.plot(placebos, discard.extreme = TRUE, mspe.limit = 2) %>% plot()
  
  # Compute p-value of the placebo test
  p_val <- mspe.test(placebos)$p.val
  cat("Placebo Test p-value:", p_val, "\n\n")
  
  # Return results for further processing if needed
  list(
    TE = TE,
    p_value = p_val,
    synth_out = synth.out,
    dataprep_out = dataprep.out
  )
}

# ---------------------------------------------------------------
# 3.3 Running Synthetic Control on Multiple Outcomes (CPS)
# ---------------------------------------------------------------
variables_for_scm <- list(
  "avg_union_member" = "Union Membership Rate",
  "avg_union" = "Union Representation Rate",
  "avg_log_wage" = "Average Log of Hourly Earnings",
  "avg_log_wageotc" = "Average Log of Hourly Earnings (OTC)",
  "avg_log_weekpay" = "Average Log of Weekly Earnings (OTC)"
)

# Predictors
predictors_list <- c(
  "avg_age", 
  "pct_female", 
  "pct_white", 
  "pct_college_educated", 
  "pct_married", 
  "pct_metropolitan",
  "avg_union_member",
  "avg_log_wage",
  "avg_log_weekpay"
)

for (var in names(variables_for_scm)) {
  y_label <- variables_for_scm[[var]]
  run_scm(
    outcome_var = var,
    pred_vars = predictors_list,
    treated_unit = treated_unit,
    donor_units = donor_units,
    pre_years = pre_years,
    treatment_year = treatment_year,
    post_years = post_years,
    data = scm_data,
    y_label = y_label
  )
}

# ---------------------------------------------------------------
# Synthetic Control for Benefits Data
# ---------------------------------------------------------------
scm_benefits_data <- benefits_data %>%
  filter(public == 1) %>%
  group_by(state, year) %>%
  summarise(
    avg_log_pension_plan_at_work = weighted.mean(pension_plan_at_work, weight),
    avg_log_employer_insurance_contribution = weighted.mean(log(employer_insurance_contribution), weight),
    avg_health_insurance_included = weighted.mean(health_insurance_included, weight),
    avg_health_paid_by_employer = weighted.mean(health_paid_by_employer, weight),
    avg_health_paid_by_employer_full = weighted.mean(health_paid_by_employer_full, weight),
    
    avg_age = weighted.mean(age, weight),
    pct_female = weighted.mean(sex == "female", weight),
    pct_white = weighted.mean(race == "white", weight),
    pct_college_educated = weighted.mean(education_level %in% c("bachelor's degree", "master's degree", "doctorate"), weight),
    pct_married = weighted.mean(married == 1, weight),
    pct_metropolitan = weighted.mean(metropolitan == 1, weight),
    .groups = "drop"
  ) %>% 
  mutate(
    state_id = as.numeric(factor(state)),
    state = as.character(state)
  )

benefits_variables_for_scm <- list(
  "avg_log_pension_plan_at_work" = "Participation in Workplace Pension",
  "avg_log_employer_insurance_contribution" = "Log of Employer Health Insurance Contribution",
  "avg_health_insurance_included" = "Group Health Insurance Coverage"
)

predictors_list_benefits <- c(
  "avg_age", 
  "pct_female", 
  "pct_white", 
  "pct_college_educated", 
  "pct_married", 
  "pct_metropolitan",
  "avg_log_pension_plan_at_work",
  "avg_log_employer_insurance_contribution",
  "avg_health_insurance_included"
)

pre_years_benefits <- 1994:2011
treatment_year_benefits <- 2012
post_years_benefits <- 2013:2018

treated_state <- "WI"
treated_unit <- scm_benefits_data %>%
  distinct(state_id, state) %>%
  filter(state == treated_state) %>%
  pull(state_id)

# Define donor pool (control states)
union_friendly_states <- c("SD", "NM", "MO", "KS", "NE", "DE", "MT", 
                           "OH", "MI", "IL", "VT", "NH", "MA", "OR", "HI",
                           "MN", "WA", "ME", "CA", "RI", "PA", "NJ", "NY")

donor_units <- scm_benefits_data %>%
  distinct(state_id, state) %>%
  filter(state %in% union_friendly_states) %>%
  pull(state_id)

for (var in names(benefits_variables_for_scm)) {
  y_label <- benefits_variables_for_scm[[var]]
  
  run_scm(
    data = scm_benefits_data,
    outcome_var = var,
    pred_vars = predictors_list_benefits,
    treated_unit = treated_unit,
    donor_units = donor_units,
    pre_years = pre_years_benefits,
    treatment_year = treatment_year_benefits,
    post_years = post_years_benefits,
    y_label = y_label
  )
}
