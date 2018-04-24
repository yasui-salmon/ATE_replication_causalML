library(tidyverse)
library(glmnet)
library(randomForest)
library(grf)

source("ate_functions.R")

# Set seed for reproducibility
set.seed(1991)

# Load data
data_raw <- read.csv('socialpresswgeooneperhh_NEIGH.csv')

# These are the covariates we'll use
cts_variables_names <- c("yob", "city", "hh_size", "totalpopulation_estimate",
                         "percent_male", "median_age",
                         "percent_62yearsandover",
                         "percent_white", "percent_black",
                         "percent_asian", "median_income",
                         "employ_20to64", "highschool", "bach_orhigher",
                         "percent_hispanicorlatino")
binary_variables_names <- c("sex","g2000", "g2002", "p2000", "p2002", "p2004")
covariates <- c(cts_variables_names, binary_variables_names)
all_variables_names <- c(covariates, "outcome_voted", "treat_neighbors")

# We will not use all observations -- it would take too long to run all the methods below
n_obs <- 50000

# Selecting only desired covariates
data_subset <- data_raw %>%
  sample_n(size = n_obs) %>%
  dplyr::select(all_variables_names)


# Extracting and scaling continuous variables
scaled_cts_covariates <- data_subset %>%
  dplyr::select(cts_variables_names) %>%
  scale()

# Extracting indicator variables
binary_covariates <- data_subset %>%
  dplyr::select(binary_variables_names)

# Extracting outcome and treatment
outcome <- data_subset %>% dplyr::select(outcome_voted)
treatment <- data_subset %>% dplyr::select(treat_neighbors)

# Setting up the data, renaming columns and discarding rows with NA (if any)
df <- data.frame(scaled_cts_covariates, binary_covariates, outcome, treatment) %>%
  plyr::rename(c(treat_neighbors = "W",
                 outcome_voted = "Y")) %>%
  na.omit()


#ground truth
ate_oracle <- naive_ate(dataset = df, treatment_var = "W", outcome_var = "Y", method = "oracle")
ggplot(ate_oracle, aes(y = ATE, x = Method, color = Method)) + 
  geom_pointrange(aes(ymax = upper_ci, ymin = lower_ci))

#introduce sampling bias
pt <- .85 # Drop p% of voters who satisfy the following condition
pc <- .85

# These individuals are likely TO GO voting: drop from TREATMENT
drop_from_treat <-  (df[,"g2000"]==1 | df[,"g2002"]==1) |
  (df[,"p2000"]==1 | df[,"p2002"]==1 | df[,"p2002"] == 1) |
  (df[,"city"] > 2) | (df[,"yob"] > 2)

# These individuals are likely NOT TO GO voting: drop from CONTROL
drop_from_control <-(df[,"g2000"]==0 | df[,"g2002"] == 0) |
  (df[,"p2000"]==0 | df[,"p2002"]==0 | df[,"p2004"]==0) |
  (df[,"city"] < -2 | df[,"yob"] < -2) 


drop_treat_idx <- which(df[,"W"] == 1 & drop_from_treat)
drop_control_idx <- which(df[,"W"] == 0 & drop_from_control)

drop_idx <- unique(c(drop_treat_idx[1:round(pt*length(drop_treat_idx))],
                     drop_control_idx[1:round(pc*length(drop_control_idx))]))

print(length(drop_idx))

df_mod <- df[-drop_idx,]



# Computing the propensity score by logistic regression of W on X.
p_logistic <- df_mod %>% 
  dplyr::select(covariates, W) %>%   
  glm(W ~ ., data = ., family = binomial(link = "logit")) %>%   
  predict(type= "response")

hist(p_logistic)

#Observational Method
naive_ate <- naive_ate(dataset = df_mod, treatment_var = "W", outcome_var = "Y")

ggplot(rbind(ate_oracle, naive_ate), aes(y = ATE, x = Method, color = Method)) + 
  geom_pointrange(aes(ymax = upper_ci, ymin = lower_ci)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




#Traditional Econometrics Method
tauhat_naive_mod <- ate_condmean_ols(df_mod, treatment_var = "W", outcome_var = "Y")
ggplot(rbind(ate_oracle, naive_ate, tauhat_naive_mod), aes(y = ATE, x = Method, color = Method)) + 
  geom_pointrange(aes(ymax = upper_ci, ymin = lower_ci)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Propensity Score Weighting
tauhat_psw <- prop_score_weight(dataset = df_mod, p = p_logistic, treatment_var = "W", outcome_var = "Y")
ggplot(rbind(ate_oracle, naive_ate, tauhat_naive_mod, tauhat_psw), aes(y = ATE, x = Method, color = Method)) + 
  geom_pointrange(aes(ymax = upper_ci, ymin = lower_ci)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Propensity Regression
tauhat_ols <- prop_score_ols(dataset = df_mod, p = p_logistic, treatment_var = "W", outcome_var = "Y")

#LASSO Propensity Score
p_lasso <- prop_score_lasso(df_mod, treatment_var = "W")
tauhat_prop_lasso <- prop_score_weight(dataset = df_mod, 
                                       p = p_lasso[,1], 
                                       treatment_var = "W", 
                                       outcome_var = "Y",
                                       method = "Propensity_Weighting_LASSOPS")

ggplot(rbind(ate_oracle, naive_ate, tauhat_naive_mod, tauhat_psw, tauhat_ols, tauhat_prop_lasso), aes(y = ATE, x = Method, color = Method)) + 
  geom_pointrange(aes(ymax = upper_ci, ymin = lower_ci)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Single-equation LASSO
tauhat_lasso <- ate_condmean_lasso(df_mod, treatment_var = "W", outcome_var = "Y")
ggplot(rbind(ate_oracle, naive_ate, tauhat_naive_mod, tauhat_psw, tauhat_ols, tauhat_prop_lasso, tauhat_lasso), aes(y = ATE, x = Method, color = Method)) + 
  geom_pointrange(aes(ymax = upper_ci, ymin = lower_ci)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#all reguralized
tauhat_lasso_all <- ate_lasso(df_mod, treatment_var = "W", outcome_var = "Y")

ggplot(rbind(ate_oracle, naive_ate, tauhat_naive_mod, tauhat_psw, tauhat_ols, tauhat_prop_lasso, tauhat_lasso, tauhat_lasso_all), aes(y = ATE, x = Method, color = Method)) + 
  geom_pointrange(aes(ymax = upper_ci, ymin = lower_ci)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##ML Method
#Doubly Robust with Random Forest
tauhat_doubly <- doubly_robust(df_mod, treatment_var = "W", outcome_var = "Y", 2500)
tauhat_doubly_glm <- doubly_robust_glm(df_mod, treatment_var = "W", outcome_var = "Y")
ggplot(rbind(ate_oracle, naive_ate, tauhat_naive_mod, tauhat_psw, 
             tauhat_ols, tauhat_lasso, tauhat_prop_lasso, tauhat_doubly,
             tauhat_doubly_glm), aes(y = ATE, x = Method, color = Method)) + 
  geom_pointrange(aes(ymax = upper_ci, ymin = lower_ci)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Belloni et al (2013)
tauhat_belloni <- belloni(df_mod, treatment_var = "W", outcome_var = "Y")
ggplot(rbind(ate_oracle, naive_ate, tauhat_naive_mod, tauhat_psw, 
             tauhat_ols, tauhat_lasso, tauhat_prop_lasso,tauhat_lasso_all,  tauhat_doubly,
             tauhat_doubly_glm, tauhat_belloni), aes(y = ATE, x = Method, color = Method)) + 
  geom_pointrange(aes(ymax = upper_ci, ymin = lower_ci)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#double machine learning
tauhat_double_ml <- double_ml(df_mod, treatment_var = "W", outcome_var = "Y", num_tree = 2000)
ggplot(rbind(ate_oracle, naive_ate, tauhat_naive_mod, tauhat_psw, 
             tauhat_ols, tauhat_lasso, tauhat_prop_lasso, tauhat_doubly,
             tauhat_doubly_glm, tauhat_belloni, tauhat_double_ml), aes(y = ATE, x = Method, color = Method)) + 
  geom_pointrange(aes(ymax = upper_ci, ymin = lower_ci)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#balance HD
library(balanceHD)
tauhat_balance <- balanceHD::residualBalance.ate(X = df_mod[,covariates],
                                                 Y = df_mod[,"Y"],
                                                 W = df_mod[,"W"],
                                                 estimate.se = T,
                                                 optimizer = "pogs")
print(tauhat_balance)
tauhat_balancdHD <- data.frame(Method = "residual_balancing",
                               ATE = tauhat_balance[1],
                               lower_ci = tauhat_balance[1] - 1.96*tauhat_balance[2],
                               upper_ci = tauhat_balance[1] + 1.96*tauhat_balance[2])

ggplot(rbind(ate_oracle, naive_ate, tauhat_naive_mod, tauhat_psw, 
             tauhat_ols, tauhat_lasso, tauhat_prop_lasso, tauhat_doubly,
             tauhat_doubly_glm, tauhat_belloni, tauhat_double_ml,
             tauhat_balancdHD), aes(y = ATE, x = Method, color = Method)) + 
  geom_pointrange(aes(ymax = upper_ci, ymin = lower_ci)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#grf
# Fitting a causal forest
forest <- grf::causal_forest(X=as.matrix(df_mod[,covariates]),
                             Y=as.matrix(df_mod[,"Y"]),
                             W=as.matrix(df_mod[,"W"]), 
                             num.trees = 2000,
                             honesty=TRUE,
                             seed=12345)


# Incorrect way to derive ATE and its standard errors
pred <- predict(forest, estimate.variance = TRUE)
ate_bad <- mean(pred$predictions)
se_bad <- sqrt(mean(pred$variance.estimates))
cat(sprintf("Incorrect ATE: %1.3f (SE: %1.3f)", ate_bad, se_bad))

# Doubly-robust ATE 
ate_cf_robust <- grf::estimate_average_effect(forest)
print(ate_cf_robust)
cat(sprintf("Doubly robust ATE: %1.3f (SE: %1.3f)", ate_cf_robust["estimate"],
            ate_cf_robust["std.err"]))

tauhat_cf <- data.frame(Method = "Causal Forest(GRF)", 
                        ATE = ate_cf_robust["estimate"], 
                        lower_ci = ate_cf_robust["estimate"] - (ate_cf_robust["std.err"]*1.96),
                        upper_ci = ate_cf_robust["estimate"] + (ate_cf_robust["std.err"]*1.96))

ggplot(rbind(ate_oracle, naive_ate, tauhat_naive_mod, tauhat_psw, 
             tauhat_ols, tauhat_lasso, tauhat_prop_lasso, tauhat_doubly,
             tauhat_doubly_glm, tauhat_belloni,
             tauhat_balancdHD, tauhat_cf), aes(y = ATE, x = Method, color = Method)) + 
  geom_pointrange(aes(ymax = upper_ci, ymin = lower_ci)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
