
# ate for rct data
naive_ate <- function(dataset, treatment_var, outcome_var, method = "naive"){
  mean_df <- dataset %>%
    group_by_(treatment_var) %>%
    summarise_(y = paste("mean(", outcome_var, ")"),
               y_var = paste("var(", outcome_var, ")"),
               count = "n()") %>%
    mutate(y_var_weight = y_var/(count - 1))
  
  E_y0 = mean_df$y[mean_df$W == 0]
  E_y1 = mean_df$y[mean_df$W == 1]
  
  tau_hat <- E_y1 - E_y0
  se_hat <- sum(mean_df$y_var_weight) %>% sqrt()
  
  upper_ci <- tau_hat + se_hat*1.96
  lower_ci <- tau_hat - se_hat*1.96
  
  return(data.frame(Method = method, ATE = tau_hat, lower_ci = lower_ci, upper_ci = upper_ci))
}


# multiple regression
ate_condmean_ols <- function(dataset, treatment_var, outcome_var, method = "Direct Method") {
  reg_formula <- as.formula(paste(outcome_var, "~ ."))
  reg_result <- dataset %>% 
    lm(formula = reg_formula, data = .) %>% 
    summary()
  
  reg_coef <- reg_result %>% coef()
  
  tau_hat <- reg_coef[treatment_var,"Estimate"]
  se_hat <- reg_coef[treatment_var, "Std. Error"]
  upper_ci <- tau_hat + se_hat*1.96
  lower_ci <- tau_hat - se_hat*1.96
  
  return(data.frame(Method = method, ATE = tau_hat, lower_ci = lower_ci, upper_ci = upper_ci))
}



# Propensity Weighting
prop_score_weight <- function(dataset, p, treatment_var, outcome_var, covariates, method = "Propensity_Weighting"){
  wyp_df <- dataset %>% 
    mutate(p = p) %>%
    mutate_(tau_hat = paste("((", treatment_var, "- p)*", outcome_var,")/(p * (1-p))"),
            ps_er = paste(treatment_var, " - p"))
  
  d <- wyp_df[,covariates] * wyp_df$ps_er
  
  e <- data.frame(tau_hat = wyp_df$tau_hat, d) %>%
    lm(formula = tau_hat ~ ., data = .) %>%
    summary() %>%
    residuals()
  
  se_hat <- sqrt(mean(e^2))/sqrt(NROW(wyp_df))
  tau_hat <- wyp_df$tau_hat %>% mean
  upper_ci <- tau_hat + se_hat*1.96
  lower_ci <- tau_hat - se_hat*1.96
  
  return(data.frame(Method = method, ATE = tau_hat, lower_ci = lower_ci, upper_ci = upper_ci))
}


# Propensity Regression
prop_score_ols <- function(dataset, p, treatment_var, outcome_var, method = "Propensity_Regression") {
  
  reg_formula <- as.formula(paste(outcome_var, "~", treatment_var))
  
  ps_reg_res <- dataset %>%
    mutate(p = p) %>%
    mutate_(weights = paste("(", treatment_var, "/p ) + (( 1 -", treatment_var,") / (1 - p))")) %>%
    lm(data = ., formula = reg_formula, weights = weights) %>%
    summary()
  
  reg_coef <- ps_reg_res %>% coef()
  
  tau_hat <- reg_coef[treatment_var, "Estimate"]
  se_hat <- reg_coef[treatment_var, "Std. Error"]
  
  upper_ci <- tau_hat + se_hat*1.96
  lower_ci <- tau_hat - se_hat*1.96
  
  return(data.frame(Method = method, ATE = tau_hat, lower_ci = lower_ci, upper_ci = upper_ci))
}

# single equation lasso
ate_condmean_lasso <- function(dataset, treatment_var, outcome_var) {
  # Covariate names
  regs <- c(covariates, treatment_var)
  
  # glmnet requires inputs as matrices
  x <- as.matrix(dataset[regs])  
  y <- as.matrix(dataset[,outcome_var])
  
  # Set the penalty to betaw to be zero
  pfac <- c(rep(1,length(covariates)), 0) 
  
  # Call glmnet with alpha=1 is LASSO penalty
  model <- cv.glmnet(x, y, 
                     alpha=1, # 
                     penalty.factor=pfac) 
  
  # Automatically performs CV!
  betaw <- coef(model)[treatment_var,]
  return(data.frame(Method = "Single-equation LASSO", ATE = betaw, lower_ci = betaw, upper_ci = betaw))
}

# usual lasso
ate_lasso <- function(dataset, treatment_var, outcome_var) {
  # Covariate names
  regs <- c(covariates, treatment_var)
  
  # glmnet requires inputs as matrices
  x <- as.matrix(dataset[regs])  
  y <- as.matrix(dataset[,outcome_var])
  
  # Set the penalty to betaw to be zero
  pfac <- c(rep(1,length(covariates)), 1) 
  
  # Call glmnet with alpha=1 is LASSO penalty
  model <- cv.glmnet(x, y, 
                     alpha=1, # 
                     penalty.factor=pfac) 
  
  # Automatically performs CV!
  betaw <- coef(model)[treatment_var,]
  return(data.frame(Method = "Usual LASSO", ATE = betaw, lower_ci = betaw, upper_ci = betaw))
}

# propensity score by lasso logistic regression
prop_score_lasso <- function(dataset, treatment_var) {
  # glmnet requires inputs as matrices
  x <- as.matrix(dataset[covariates])  
  w <- as.matrix(dataset[,treatment_var])
  
  # Call glmnet with alpha=1 is LASSO penalty
  model <- cv.glmnet(x, w, 
                     alpha=1, 
                     family="binomial") 
  
  # Automatically performs CV
  p <- predict(model, newx=x, type="response")
  return(p)
}

#Doubly Robust Estimator with Random Forest
doubly_robust <- function(dataset, treatment_var, outcome_var, num_trees = 100, bootstrap_se = F) {
  
  #formula
  condmean_formula <- as.formula(paste(outcome_var, "~ ."))
  propensity_formula  <- as.formula(paste("I(factor(",treatment_var, ")) ~ . -", outcome_var))
  
  # Conditional mean 
  condmean <- glm(formula= condmean_formula, 
                  data=dataset,
                  family=binomial("logit"))
  tauhat1x <- dataset %>%
    mutate_(paste(treatment_var, "= 1")) %>%
    predict(condmean, type="response", newdata=.) %>%
    as.numeric()
  tauhat0x <- dataset %>%
    mutate_(paste(treatment_var, "= 0")) %>%
    predict(condmean, type="response", newdata=.) %>%
    as.numeric()
  
  # Propensity score (Will take ~1min to run)
  p <-randomForest(formula= propensity_formula, 
                   data=dataset,
                   ntree=num_trees,
                   type="classification",
                   seed=12325) %>%
    predict(., type="prob") %>% .[,2] %>% as.numeric()
  
  # Double robust estimator
  w <- dataset[,treatment_var]
  y <- dataset[,outcome_var]
  
  #clip p
  p <- ifelse( p == 0, min(p[p>0]), p)
  p <- ifelse( p == 1, max(p[p<1]), p)
  #tau_hat
  est1 <- w*(y - tauhat1x)/p + (1-w)*(y - tauhat0x)/(1-p)
  est2 <- tauhat1x - tauhat0x
  tau_hat <- mean(est1, na.rm = TRUE) + mean(est2)
  
  if(bootstrap_se == T){
    #se_hat Bootstrap
    B <- 1000
    Boot_result <- c()
    for(i in 1:B){
      Boot_result[i] <- tau_hat_dr_est(w,y,p,tauhat0x,tauhat1x)
    }
    se_hat <- sd(Boot_result)
  }else{
    #se_hat from law of large number(sandwitch estimator)
    Ii <- (w*y)/p - tauhat1x*(w-p)/p - ( ((1-w)*y / (1-p)) + (tauhat0x*(w - p) / (1 - p))  ) - tau_hat
    se_hat <- sqrt(sum(Ii^2)*(length(Ii)^(-2)))
  }
  
  #confidence interval
  upper_ci <- tau_hat + se_hat*1.96
  lower_ci <- tau_hat - se_hat*1.96
  return(data.frame(Method = "Doubly Robust with Random Forest PS", 
                    ATE = tau_hat, lower_ci = lower_ci, upper_ci = upper_ci))
}


#Doubly Robust Estimator with glm
doubly_robust_glm <- function(dataset, treatment_var, outcome_var, bootstrap_se = F) {
  
  #formula
  condmean_formula <- as.formula(paste(outcome_var, "~ ."))
  propensity_formula  <- as.formula(paste("I(factor(",treatment_var, ")) ~ . -", outcome_var))
  
  # Conditional mean 
  condmean <- glm(formula= condmean_formula, 
                  data=dataset,
                  family=binomial("logit"))
  tauhat1x <- dataset %>%
    mutate(W = 1) %>%
    predict(condmean, type="response", newdata=.) %>%
    as.numeric()
  tauhat0x <- dataset %>%
    mutate(W = 0) %>%
    predict(condmean, type="response", newdata=.) %>%
    as.numeric()
  
  # Propensity score (Will take ~1min to run)
  p <- glm(formula= propensity_formula, 
          data=dataset,
          family=binomial("logit")) %>%
    predict(., type="response") %>% as.numeric()
  
  # Double robust estimator
  w <- dataset[,treatment_var]
  y <- dataset[,outcome_var]
  
  #tau_hat
  est1 <- w*(y - tauhat1x)/p + (1-w)*(y - tauhat0x)/(1-p)
  est2 <- tauhat1x - tauhat0x
  tau_hat <- mean(est1, na.rm = TRUE) + mean(est2)
  
  if(bootstrap_se == T){
    #se_hat Bootstrap
    B <- 1000
    Boot_result <- c()
    for(i in 1:B){
      Boot_result[i] <- tau_hat_dr_est(w,y,p,tauhat0x,tauhat1x)
    }
    se_hat <- sd(Boot_result)
  }else{
    #se_hat from law of large number(sandwitch estimator)
    Ii <- (w*y)/p - tauhat1x*(w-p)/p - ( ((1-w)*y / (1-p)) + (tauhat0x*(w - p) / (1 - p))  ) - tau_hat
    se_hat <- sqrt(sum(Ii^2)*(length(Ii)^(-2)))
  }
  
  #confidence interval
  upper_ci <- tau_hat + se_hat*1.96
  lower_ci <- tau_hat - se_hat*1.96
  return(data.frame(Method = "Doubly Robust with logistic regression PS", 
                    ATE = tau_hat, lower_ci = lower_ci, upper_ci = upper_ci))
}

# Bootstrap Standard error for doubly robust
tau_hat_dr_est <- function(w, y, p, tauhat0x, tauhat1x){
  #sampling with replacement
  bootstrap_sample_idx <- sample(length(w), length(w), replace = T)
  
  #resample each variable
  w_B <- w[bootstrap_sample_idx]
  y_B <- y[bootstrap_sample_idx]
  p_B <- p[bootstrap_sample_idx]
  tauhat0x_B <- tauhat0x[bootstrap_sample_idx]
  tauhat1x_B <- tauhat1x[bootstrap_sample_idx]
  
  #tau_hat
  est1 <- w_B*(y_B - tauhat1x_B)/p_B + (1-w_B)*(y_B - tauhat0x_B)/(1-p_B)
  est2 <- tauhat1x_B - tauhat0x_B
  tau_hat_dr <- mean(est1, na.rm = TRUE) + mean(est2)
  return(tau_hat_dr)
}

# Belloni et.al(2013)
belloni <- function(dataset, treatment_var, outcome_var, method = "Belloni et.al") {
  
  # Creating all pairwise interaction terms
  newcovs <- covariates
  for (c1 in covariates) {
    for (c2 in covariates) {
      newc_name <- paste(c1, c2, sep="")
      dataset[,newc_name] <- dataset[,c1]*dataset[,c2]
      newcovs <- c(newcovs, newc_name)
    }
  }
  
  # glmnet requires inputs as matrices
  x <- as.matrix(cbind(dataset[,newcovs]))
  w <- as.matrix(dataset[treatment_var])
  y <- as.matrix(dataset[outcome_var])
  
  # Call glmnet with alpha=1 is LASSO penalty
  model_xw <- cv.glmnet(x, w,  alpha=1)
  model_xy <- cv.glmnet(x, y,  alpha=1)
  
  # Grab coefficients
  c_xw <- coef(model_xw, s=model_xw$lambda.min)
  c_xy <- coef(model_xy, s=model_xw$lambda.min)
  
  # Nonzero coefficients
  c_xw_nonzero <- which(c_xw[2:(length(newcovs)+1)] > 0)
  c_xy_nonzero <- which(c_xy[2:(length(newcovs)+1)] > 0)
  c_union <- unique(c(c_xw_nonzero, c_xy_nonzero)) - 1
  
  # Restricted
  x_restricted <- cbind(x[,c_union], w)
  
  # OLS on resulting regressor matrix
  reg_coef <- lm(y ~ x_restricted) %>% summary() %>% coef()
  #betaw <- as.numeric(coef(post_ols)["x_restrictedW"])
  tau_hat <- reg_coef["x_restrictedW","Estimate"]
  se_hat <- reg_coef["x_restrictedW", "Std. Error"]
  upper_ci <- tau_hat + se_hat*1.96
  lower_ci <- tau_hat - se_hat*1.96
  
  return(data.frame(Method = method, ATE = tau_hat, lower_ci = lower_ci, upper_ci = upper_ci))
}


# Double Machine Learning predict W and Y
chernozhukov <- function(dataset, treatment_var, outcome_var, idx1, idx2, num_trees) {
  
  #formula
  predict_w <- as.formula(paste("I(factor(",treatment_var, ")) ~ . -", outcome_var))
  predict_y <- as.formula(paste("I(factor(",outcome_var, ")) ~ . -", treatment_var))
  
  
  # Estimate each submodel separately, on its own half of the data
  rf1 <- randomForest(formula= predict_w, 
                      data=dataset[idx1,], # Only on one half
                      ntree=num_trees,
                      type="classification",
                      seed = 123) 
  rf2 <- randomForest(formula= predict_y, 
                      data=dataset[idx2,], # Only on the other half
                      ntree=num_trees,
                      type="classification",
                      seed = 123)
  
  # Predict and residualize
  EWhat <- dataset %>% 
    predict(rf1, newdata=., type="prob") %>%
    .[,2] %>% as.numeric()
  EYhat <- dataset %>% 
    predict(rf2, newdata=., type="prob") %>%
    .[,2] %>% as.numeric()
  
  W_resid <- dataset[,treatment_var] - EWhat
  Y_resid <- dataset[,outcome_var] - EYhat
  
  # Linear regression
  reg_coef <- lm(Y_resid ~ 0 + W_resid) %>% summary() %>% coef
  
  tau_hat <- reg_coef["W_resid","Estimate"]
  se_hat <- reg_coef["W_resid", "Std. Error"]
  
  return(list(tau_hat = tau_hat, se_hat = se_hat))
}

# Double Machine Learning
double_ml <- function(dataset, treatment_var, outcome_var, num_trees = 100, method = "Double Machine Learning") {
  # Splits sample
  N <- dim(dataset)[1]
  idx1 <- 1:floor(N/2)
  idx2 <- (floor(N/2)+1):N
  
  # Apply the algorithm in each half, then swaps them
  betaw1 <- chernozhukov(dataset, treatment_var, outcome_var, idx1, idx2, num_trees)
  betaw2 <- chernozhukov(dataset, treatment_var, outcome_var, idx2, idx1, num_trees) # Swaps halves
  
  tau_hat <- mean(c(betaw1$tau_hat, betaw2$tau_hat))
  se_hat <- mean(c(betaw1$se_hat, betaw2$se_hat))
  
  upper_ci <- tau_hat + se_hat*1.96
  lower_ci <- tau_hat - se_hat*1.96
  
  return(data.frame(Method = method, ATE = tau_hat, lower_ci = lower_ci, upper_ci = upper_ci))
}


# approximate residual balancing
residual_balance_ATE <- function(dataset, treatment_var, outcome_var, optimizer = "quadprog", method = "residual_balancing"){
  tauhat_balance <- balanceHD::residualBalance.ate(X = df_mod[,covariates],
                                                   Y = df_mod[,outcome_var],
                                                   W = df_mod[,treatment_var],
                                                   estimate.se = T,
                                                   optimizer = optimizer)
  
  return(data.frame(Method = "residual_balancing",
                                 ATE = tauhat_balance[1],
                                 lower_ci = tauhat_balance[1] - 1.96*tauhat_balance[2],
                                 upper_ci = tauhat_balance[1] + 1.96*tauhat_balance[2]))
  
}
