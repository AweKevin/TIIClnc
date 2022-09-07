library(survival)
library(randomForestSRC)
library(CoxBoost)
library(glmnet)
library(MASS)
library(plsRcox)
library(gbm)
library(superpc)
library(survivalsvm)


# RSF --------------------------------------------------------------


RSF.mod <- function(train) {
  mod.RSF <- rfsrc(Surv(time, status) ~ .,
    data = train,
    ntree = 1000,
    splitrule = "logrank",
    importance = T,
    proximity = T,
    forest = T,
    seed = 123456
  )
  return(mod.RSF)
}


# CoxBoost --------------------------------------------------------------


CoxBoost.mod <- function(train) {
  time <- train$time
  status <- train$status
  x <- train %>%
    dplyr::select(-c(1:2)) %>%
    as.matrix()
  # determine penalty parameter
  
  optim.CoxBoost <- optimCoxBoostPenalty(
    time = time, status = status, x = x,
    trace = T, start.penalty = 500
  )
  mod.CoxBoost <- CoxBoost(
    time = time, status = status, x = x,
    stepno = optim.CoxBoost$cv.res$optimal.step,
    penalty = optim.CoxBoost$penalty
  )
  return(mod.CoxBoost)
}


# Ridge --------------------------------------------------------------


Ridge.mod <- function(train) {
  x2 <- as.matrix(train %>% dplyr::select(-c(1:2)))
  y2 <- data.matrix(Surv(train$time, as.factor(train$status)))
  
  cv.fit <- cv.glmnet(x2, y2,
    nfolds = 10,
    family = "cox", 
    grouped = FALSE, 
    alpha = 0, 
    type.measure = "mse"
  )
  # plot(cv.fit)
  fit <- glmnet(x2, y2, family = "cox", alpha = 0)
  # plot(fit)
  mod.Ridge <- list(cv.fit, fit)
  return(mod.Ridge)
}


# Lasso --------------------------------------------------------------


Lasso.mod <- function(train) {
  x2 <- as.matrix(train %>% dplyr::select(-c(1:2)))
  y2 <- data.matrix(Surv(train$time, as.factor(train$status)))
  
  cv.fit <- cv.glmnet(x2, y2,
    nfolds = 10,
    family = "cox", 
    grouped = FALSE, 
    alpha = 1, 
    type.measure = "mse"
  )
  # plot(cv.fit)
  fit <- glmnet(x2, y2, family = "cox", alpha = 1)
  # plot(fit)
  mod.Lasso <- list(cv.fit, fit)
  return(mod.Lasso)
}


# Enet --------------------------------------------------------------


Enet.mod <- function(train) {
  alpha <- 0.1 # 0-1
  x2 <- as.matrix(train %>% dplyr::select(-c(1:2)))
  y2 <- data.matrix(Surv(train$time, as.factor(train$status)))
  
  cv.fit <- cv.glmnet(x2, y2,
                      nfolds = 10,
                      family = "cox", 
                      grouped = FALSE, 
                      alpha = alpha, 
                      type.measure = "mse"
  )
  # plot(cv.fit)
  fit <- glmnet(x2, y2, family = "cox", alpha = alpha)
  # plot(fit)
  mod.Enet <- list(cv.fit, fit)
  return(mod.Enet)
}


# stepwiseCox.forward --------------------------------------------------------------


stepwiseCox.forward.mod <- function(train) {
  fit0 <- coxph(Surv(time, status) ~ ., data = train)
  fit_forward <- stepAIC(fit0, direction = "forward")
  AIC_forward <- extractAIC(fit_forward) %>% .[2]
  mod.stepwiseCox.forward <- fit_forward
  return(mod.stepwiseCox.forward)
}


# stepwiseCox.both --------------------------------------------------------------


stepwiseCox.both.mod <- function(train) {
  fit0 <- coxph(Surv(time, status) ~ ., data = train)
  fit_both <- stepAIC(fit0, direction = "both")
  AIC_both <- extractAIC(fit_both) %>% .[2]
  mod.stepwiseCox.both <- fit_both
  return(mod.stepwiseCox.both)
}


# stepwiseCox.backward --------------------------------------------------------------


stepwiseCox.backward.mod <- function(train) {
  fit0 <- coxph(Surv(time, status) ~ ., data = train)
  fit_backward <- stepAIC(fit0, direction = "backward")
  AIC_backward <- extractAIC(fit_backward) %>% .[2]
  mod.stepwiseCox.backward <- fit_backward
  return(mod.stepwiseCox.backward)
}


# plsRcox --------------------------------------------------------------


plsRcox.mod <- function(train) {
  time <- train$time
  status <- train$status
  x <- train %>%
    dplyr::select(-c(1:2)) %>%
    as.matrix()
  gene.names <- colnames(x)
  
  plsRcox.cv <- cv.plsRcox(data = list(x = x, time = time, status = status), nfold = 10, nt = 10)
  nt <- which.max(plsRcox.cv$cv.error5[-1])
  mod.plsRcox <- plsRcox(x, time = time, event = status, nt = nt, sparse = F)
  mod.plsRcox <- list(mod.plsRcox, gene.names)
  return(mod.plsRcox)
}


# GBM --------------------------------------------------------------


GBM.mod <- function(train) {
  
  mod.GBM <- gbm(Surv(time, status) ~ .,
    distribution = "coxph",
    data = train,
    n.trees = 1000, shrinkage = 0.01, cv.folds = 10, verbose = FALSE
  )
  # which.min(mod.GBM$cv.error)
  return(mod.GBM)
}


# SuperPC --------------------------------------------------------------


SuperPC.mod <- function(train) {
  x2 <- train %>%
    dplyr::select(-c(1:2)) %>%
    t()
  y2 <- train$time
  censoring.status <- train$status
  featurenames <- colnames(x2)
  data.train <- list(
    x = x2,
    y = y2,
    censoring.status = censoring.status,
    featurenames = featurenames
  )
  obj.SuperPC <- superpc.train(data.train, type = "survival")
  
  cv.SuperPC <- superpc.cv(obj.SuperPC,
    data.train,
    n.components = 1,
    n.threshold = 20,
    n.fold = 10
  )
  sel <- cv.SuperPC$scor %>% which.max()
  threshold <- cv.SuperPC$thresholds[sel]
  mod.SuperPC <- list(obj.SuperPC, threshold, data.train)
  return(mod.SuperPC)
}


# survivalSVM --------------------------------------------------------------


survivalSVM.mod <- function(train) {
  mod.survivalSVM <- survivalsvm(Surv(time, status) ~ .,
    data = train,
    gamma.mu = 0.5,
    type = "regression",
    opt.meth = "ipop", 
  )
  return(mod.survivalSVM)
}
