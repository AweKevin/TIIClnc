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
  set.seed(123456)
  optim.CoxBoost <- optimCoxBoostPenalty(
    time = time, status = status, x = x,
    trace = T, start.penalty = 500
  )
  # Fit with obtained penalty parameter and optimal number of boosting
  # steps obtained by cross-validation
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
  set.seed(123456)
  cv.fit <- cv.glmnet(x2, y2,
    # nfolds = nrow(x2), 
    nfolds = 10,
    family = "cox", # cox
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
  set.seed(123456)
  cv.fit <- cv.glmnet(x2, y2,
    # nfolds = nrow(x2), 
    nfolds = 10,
    family = "cox", # cox
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


# Enet.1 --------------------------------------------------------------


Enet.1.mod <- function(train) {
  alpha <- 0.1
  x2 <- as.matrix(train %>% dplyr::select(-c(1:2)))
  y2 <- data.matrix(Surv(train$time, as.factor(train$status)))
  set.seed(123456)
  cv.fit <- cv.glmnet(x2, y2,
                      # nfolds = nrow(x2), 
                      nfolds = 10,
                      family = "cox", # cox
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


# Enet.2 --------------------------------------------------------------


Enet.2.mod <- function(train) {
  alpha <- 0.2
  x2 <- as.matrix(train %>% dplyr::select(-c(1:2)))
  y2 <- data.matrix(Surv(train$time, as.factor(train$status)))
  set.seed(123456)
  cv.fit <- cv.glmnet(x2, y2,
                      # nfolds = nrow(x2), 
                      nfolds = 10,
                      family = "cox", # cox
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


# Enet.3 --------------------------------------------------------------


Enet.3.mod <- function(train) {
  alpha <- 0.3
  x2 <- as.matrix(train %>% dplyr::select(-c(1:2)))
  y2 <- data.matrix(Surv(train$time, as.factor(train$status)))
  set.seed(123456)
  cv.fit <- cv.glmnet(x2, y2,
                      # nfolds = nrow(x2), 
                      nfolds = 10,
                      family = "cox", # cox
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


# Enet.4 --------------------------------------------------------------


Enet.4.mod <- function(train) {
  alpha <- 0.4
  x2 <- as.matrix(train %>% dplyr::select(-c(1:2)))
  y2 <- data.matrix(Surv(train$time, as.factor(train$status)))
  set.seed(123456)
  cv.fit <- cv.glmnet(x2, y2,
                      # nfolds = nrow(x2), 
                      nfolds = 10,
                      family = "cox", # cox
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


# Enet.5 --------------------------------------------------------------


Enet.5.mod <- function(train) {
  alpha <- 0.5
  x2 <- as.matrix(train %>% dplyr::select(-c(1:2)))
  y2 <- data.matrix(Surv(train$time, as.factor(train$status)))
  set.seed(123456)
  cv.fit <- cv.glmnet(x2, y2,
                      # nfolds = nrow(x2), 
                      nfolds = 10,
                      family = "cox", # cox
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


# Enet.6 --------------------------------------------------------------


Enet.6.mod <- function(train) {
  alpha <- 0.6
  x2 <- as.matrix(train %>% dplyr::select(-c(1:2)))
  y2 <- data.matrix(Surv(train$time, as.factor(train$status)))
  set.seed(123456)
  cv.fit <- cv.glmnet(x2, y2,
                      # nfolds = nrow(x2), 
                      nfolds = 10,
                      family = "cox", # cox
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


# Enet.7 --------------------------------------------------------------


Enet.7.mod <- function(train) {
  alpha <- 0.7
  x2 <- as.matrix(train %>% dplyr::select(-c(1:2)))
  y2 <- data.matrix(Surv(train$time, as.factor(train$status)))
  set.seed(123456)
  cv.fit <- cv.glmnet(x2, y2,
                      # nfolds = nrow(x2), 
                      nfolds = 10,
                      family = "cox", # cox
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


# Enet.8 --------------------------------------------------------------


Enet.8.mod <- function(train) {
  alpha <- 0.8
  x2 <- as.matrix(train %>% dplyr::select(-c(1:2)))
  y2 <- data.matrix(Surv(train$time, as.factor(train$status)))
  set.seed(123456)
  cv.fit <- cv.glmnet(x2, y2,
                      # nfolds = nrow(x2), 
                      nfolds = 10,
                      family = "cox", # cox
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


# Enet.9 --------------------------------------------------------------


Enet.9.mod <- function(train) {
  alpha <- 0.9
  x2 <- as.matrix(train %>% dplyr::select(-c(1:2)))
  y2 <- data.matrix(Surv(train$time, as.factor(train$status)))
  set.seed(123456)
  cv.fit <- cv.glmnet(x2, y2,
                      # nfolds = nrow(x2), 
                      nfolds = 10,
                      family = "cox", # cox
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
  set.seed(123456)
  plsRcox.cv <- cv.plsRcox(data = list(x = x, time = time, status = status), nfold = 10, nt = 10)
  nt <- which.max(plsRcox.cv$cv.error5[-1])
  mod.plsRcox <- plsRcox(x, time = time, event = status, nt = nt, sparse = F)
  mod.plsRcox <- list(mod.plsRcox, gene.names)
  return(mod.plsRcox)
}


# GBM --------------------------------------------------------------


GBM.mod <- function(train) {
  set.seed(123456)
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
  set.seed(123456)
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
    opt.meth = "ipop", kernel = "add_kernel"
  )
  return(mod.survivalSVM)
}



ML <- c("RSF", "CoxBoost", "Ridge", "Lasso", "Enet.1", "Enet.2", "Enet.3", "Enet.4", "Enet.5", "Enet.6", "Enet.7", "Enet.8", "Enet.9", "stepwiseCox.forward", "stepwiseCox.both", "stepwiseCox.backward", "plsRcox", "GBM", "SuperPC", "survivalSVM")
fs <- c("Lasso", "stepwiseCox.forward", "stepwiseCox.both", "stepwiseCox.backward", "CoxBoost", "RSF")
pm <- c("Ridge", "Enet.1", "Enet.2", "Enet.3", "Enet.4", "Enet.5", "Enet.6", "Enet.7", "Enet.8", "Enet.9", "plsRcox", "GBM", "SuperPC", "survivalSVM")

comb <- expand.grid(ML, ML)
names(comb) <- c("ML1", "ML2")
comb$ML1 %<>% as.character()
comb$ML2 %<>% as.character()
comb <- comb %>% dplyr::filter(ML1 %in% fs) %>% dplyr::filter(ML1 != ML2)
comb <- comb %>% dplyr::filter(!(ML1 %>% str_detect("stepwiseCox") & ML2 %>% str_detect("stepwiseCox")))
comb <- comb %>% dplyr::filter(ML1 != "stepwiseCox.forward")
comb1 <- comb %>% dplyr::filter(ML1 != "Lasso")
comb2 <- comb %>% dplyr::filter(ML1 == "Lasso")
comb2 <- comb2 %>% dplyr::filter(!ML2 %>% str_detect("Enet")) %>% dplyr::filter(ML2 != "Ridge")
comb <- rbind(comb1, comb2) %>% arrange(ML1)

FUN_comb_ML <- function(i) {
  ML1 <- comb$ML1[i]
  ML2 <- comb$ML2[i]
  function1 <- get(str_c(ML1, ".fs"))
  function2 <- get(str_c(ML2, ".mod"))
  function3 <- get(str_c(ML2, ".pred"))
  genes <- function1(TCGA)
  TCGA.2 <- TCGA %>% dplyr::select(c("time", "status", genes))
  CGGA.2 <- CGGA %>% dplyr::select(c("time", "status", genes))
  xiangya.2 <- xiangya %>% dplyr::select(c("time", "status", genes))
  GSE108474.2 <- GSE108474 %>% dplyr::select(c("time", "status", genes))
  model <- function2(TCGA.2)
  cindex.TCGA <- model %>% function3(test = TCGA.2) %>% cindex()
  cindex.CGGA <- model %>% function3(test = CGGA.2) %>% cindex()
  cindex.xiangya <- model %>% function3(test = xiangya.2) %>% cindex()
  cindex.GSE108474 <- model %>% function3(test = GSE108474.2) %>% cindex()
  tmp_df <- data.frame(method = str_c(ML1, ML2, sep = " + "), TCGA = cindex.TCGA, CGGA = cindex.CGGA, xiangya = cindex.xiangya, GSE108474 = cindex.GSE108474)
  return(tmp_df)
}

res_l1 <- pbapply::pblapply(1:nrow(comb), FUN = function(i) {try(FUN_comb_ML(i), TRUE)})
res_df1 <- do.call(rbind, res_l1)


FUN_single_ML <- function(i) {
  function1 <- get(str_c(i, ".mod"))
  function2 <- get(str_c(i, ".pred"))
  model <- function1(TCGA)
  cindex.TCGA <- model %>% function2(test = TCGA) %>% cindex()
  cindex.CGGA <- model %>% function2(test = CGGA) %>% cindex()
  cindex.xiangya <- model %>% function2(test = xiangya) %>% cindex()
  cindex.GSE108474 <- model %>% function2(test = GSE108474) %>% cindex()
  data.frame(method = i, TCGA = cindex.TCGA, CGGA = cindex.CGGA, xiangya = cindex.xiangya, GSE108474 = cindex.GSE108474)
}
res_l2 <- pbapply::pblapply(ML, FUN = function(i) {try(FUN_single_ML(i), TRUE)})
res_df2 <- do.call(rbind, res_l2)

res_df <- rbind(res_df1, res_df2)

