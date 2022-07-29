library(survival)
comp <- function(model) {
  coeff <- df2 %>% dplyr::filter(Model == model) %>% dplyr::select(ENSEMBL, Coef)
  sur <- surv_expr_input %>% dplyr::select(1:2)
  expr <- surv_expr_input %>% dplyr::select(-c(1:2))
  expr <- zscore(expr) %>% as.data.frame()
  inter <- intersect(names(expr), coeff$ENSEMBL)
  expr <- expr[, inter]
  coeff <- coeff[match(inter, coeff$ENSEMBL), ]
  expr <- expr %>% mutate_all(~replace(., is.nan(.), 0))
  surv_expr <- cbind(sur, expr)
  surv_expr_train <- surv_expr %>% dplyr::select(c("time", "status", coeff$ENSEMBL))
  tmp_expr <- surv_expr_train %>% dplyr::select(coeff$ENSEMBL)
  score <- c()
  for (i in 1:nrow(tmp_expr)) {
    score[i] <- sum(coeff$Coef * tmp_expr[i, ])
  }
  temp_surv_df <- cbind(surv_expr_train[1:2], score = score %>% as.numeric())# %>% rownames_to_column("sample")
  temp_surv_df %<>% rownames_to_column("sample")
  fit <- coxph(Surv(time, status) ~ score, data = temp_surv_df)
  sum.surv<- summary(fit)
  c_index <- sum.surv$concordance[1] %>% as.numeric()
  res <- data.frame(model = model, c_index = c_index)
  return(res)
}
l <- pbapply::pblapply(unique(df2$Model), FUN = function(i) {try(comp(i), TRUE)})
df <- do.call(rbind, l)  
