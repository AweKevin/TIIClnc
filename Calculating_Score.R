ML1 <- "CoxBoost"
ML2 <- "RSF"
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


score.TCGA <- model %>% function3(test = TCGA.2) %>% rownames_to_column("sample")
score.CGGA <- model %>% function3(test = CGGA.2) %>% rownames_to_column("sample")
score.xiangya <- model %>% function3(test = xiangya.2) %>% rownames_to_column("sample")
score.GSE108474 <- model %>% function3(test = GSE108474.2) %>% rownames_to_column("sample")
