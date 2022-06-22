group <- c(rep("Immune", ncol(expr1)), rep("Cancer", ncol(expr2)))
group_list <- factor(group, levels = c("Immune", "Cancer"))

library(limma)
design <- model.matrix(~ 0 + group_list)
colnames(design) <- levels(group_list)
fit <- lmFit(exprSet, design)
contrast.matrix <- makeContrasts(Immune - Cancer,
                                 levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

tempOutput <- topTable(fit2, coef = 1, n = Inf)
nrDEG <- na.omit(tempOutput)
nrDEG <- nrDEG %>%
  rownames_to_column("symbol") %>%
  distinct()
head(nrDEG)
write.table(nrDEG, file = "nrDEG_Immune_vs_Cancer.txt", row.names = F, quote = F, sep = "\t")

nrDEG <- fread("nrDEG_Immune_vs_Cancer.txt", data.table = F)
finalDEG <- nrDEG %>% dplyr::filter(logFC > log2(1.5) & P.Value < 0.05)
dim(finalDEG)
write.table(finalDEG, file = "finalDEG_Immune_vs_Cancer.txt", row.names = F, quote = F, sep = "\t")
