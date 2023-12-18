library(pROC)
library(stringr)
library(data.table)
setwd('~/machine_learning/pROC')


gene_name = 'LSP1'

temp_exp <- fread('laml_counts.txt', sep = '\t', header = T, stringsAsFactors = FALSE)
temp_exp <- temp_exp[which(temp_exp[,c(1)]==gene_name),]

exp <- t(temp_exp)

normal_expression <- as.numeric(as.vector(exp[which(str_extract(rownames(exp),'GTEX')=='GTEX')]))
disease_expression <- as.numeric(as.vector(exp[which(str_extract(rownames(exp),'TCGA')=='TCGA')]))


all_expression <- c(normal_expression, disease_expression)
labels <- factor(c(rep(0, length(normal_expression)), rep(1, length(disease_expression))))


roc_curve <- roc(labels, all_expression)


pdf(paste('Gene/',gene_name,'.pdf',sep=''))
plot(roc_curve, main = paste("ROC Curve for Gene Expression of ",gene_name,sep=''), col = "#FFC107", lwd = 4,print.auc = TRUE, auc.polygon = TRUE,auc.polygon.col = "#FFECB3")
dev.off()