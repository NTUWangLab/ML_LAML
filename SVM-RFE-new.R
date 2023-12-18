setwd('~/machine_learning/single_algorithm/SVM')
library(dplyr)
all_sur_data <- read.table('~/machine_learning/data/tcga/LAML_patient.txt',header=T,check.names = F)

gene_exp <- read.csv('~/machine_learning/data/tcga/survival_out.csv',header=T,check.names = F)

mixed <- merge(gene_exp, all_sur_data, by = "SampleName")
rownames(mixed) <- mixed$SampleName
mixed <- na.omit(mixed)
gene_list <- colnames(mixed)

mixed <- as.data.frame(apply(mixed,2,function(x) as.numeric(as.character(x))))
colnames(mixed) <- gene_list
mixed <- cbind(data.frame(status = mixed$Status),mixed)
mixed <- select(mixed,-c(SampleName,Time,Status))

library(e1071)

source("~/machine_learning/single_algorithm/SVM/msvmRFE.R")

nfold = 10 #10倍交叉验证
nrows = nrow(mixed)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))


results = lapply(folds, svmRFE.wrap, mixed, k=10, halve.above=100)
top.features = WriteFeatures(results, mixed, save=F)
featsweep = lapply(1:5, FeatSweep.wrap, results, mixed)

Fun <- function(x){
  return ((x - min(x))/(max(x)-min(x)))
}

SVM_order <- data.frame(SVM=Fun(length(top.features$AvgRank)-top.features$AvgRank),gene=top.features$FeatureName)
save(SVM_order,file='~/machine_learning/single_algorithm/SVM_order.Rdata')

no.info = min(prop.table(table(mixed[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

pdf("svm_rfe.pdf", height = 8, width = 10)
plot(errors,type=c("b"),col = 4, lty = 2)
text(x = errors,y = sprintf("%.5f", errors),label  = sprintf("%.5f", errors), pos = 3, offset = 0.5, col = "red")
dev.off()


library(ggplot2)
pdf("rank.pdf",width=10)
p <- ggplot(top.features, aes(x = reorder(FeatureName, AvgRank), y = AvgRank)) +
  geom_bar(stat = "identity",   
           show.legend = FALSE,   
           width = .8,fill = 'lightblue') +
  labs(x = "FeatureName", y = "rank") + theme_classic()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
dev.off()
