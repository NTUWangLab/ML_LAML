
library(neuralnet)
library(dplyr)

sample_name <- 'tcga'
setwd('~/machine_learning/single_algorithm/ANN')
all_sur_data <- read.table(paste('~/machine_learning/data/',sample_name,'/LAML_patient.txt',sep=''),header=T)
gene_exp <- read.csv(paste('~/machine_learning/data/',sample_name,'/survival_out.csv',sep=''),header=T,check.names = F)
mixed <- merge(gene_exp, all_sur_data, by = "SampleName")
rownames(mixed) <- mixed$SampleName
mixed<- subset(mixed, select = -c(SampleName))
data <- na.omit(as.data.frame(lapply(mixed,as.numeric)))




set.seed(123) 
train_idx <- sample(1:nrow(data), 0.8 * nrow(data))
train_data <- data[train_idx, ]
test_data <- data[-train_idx, ]

set.seed(123)
#100%
nn <- neuralnet(Status+Time ~ ., data=data, hidden=15)
if(length(nn)<11) {nn <- neuralnet(Status+Time ~ ., data=data, hidden=15)}
#%80
#nn <- neuralnet(target ~ ., data=train_data, hidden=9)
summary(nn)

# train_pred <- compute(nn, train_data[, -ncol(train_data)])$net.result
# test_pred <- compute(nn, test_data[, -ncol(test_data)])$net.result
# 
# train_accuracy <- sum(train_data$target == round(train_pred)) / nrow(train_data)
# test_accuracy <- sum(test_data$target == round(test_pred)) / nrow(test_data)

# cat("Train Accuracy:", train_accuracy, "\n")
# cat("Test Accuracy:", test_accuracy, "\n")
weight <- data.frame(sum = rowSums(abs(nn$weights[[1]][[1]]))[-1],row.names = colnames(gene_exp)[-1])

#排名
# ANN_order <- data.frame(ANN=rep(length(weight[,1]):1),gene=rownames(weight)[order(weight$sum)])
#权重
Fun <- function(x){
  return ((x - min(x))/(max(x)-min(x)))
}

ANN_order <- data.frame(ANN=Fun(weight[order(weight$sum),]),gene=rownames(weight)[order(weight$sum)])
save(ANN_order,file='~/machine_learning/single_algorithm/ANN_order.Rdata')

library(ggplot2)
pdf("rank.pdf",width=12)
p <- ggplot(weight, aes(x = reorder(rownames(weight), -sum), y = sum)) +
  geom_bar(stat = "identity",   
           show.legend = FALSE,   
           width = .8,fill = 'darkblue') +
  labs(x = "FeatureName", y = "weight") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
dev.off()
