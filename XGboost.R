library(survival)
library(xgboost)
library(Matrix)
library(xgboostExplainer)
library(timeROC)
library(dplyr)
library(survminer)
setwd('~/machine_learning/single_algorithm/XGboost')
set.seed(2)
all_sur_data <- read.table('~/machine_learning/data/tcga/LAML_patient.txt',header=T,check.names = F)

gene_exp <- read.csv('~/machine_learning/data/tcga/survival_out.csv',header=T,check.names = F)

mixed <- merge(gene_exp, all_sur_data, by = "SampleName")
rownames(mixed) <- mixed$SampleName
mixed <- na.omit(mixed)
gene_list <- colnames(mixed)

mixed <- as.data.frame(apply(mixed,2,function(x) as.numeric(as.character(x))))
colnames(mixed) <- gene_list
mixed <- cbind(data.frame(status = mixed$Status),mixed)
mixed <- select(mixed,-c(SampleName))

mixed_train <- mixed#[ 1:86, ]
mixed_test <- mixed#[ 87:nrow(mixed), ]
mixed_full <- mixed

train_status <- mixed_train$Status
test_status <- mixed_test$Status
full_status <- mixed_full$Status

train_time <- mixed_train$Time
test_time <- mixed_test$Time
full_time <- mixed_full$Time

mixed_train <- select(mixed_train,-c(Status,Time))
mixed_test <- select(mixed_test,-c(Status,Time))
mixed_full <- select(mixed_full,-c(Status,Time))
#xgb_data <- xgb.DMatrix(data = as.matrix(features), label = as.matrix(status))

dtrain <- xgb.DMatrix(as.matrix(mixed_train), label = train_time,weight = train_status)
dtest <- xgb.DMatrix(as.matrix(mixed_test), label = test_time, weight = test_status)
dfull <- xgb.DMatrix(as.matrix(mixed_full), label = full_time, weight = full_status)


params <- list(
  objective = "survival:cox",
  eval_metric = "cox-nloglik",
  eta = 0.01,
  max_depth = 3,
  subsample = 0.8,
  colsample_bytree = 0.8
)


set.seed(123)
xgb_model <- xgb.train(params = params,nrounds = 100,dtrain)

# pre_xgb = round(predict(xgb_model,newdata = dtest))
# table(test_status,pre_xgb,dnn=c("true","pre"))
# xgboost_roc <- roc(test_status,as.numeric(pre_xgb))
# plot(xgboost_roc, print.auc=TRUE, auc.polygon=TRUE, 
#      grid=c(0.1, 0.2),grid.col=c("green", "red"), 
#      max.auc.polygon=TRUE,auc.polygon.col="skyblue", 
#      print.thres=TRUE,main='ROC curve')

feature_importance <- xgb.importance(model = xgb_model)

predict_test <- predict(xgb_model,dtest)
test_time <- test_time/365
ROC_rt=timeROC(T=test_time,delta=test_status,
               marker=predict_test,cause=1,
               weighting='aalen',
               times=c(1,2,3),ROC=TRUE)
pdf(file='roc_test.pdf',width=5,height=5)
plot(ROC_rt,time=1,col='green3',title=FALSE,lwd=2)
plot(ROC_rt,time=2,col='blue4',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='darkred',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=c("green3",'blue4','darkred'),lwd=2,bty = 'n')
dev.off()

predict_full <- predict(xgb_model,dfull)
full_time <- full_time/365
ROC_rt=timeROC(T=full_time,delta=full_status,
               marker=predict_full,cause=1,
               weighting='aalen',
               times=c(1,2,3),ROC=TRUE)
pdf(file='roc_full.pdf',width=5,height=5)
plot(ROC_rt,time=1,col='green3',title=FALSE,lwd=2)
plot(ROC_rt,time=2,col='blue4',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='darkred',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=c("green3",'blue4','darkred'),lwd=2,bty = 'n')
dev.off()

Fun <- function(x){
  return ((x - min(x))/(max(x)-min(x)))
}

feature_importance <- feature_importance[order(feature_importance$Gain),]
XGBOOST_order <- data.frame(XGBOOST=Fun(feature_importance$Gain),gene=feature_importance$Feature)
save(XGBOOST_order,file='~/machine_learning/single_algorithm/XGBOOST_order.Rdata')
print(colnames(mixed_full)[which(!colnames(mixed_full)%in%XGBOOST_order$gene)])
pdf("rank.pdf",width=10)
p <- ggplot(feature_importance, aes(x = reorder(feature_importance$Feature, -feature_importance$Gain), y = feature_importance$Gain)) +
  geom_bar(stat = "identity",   
           show.legend = FALSE,   
           width = .8,fill = 'orange2') +
  labs(x = "FeatureName", y = "Gain") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
dev.off()