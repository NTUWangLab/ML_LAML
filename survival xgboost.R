library(survival)
library(xgboost)
library(Matrix)
library(xgboostExplainer)
library(timeROC)
library(dplyr)
library(survminer)
setwd('~/machine_learning/XGboost')

name_fold = 'tcga'

all_sur_data <- read.table(paste('~/machine_learning/data/',name_fold,'/LAML_patient.txt',sep=''),header=T,check.names = F)

gene_exp <- read.csv(paste('~/machine_learning/data/',name_fold,'/survival_out.csv',sep=''),header=T,check.names = F)

mixed <- merge(gene_exp, all_sur_data, by = "SampleName")
rownames(mixed) <- mixed$SampleName
mixed <- na.omit(as.matrix(mixed))

gene_list <- colnames(mixed)

mixed <- as.data.frame(apply(mixed,2,function(x) as.numeric(as.character(x))))
colnames(mixed) <- gene_list
mixed <- cbind(data.frame(status = mixed$Status),mixed)
mixed <- select(mixed,-c(SampleName))

mixed_train <- mixed[ 1:86, ]
mixed_test <- mixed[ 87:nrow(mixed), ]
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
  max_depth = 5,
  subsample = 0.8,
  colsample_bytree = 0.8
)


set.seed(123)
xgb_model <- xgb.train(params = params,nrounds = 100,dfull)



####KM curve

if (!dir.exists(name_fold)){
  dir.create(name_fold)
} else {
  print("Dir already exists!")
}

threshold <- median(predict(xgb_model,dfull))
risk <- data.frame(riskscore = predict(xgb_model,dfull),group = risk_group <- ifelse(predict(xgb_model,dfull) > threshold, "High", "Low"),Status=mixed$Status,Time=mixed$Time/365)

rt=risk[order(risk$riskscore),]
riskClass=rt$group
lowLength=length(which(riskClass=="Low"))
highLength=length(which(riskClass=="High"))
lowMax=max(rt[which(rt$group=="Low"),]$riskscore)
line=rt$riskscore
line[line>1000]=1000
pdf(file=paste(name_fold,'/risk line.pdf',sep=''),width = 8,height = 6)
plot(line, type="p", pch=20,
     xlab="Patients (increasing risk socre)", ylab="Risk score",
     col=c(rep("lightblue",lowLength),rep("red",highLength)) )
abline(h=lowMax,v=lowLength,lty=2)
legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","lightblue"),cex=1.2)
dev.off()
color=as.vector(rt$Status)
color[color==1]="red"
color[color==0]="lightblue"
pdf(file=paste(name_fold,'/risk point.pdf',sep=''),width = 8,height = 6)
plot(rt$Time, pch=19,
     xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","lightblue"),cex=1.2)
abline(v=lowLength,lty=2)
dev.off()

library(survminer)
rt=risk[order(risk$riskscore),]
diff=survdiff(Surv(Time, Status) ~ group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(Time, Status) ~ group, data = rt)
surPlot=ggsurvplot(fit,
                   conf.int = TRUE,
                   data=rt,
                   pval=paste0("p=",pValue),
                   pval.size=5,
                   legend.labs=c("High risk", "Low risk"),
                   legend.title="Risk",
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   risk.table=F,
                   risk.table.height=.25)
pdf(file=paste(name_fold,'/survival.pdf',sep=''),onefile = FALSE,width = 5,height =4.5)
print(surPlot)
dev.off()


######roc
pre_xgb = round(predict(xgb_model,newdata = dfull))
table(full_status,pre_xgb,dnn=c("true","pre"))
xgboost_roc <- roc(full_status,as.numeric(pre_xgb))
pdf(paste(name_fold,'/model_roc.pdf',sep=''))
plot(xgboost_roc, print.auc=TRUE, auc.polygon=TRUE,
     grid=c(0.1, 0.2),grid.col=c("green", "red"),
     max.auc.polygon=TRUE,auc.polygon.col="skyblue",
     print.thres=TRUE,main='ROC curve')
dev.off()

######time_roc
predict_test <- predict(xgb_model,dtest)
test_time <- test_time/365
ROC_rt=timeROC(T=test_time,delta=test_status,
               marker=predict_test,cause=1,
               weighting='aalen',
               times=c(1,2,3),ROC=TRUE)
pdf(file=paste(name_fold,'/roc_test.pdf',sep=''),width=5,height=5)
plot(ROC_rt,time=1,col='green3',title=FALSE,lwd=2)
plot(ROC_rt,time=2,col='blue4',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='darkred',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=c("green3",'blue4','darkred'),lwd=2,bty = 'n')
dev.off()

predict_train <- predict(xgb_model,dtrain)
train_time <- train_time/365
ROC_rt=timeROC(T=train_time,delta=train_status,
               marker=predict_train,cause=1,
               weighting='aalen',
               times=c(1,2,3),ROC=TRUE)
pdf(file=paste(name_fold,'/roc_train.pdf',sep=''),width=5,height=5)
plot(ROC_rt,time=1,col='green3',title=FALSE,lwd=2)
plot(ROC_rt,time=2,col='blue4',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='darkred',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=c("green3",'blue4','darkred'),lwd=2,bty = 'n')
dev.off()

predict_full <- predict(xgb_model,dfull,type = "risk")
full_time <- full_time/365
ROC_rt=timeROC(T=full_time,delta=full_status,
               marker=predict_full,cause=1,
               weighting='aalen',
               times=c(1,2,3),ROC=TRUE)
pdf(file=paste(name_fold,'/roc_full.pdf',sep=''),width=5,height=5)
plot(ROC_rt,time=1,col='green3',title=FALSE,lwd=2)
plot(ROC_rt,time=2,col='blue4',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='darkred',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=c("green3",'blue4','darkred'),lwd=2,bty = 'n')
dev.off()




######importance
# feature_importance <- xgb.importance(model = xgb_model)
# pdf("rank.pdf",width=10)
# p <- ggplot(feature_importance, aes(x = reorder(feature_importance$Feature, -feature_importance$Gain), y = feature_importance$Gain)) +
#   geom_bar(stat = "identity",   
#            show.legend = FALSE,   
#            width = .8,fill = 'orange2') +
#   labs(x = "FeatureName", y = "Gain") +
#   theme(axis.text.x = element_text(angle = 60, hjust = 1))
# p
# dev.off()