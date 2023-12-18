##设置目录
setwd('~/machine_learning/new diagram')

##载入R包
library(dplyr)
library(tibble)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(BART)
library(neuralnet)
library(xgboost)
library(Matrix)
library(xgboostExplainer)
library(ggplot2)
library(ggrepel)

##因为有些GSE数据库问题，cox的最大深度以及SPC的fold设置为3
coxphitermax = 3
SPCnfold = 3
final_result <- data.frame()

for(sample_name in c('tcga','GSE12417','GSE37642','GSE146173','GSE106291'))
{
  
  #####准备数据集
  result <- data.frame()
  all_sur_data <- read.table(paste('~/machine_learning/data/',sample_name,'/LAML_patient.txt',sep=''),header=T)
  gene_exp <- read.csv(paste('~/machine_learning/data/',sample_name,'/survival_out.csv',sep=''),header=T)
  pre_var <- colnames(gene_exp)
  pre_var <- pre_var[-1]
  mixed <- merge(gene_exp, all_sur_data, by = "SampleName")
  colnames(mixed) <- sub('Time','OS.time',colnames(mixed))
  colnames(mixed) <- sub('Status','OS',colnames(mixed))
  temp_row <- mixed$SampleName
  mixed<- mixed %>% 
    subset(select = -c(SampleName)) %>% 
    lapply(as.numeric) %>% 
    as.data.frame()
  mixed$SampleName <- temp_row
  mixed <- mixed %>% 
    select( SampleName,OS, OS.time, everything())
  seed <- 123456
  mixed <- na.omit(mixed)
  
  #70% 30%划分
  # train_idx <- sample(1:nrow(mixed), 0.7 * nrow(mixed)) # 70% 的数据作为训练集
  # mixed_test <- mixed[-train_idx, ]
  # mixed <- mixed[train_idx, ]
  
  ##所有预测变量
  #RS_COXBOOST
  #RS_STEPWISE
  #RS_RIDGE
  #RS_LASSO
  #RS_SVM
  #RS_GBDT
  #RS_SPC
  #RS_PLS
  #RS_ANN
  #RS_XGBOOST
  predict_list <- c("RS_COXBOOST","RS_STEPWISE","RS_RIDGE","RS_LASSO","RS_SVM","RS_GBDT","RS_SPC","RS_PLS","RS_ANN","RS_XGBOOST")
  comb <- combn(predict_list,2)
  ### CoxBoost ####
  set.seed(seed)
  pen <- optimCoxBoostPenalty(mixed$OS.time,mixed$OS,as.matrix(mixed[,-c(1,2,3)]),
                              trace=TRUE,start.penalty=500,parallel = T)
  cv.res <- cv.CoxBoost(mixed$OS.time,mixed$OS,as.matrix(mixed[,-c(1,2,3)]),
                        maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
  fit <- CoxBoost(mixed$OS.time,mixed$OS,as.matrix(mixed[,-c(1,2,3)]),
                  stepno=cv.res$optimal.step,penalty=pen$penalty)
  RS_COXBOOST <- data.frame(RS = as.numeric(predict(fit,newdata=mixed[,-c(1,2,3)],newtime=mixed[,3], newstatus=mixed[,2], type="lp")),name='CoxBoost')
  rs <- cbind(mixed[,2:3],RS = as.numeric(predict(fit,newdata=mixed[,-c(1,2,3)],newtime=mixed[,3], newstatus=mixed[,2], type="lp")))
  cc <- t(data.frame('CoxBoost'=summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]))
  result <- rbind(result,cc)
  rm(list=c('fit'))
  
  
  #### Stepwise Cox ####
  for (direction in c("both", "backward", "forward")) {coxphitermax
    fit <- step(coxph(Surv(OS.time,OS)~.,mixed[,-c(1)],iter.max=coxphitermax),direction = direction)
    if(direction == 'forward') {RS_STEPWISE <- data.frame(RS=as.numeric(predict(fit,type = 'risk',newdata = mixed[,-c(1)])),name='Stepwise Cox')}
    rs <- cbind(mixed[,c(2,3)],RS = as.numeric(predict(fit,type = 'risk',newdata = mixed[,-c(1)])))
    cc <- data.frame(name = summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
    colnames(cc) <- c(paste('Stepwise Cox','[',direction,']',sep=''))
    cc <- t(cc)
    result <- rbind(result,cc)
    rm(list=c('fit'))
  }
  
  
  #### Lasso,Ridge,Enet ####
  for (alpha in seq(0,1,0.1)) {
    set.seed(seed)
    fit <- cv.glmnet(as.matrix(mixed[,-c(1,2,3)]), as.matrix(Surv(mixed$OS.time,mixed$OS)),family = "cox",alpha=alpha,nfolds = 10)
    rs <- cbind(mixed[,2:3],RS = as.numeric(predict(fit,type='link',newx=as.matrix(mixed[,-c(1,2,3)]),s=fit$lambda.min)))
    cc <- data.frame(name = as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]))
    colnames(cc) <- c(paste0('Enet','[a=',alpha,']',sep=''))
    if(alpha == 0)
    {
      colnames(cc) <- c(paste0('Ridge',sep=''))
      RS_RIDGE <- data.frame(RS=as.numeric(predict(fit,type='link',newx=as.matrix(mixed[,-c(1,2,3)]),s=fit$lambda.min)),name='Ridge')
    }
    if(alpha == 1) 
    {
      colnames(cc) <- c(paste0('Lasso',sep=''))
      RS_LASSO <- data.frame(RS=as.numeric(predict(fit,type='link',newx=as.matrix(mixed[,-c(1,2,3)]),s=fit$lambda.min)),name='Lasso')  
    }
    cc <- t(cc)
    colnames(cc) <- c('C')
    result <- rbind(result,cc)
    rm(list=c('fit'))
  }
  
  
  #### Survival SVM ####
  fit <- survivalsvm(Surv(OS.time,OS)~., data= mixed[,-c(1)], gamma.mu = 1)
  RS_SVM <- data.frame(RS=as.numeric(predict(fit, mixed[,-c(1)])$predicted),name='Survival SVM')
  rs <- cbind(mixed[,2:3],RS=as.numeric(predict(fit, mixed[,-c(1)])$predicted))
  as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  cc <- data.frame(name = as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]))
  colnames(cc) <- c('Survival SVM')
  cc <- t(cc)
  colnames(cc) <- c('C')
  result <- rbind(result,cc)
  rm(list=c('fit'))
  
  
  #### GBDT ####
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time,OS)~.,data = mixed[,-c(1)],distribution = 'coxph',
             n.trees = 10000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 6)
  best <- which.min(fit$cv.error)
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time,OS)~.,data = mixed[,-c(1)],distribution = 'coxph',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 8)
  RS_GBDT <- data.frame(RS = as.numeric(predict(fit,mixed[,-c(1)],n.trees = best,type = 'link')),name='GBDT')
  rs <- cbind(mixed[,2:3],RS=as.numeric(predict(fit,mixed[,-c(1)],n.trees = best,type = 'link')))
  cc <- data.frame(name = as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]))
  colnames(cc) <- c('GBDT')
  cc <- t(cc)
  colnames(cc) <- c('C')
  result <- rbind(result,cc)
  rm(list=c('fit'))
  
  
  #### Supervised principal components ####
  data <- list(x=t(mixed[,-c(1,2,3)]),y=mixed$OS.time,censoring.status=mixed$OS,featurenames=colnames(mixed)[-c(1,2,3)])
  set.seed(seed)
  fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
  set.seed(seed)
  cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default
                       n.fold = SPCnfold,
                       n.components=3,
                       min.features=5,
                       max.features=nrow(data$x),
                       compute.fullcv= TRUE,
                       compute.preval=TRUE)
  RS_SPC <- data.frame(RS = as.numeric(superpc.predict(fit,data,list(x=t(mixed[,-c(1,2,3)]),y=mixed$OS.time,censoring.status=mixed$OS,featurenames=colnames(mixed)[-c(1,2,3)]),threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)$v.pred),name = "Supervised principal components")
  rs <- cbind(mixed[,2:3],RS=as.numeric(superpc.predict(fit,data,list(x=t(mixed[,-c(1,2,3)]),y=mixed$OS.time,censoring.status=mixed$OS,featurenames=colnames(mixed)[-c(1,2,3)]),threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)$v.pred))
  cc <- data.frame(name=as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]))
  colnames(cc) <- c('Survival PCA')
  cc <- t(cc)
  colnames(cc) <- c('C')
  result <- rbind(result,cc)
  rm(list=c('fit'))
  
  
  #### plsRcox ####
  set.seed(seed)
  
  pdf('plsRcox.pdf')
  cv.plsRcox.res=cv.plsRcox(list(x=mixed[,-c(1,2,3)],time=mixed$OS.time,status=mixed$OS),nt=10,verbose = FALSE)
  fit <- plsRcox(mixed[,-c(1,2,3)],time=mixed$OS.time,event=mixed$OS,nt=as.numeric(cv.plsRcox.res[5]))
  RS_PLS <- data.frame(RS=as.numeric(predict(fit,type="lp",newdata=mixed[,-c(1,2,3)])),name='plsRcox')
  rs <- cbind(mixed[,2:3],RS=as.numeric(predict(fit,type="lp",newdata=mixed[,-c(1,2,3)])))
  cc <- data.frame(namx=as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]))
  colnames(cc) <- c('plsRcox')
  cc <- t(cc)
  colnames(cc) <- c('C')
  result <- rbind(result,cc)
  rm(list=c('fit'))
  dev.off()
  #ANN
  set.seed(seed)
  for (max in seq(10,80,10)) {
    fit <- neuralnet(OS+OS.time ~ ., data=mixed[,-c(1)], hidden=max)
    while(length(fit)<14)
    {
      fit <- neuralnet(OS+OS.time ~ ., data=mixed[,-c(1)], hidden=max)
    }
    if(max==40) {
      RS_ANN<-data.frame(RS=predict(fit,newdata = mixed[,-c(1)]),name='ANN')
      colnames(RS_ANN) <- c('RS','RS.1','name')
    }
    rs <- cbind(mixed[,2:3],RS=predict(fit,newdata = mixed[,-c(1)]))
    cc <- data.frame(name=as.numeric(summary(coxph(Surv(OS.time,OS)~RS.1+RS.2,rs))$concordance[1]))
    colnames(cc) <- c(paste('ANN [hidden=',max,']',sep=''))
    cc <- t(cc)
    colnames(cc) <- c('C')
    result <- rbind(result,cc)
    rm(list=c('fit'))
  }
  
  #XGboost
  set.seed(seed)
  dtrain <- xgb.DMatrix(as.matrix(mixed[,-c(1,2,3)]), label = mixed$OS.time,weight = mixed$OS)
  for (max in seq(1,10,1)) {
    params <- list(
      objective = "survival:cox",
      eval_metric = "cox-nloglik",
      eta = 0.01,
      max_depth = max,
      subsample = 0.8,
      colsample_bytree = 0.8
    )
    xgb_model <- xgb.train(params = params,nrounds = 100,dtrain)
    if(max==10) {RS_XGBOOST<-data.frame(RS=predict(xgb_model,dtrain),name='XGboot')}
    rs <- cbind(mixed[,2:3],RS=predict(xgb_model,dtrain))
    cc <- data.frame(name=as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]))
    colnames(cc) <- c(paste('XGboost [max_depth=',max,']',sep=''))
    cc <- t(cc)
    colnames(cc) <- c('C')
    result <- rbind(result,cc)
  }
  
  for(i in 1:ncol(comb))
  {
    eval(parse(text = paste('rs <- cbind(mixed[,2:3],RS1=',comb[1,i],'$RS,RS2=',comb[2,i],'$RS)',sep='')))
    cc <- data.frame(name=as.numeric(summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1]))
    eval(parse(text = paste('colnames(cc) <- c(paste(',comb[1,i],'$name[1],',comb[2,i],'$name[1],sep=" + "))',sep='')))
    cc <- t(cc)
    colnames(cc) <- c('C')
    result <- rbind(result,cc)
  }
  
  #### Random survival forest ####
  set.seed(seed)
  fit_rf <- rfsrc(Surv(OS.time,OS)~.,data = mixed[,-c(1)],
                  ntree = 1000,nodesize = 5,
                  splitrule = 'logrank',
                  importance = T,
                  proximity = T,
                  forest = T,
                  seed = seed)
  rs <- cbind(mixed[,2:3],RS=as.numeric(predict(fit_rf,newdata = mixed[,-c(1)])$predicted))
  cc <- data.frame(name=as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]))
  colnames(cc) <- c('Survival RF')
  cc <- t(cc)
  colnames(cc) <- c('C')
  result <- rbind(result,cc)
  
  #ANN + Random survival forest
  set.seed(seed)
  for (max in seq(10,80,10)) {
    fit <- neuralnet(OS+OS.time ~ ., data=mixed[,-c(1)], hidden=max)
    while(length(fit)<14)
    {
      fit <- neuralnet(OS+OS.time ~ ., data=mixed[,-c(1)], hidden=max)
    }
    rs <- cbind(mixed[,2:3],RS1=as.numeric(predict(fit_rf,newdata = mixed[,-c(1)])$predicted),RS2=predict(fit,newdata = mixed[,-c(1)]))
    cc <- data.frame(name=as.numeric(summary(coxph(Surv(OS.time,OS)~RS1+RS2.1+RS2.2,rs))$concordance[1]))
    colnames(cc) <- c(paste('Survival RF + ANN [hidden=',max,']',sep=''))
    cc <- t(cc)
    colnames(cc) <- c('C')
    result <- rbind(result,cc)
    rm(list=c('fit'))
  }
  
  
  #GBDT + Random survival forest
  set.seed(seed)
  fit_GB <- gbm(formula = Surv(OS.time,OS)~.,data = mixed[,-c(1)],distribution = 'coxph',
                n.trees = 10000,
                interaction.depth = 3,
                n.minobsinnode = 10,
                shrinkage = 0.001,
                cv.folds = 10,n.cores = 6)
  best <- which.min(fit_GB$cv.error)
  set.seed(seed)
  fit_GB <- gbm(formula = Surv(OS.time,OS)~.,data = mixed[,-c(1)],distribution = 'coxph',
                n.trees = best,
                interaction.depth = 3,
                n.minobsinnode = 10,
                shrinkage = 0.001,
                cv.folds = 10,n.cores = 8)
  rs <- cbind(mixed[,2:3],RS1=as.numeric(predict(fit_rf,newdata = mixed[,-c(1)])$predicted),RS2=as.numeric(predict(fit_GB,mixed[,-c(1)],n.trees = best,type = 'link')))
  cc <- data.frame(name = as.numeric(summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1]))
  colnames(cc) <- c('Survival RF + GBDT')
  cc <- t(cc)
  colnames(cc) <- c('C')
  result <- rbind(result,cc)
  rm(list=c('fit_GB'))
  
  
  
  
  #Stepwise Cox + Random survival forest
  for (direction in c("both", "backward", "forward")) {
    fit <- step(coxph(Surv(OS.time,OS)~.,mixed[,-c(1)],iter.max=coxphitermax),direction = direction)
    rs = cbind(mixed[,c(2,3)],RS1=as.numeric(predict(fit_rf,newdata = mixed[,-c(1)])$predicted),RS2=as.numeric(predict(fit,type = 'risk',newdata = mixed[,-c(1)])))
    cc <- data.frame(name = summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1])
    colnames(cc) <- c(paste('Survival RF + Stepwise Cox',' [',direction,']',sep=''))
    cc <- t(cc)
    result <- rbind(result,cc)
    rm(list=c('fit'))
  }
  
  #### Supervised principal components ####
  data <- list(x=t(mixed[,-c(1,2,3)]),y=mixed$OS.time,censoring.status=mixed$OS,featurenames=colnames(mixed)[-c(1,2,3)])
  set.seed(seed)
  fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
  set.seed(seed)
  cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default
                       n.fold = SPCnfold,
                       n.components=3,
                       min.features=5,
                       max.features=nrow(data$x),
                       compute.fullcv= TRUE,
                       compute.preval=TRUE)
  RS_SPC <- data.frame(RS = as.numeric(superpc.predict(fit,data,list(x=t(mixed[,-c(1,2,3)]),y=mixed$OS.time,censoring.status=mixed$OS,featurenames=colnames(mixed)[-c(1,2,3)]),threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)$v.pred),name = "Supervised principal components")
  rs <- cbind(mixed[,2:3],RS1=as.numeric(predict(fit_rf,newdata = mixed[,-c(1)])$predicted),RS2=as.numeric(superpc.predict(fit,data,list(x=t(mixed[,-c(1,2,3)]),y=mixed$OS.time,censoring.status=mixed$OS,featurenames=colnames(mixed)[-c(1,2,3)]),threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)$v.pred))
  cc <- data.frame(name=as.numeric(summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1]))
  colnames(cc) <- c('Survival RF + Survival PCA')
  cc <- t(cc)
  colnames(cc) <- c('C')
  result <- rbind(result,cc)
  rm(list=c('fit'))
  
  #plsRcox + Random survival forest
  set.seed(seed)
  pdf('pls.pdf')
  cv.plsRcox.res=cv.plsRcox(list(x=mixed[,-c(1,2,3)],time=mixed$OS.time,status=mixed$OS),nt=10,verbose = FALSE)
  fit <- plsRcox(mixed[,-c(1,2,3)],time=mixed$OS.time,event=mixed$OS,nt=as.numeric(cv.plsRcox.res[5]))
  
  rs <- cbind(mixed[,2:3],RS1=as.numeric(predict(fit_rf,newdata = mixed[,-c(1)])$predicted),RS2=as.numeric(predict(fit,type="lp",newdata=mixed[,-c(1,2,3)])))
  
  cc <- data.frame(namx=as.numeric(summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1]))
  colnames(cc) <- c('Survival RF + plsRcox')
  cc <- t(cc)
  colnames(cc) <- c('C')
  result <- rbind(result,cc)
  rm(list=c('fit'))
  dev.off()
  
  #XGboost + Random survival forest
  set.seed(seed)
  dtrain <- xgb.DMatrix(as.matrix(mixed[,-c(1,2,3)]), label = mixed$OS.time,weight = mixed$OS)
  for (max in seq(1,5,1)) {
    params <- list(
      objective = "survival:cox",
      eval_metric = "cox-nloglik",
      eta = 0.01,
      max_depth = max,
      subsample = 0.8,
      colsample_bytree = 0.8
    )
    xgb_model <- xgb.train(params = params,nrounds = 100,dtrain)
    rs <- cbind(mixed[,2:3],RS1=as.numeric(predict(fit_rf,newdata = mixed[,-c(1)])$predicted),RS2=predict(xgb_model,dtrain))
    cc <- data.frame(name=as.numeric(summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1]))
    colnames(cc) <- c(paste('Survival RF + XGboost [max_depth=',max,']',sep=''))
    cc <- t(cc)
    colnames(cc) <- c('C')
    result <- rbind(result,cc)
  }
  
  
  #Lasso + Random survival forest
  for (alpha in seq(0,1,0.1)) {
    set.seed(seed)
    fit <- cv.glmnet(as.matrix(mixed[,-c(1,2,3)]), as.matrix(Surv(mixed$OS.time,mixed$OS)),family = "cox",alpha=alpha,nfolds = 10)
    rs <- cbind(mixed[,2:3],RS1=as.numeric(predict(fit_rf,newdata = mixed[,-c(1)])$predicted),RS2 = as.numeric(predict(fit,type='link',newx=as.matrix(mixed[,-c(1,2,3)]),s=fit$lambda.min)))
    cc <- data.frame(name = as.numeric(summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1]))
    colnames(cc) <- c(paste0('Survival RF + Enet',' [a=',alpha,']',sep=''))
    if(alpha == 0) {colnames(cc) <- c(paste0('Survival RF + Ridge',sep=''))}
    if(alpha == 1) {colnames(cc) <- c(paste0('Survival RF + Lasso',sep=''))}
    cc <- t(cc)
    colnames(cc) <- c('C')
    result <- rbind(result,cc)
    rm(list=c('fit'))
  }
  
  
  
  colnames(result) = c(sample_name)
  if(sample_name=='tcga')
  {
    final_result <- result
  }else{
    final_result <- cbind(final_result,result)
  }
}

save(final_result,file='result_final.Rdata')


final_result$Total <- rowSums(final_result)/5

final_result <- final_result[order(final_result$Total),]

data <- data.frame(
  x = rep(1:5, each = length(final_result$tcga)),
  y = rep(1:length(final_result$tcga), times = 5),
  value = c(round(final_result$tcga, 2),
            round(final_result$GSE12417, 2),
            round(final_result$GSE37642, 2),
            round(final_result$GSE146173, 2),
            round(final_result$GSE106291, 2)),
  text = rownames(final_result),
  sample = rep(colnames(final_result)[1:5]),
  fill_column = rep(c("no_fill", "no_fill", "no_fill","no_fill", "fill"), each = length(final_result$tcga))
)

# 创建颜色标尺
#color_scale <- scale_fill_gradient(low = "#FFCDD2", high = "#B71C1C")

color_scale <- scale_fill_gradient(low = "white", high = "#B71C1C",limits = c(0.4, 1))

pdf('Q.pdf',width=20,height=22)

ggplot(data, aes(x = x, y = y, fill = value)) +
  geom_tile(width = 0.8, height = 0.9) + 
  geom_line(aes(x=1,y=1))+
  theme_void()+geom_text(aes(label = value), color = "black", size = 4)+
  theme_void()+geom_text(aes(label = text,x=x+0.01,y=y),data = subset(data, fill_column == "fill"), color = "black", size = 4)+
  #theme_void()+geom_text(aes(label = sample,x=x,y=y+0.01),data = subset(data,fill_row), color = "black", size = 6)+
  scale_x_continuous(expand = c(4, 0))+
  color_scale

dev.off()

pdf('C.pdf',width=40,height=10)

ggplot(final_result, aes(x = reorder(rownames(final_result), -Total), y = Total)) + 
  geom_bar(stat = "identity",   
           show.legend = FALSE,   
           width = 0.8) + aes(fill=Total)+theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=20),axis.text.y = element_text( hjust = 1,size=20))

dev.off()
