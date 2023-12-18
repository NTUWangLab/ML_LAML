setwd('~/machine_learning/single_algorithm/randomforest')

name_fold = 'tcga'

all_sur_data <- read.table(paste('~/machine_learning/data/',name_fold,'/LAML_patient.txt',sep=''),header=T,check.names = F)

gene_exp <- read.csv(paste('~/machine_learning/data/',name_fold,'/survival_out.csv',sep=''),header=T,check.names = F)

mixed <- merge(gene_exp, all_sur_data, by = "SampleName")

rownames(mixed) <- mixed$SampleName

mixed<- subset(mixed, select = -c(SampleName))
mixed <- as.data.frame(lapply(mixed,as.numeric))
#########################################################################

mixed_train <- mixed#[1:86,]
mixed_test <- mixed#[87:nrow(mixed),]

# mixed_train <- mixed
# mixed_test <- mixed

# 
# library(randomForest)
# library(survival)
# set.seed(123)
# 
# mixed_train.forest <- randomForest(Surv(Time,Status)~., data = mixed_train, importance = TRUE,ntree = 20000, mtry = 20)
# 
# mixed_train.forest
# 
# patient_predict <- predict(mixed_train.forest, mixed_train)
# 
# plot(mixed_train$Time, patient_predict, main = '训练集',
#      xlab = 'Survival time (days)', ylab = 'Predict')
# abline(1, 1)



library(randomForestSRC)
library(survival)

set.seed(123456)



###############
#Risk
# if (!dir.exists(name_fold)){
#   dir.create(name_fold)
# } else {
#   print("Dir already exists!")
# }
# 
# p.obj <- rfsrc(Surv(Time,Status)~.,data = mixed,
#                mtry=3,
#                nodesize=5,
#                ntree=2000
# )
# threshold <- median(p.obj$predicted)
# risk <- data.frame(riskscore = p.obj$predicted,group = risk_group <- ifelse(p.obj$predicted > threshold, "High", "Low"),Status=mixed$Status,Time=mixed$Time/365)
# 
# rt=risk[order(risk$riskscore),]   
# riskClass=rt$group
# lowLength=length(which(riskClass=="Low"))
# highLength=length(which(riskClass=="High"))
# lowMax=max(rt[which(rt$group=="Low"),]$riskscore)
# line=rt$riskscore
# line[line>1000]=1000
# pdf(file=paste(name_fold,'/risk line.pdf',sep=''),width = 8,height = 6)
# plot(line, type="p", pch=20,
#      xlab="Patients (increasing risk socre)", ylab="Risk score",
#      col=c(rep("lightblue",lowLength),rep("red",highLength)) )
# abline(h=lowMax,v=lowLength,lty=2)
# legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","lightblue"),cex=1.2)
# dev.off()
# color=as.vector(rt$Status)
# color[color==1]="red"
# color[color==0]="lightblue"
# pdf(file=paste(name_fold,'/risk point.pdf',sep=''),width = 8,height = 6)
# plot(rt$Time, pch=19,
#      xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
#      col=color)
# legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","lightblue"),cex=1.2)
# abline(v=lowLength,lty=2)
# dev.off()
# 
# library(survminer)
# rt=risk[order(risk$riskscore),]                  
# diff=survdiff(Surv(Time, Status) ~ group,data = rt)
# pValue=1-pchisq(diff$chisq,df=1)
# pValue=signif(pValue,4)
# pValue=format(pValue, scientific = TRUE)
# fit <- survfit(Surv(Time, Status) ~ group, data = rt)
# surPlot=ggsurvplot(fit, 
#                    conf.int = TRUE,
#                    data=rt,
#                    pval=paste0("p=",pValue),
#                    pval.size=5,
#                    legend.labs=c("High risk", "Low risk"),
#                    legend.title="Risk",
#                    xlab="Time(years)",
#                    break.time.by = 1,
#                    risk.table.title="",
#                    risk.table=F,
#                    risk.table.height=.25)
# pdf(file=paste(name_fold,'/survival.pdf',sep=''),onefile = FALSE,width = 5,height =4.5)
# print(surPlot)
# dev.off()
######################




v.obj <- rfsrc(Surv(Time,Status)~.,data = mixed_train,
               mtry=3,
               nodesize=5,
               ntree=2000,
               tree.err = TRUE,
               importance = TRUE
)

out.rf <- var.select(object=v.obj,conservative = "high",)

test <- predict(v.obj, mixed_test,importance = TRUE)
full <- predict(v.obj, mixed,importance = TRUE)
pdf('out_train.pdf',height = 18,width=15)
plot(v.obj)
dev.off()
pdf('out_test.pdf',height = 18,width=15)
plot(test)
dev.off()
pdf('full_matrix.pdf',height=18,width=15)
plot(full)
dev.off()

Fun <- function(x){
  return ((x - min(x))/(max(x)-min(x)))
}

RF_order <- data.frame(RF=Fun(out.rf$varselect$vimp),gene=rownames(out.rf$varselect))

save(RF_order,file='~/machine_learning/single_algorithm/RF_order.Rdata')

var_out<- data.frame(vimp=out.rf$varselect$vimp,group=rownames(out.rf$varselect))

pdf('vimp.pdf',width=10,height=5)

library(ggplot2)
library(hrbrthemes)
library(showtext)
showtext_auto()
ggplot(var_out, aes(x = reorder(group, -vimp), y = vimp)) + 
  geom_bar(stat = "identity",   
           show.legend = FALSE,   
           width = .7,fill='darkred') + aes(fill=vimp)+
  xlab("Gene") + 
  ylab("Vimp")+  theme_classic()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

dev.off()


var_out<- data.frame(depth=out.rf$varselect$depth,group=rownames(out.rf$varselect))

pdf('depth.pdf',width=20,height=10)

library(ggplot2)
library(hrbrthemes)
library(showtext)
showtext_auto()
ggplot(var_out, aes(x = reorder(group, -depth), y = depth)) + 
  geom_bar(stat = "identity",   
           show.legend = FALSE,   
           width = .7,fill='darkblue') + aes(fill=depth)+
  xlab("Gene") + 
  ylab("Depth")+  theme_classic()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

dev.off()

data <- data.frame(
  PointName = rownames(out.rf$varselect),
  XValue = out.rf$varselect$depth,
  YValue = out.rf$varselect$vimp
)

library(ggrepel)
pdf('point.pdf',height=10,width=10)
ggplot(data, aes(x = XValue, y = YValue, label = PointName)) +
  geom_point() +
  geom_text_repel(vjust = -0.5) +
  labs(x = "Depth", y = "Vimp")
dev.off()
# errRate <- c()
# for (i in 1:ncol(mixed_train)-1){ 


#   err<-mean(na.omit(v.obj$err.rate))
#   
#   errRate[i] <- err 
# }  
# 
# print(errRate)
# 
# print(which.min(errRate))
