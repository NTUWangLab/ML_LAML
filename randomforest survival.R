setwd('~/machine_learning/RF_survival')

name_fold = 'tcga'

all_sur_data <- read.table(paste('~/machine_learning/data/',name_fold,'/LAML_patient.txt',sep=''),header=T,check.names = F)

gene_exp <- read.csv(paste('~/machine_learning/data/',name_fold,'/survival_out.csv',sep=''),header=T,check.names = F)

mixed <- merge(gene_exp, all_sur_data, by = "SampleName")

samplename <- mixed$SampleName

library(dplyr)

mixed<- subset(mixed, select = -c(SampleName))

mixed <- as.data.frame(lapply(mixed,as.numeric))

rownames(mixed) <- samplename

library(randomForestSRC)

library(survival)

set.seed(123456)

if (!dir.exists(name_fold)){
  dir.create(name_fold)
} else {
  print("Dir already exists!")
}

p.obj <- rfsrc(Surv(Time,Status)~.,data = mixed,
               mtry=3,
               nodesize=5,
               ntree=2000
)




threshold <- median(p.obj$predicted)

risk <- data.frame(riskscore = p.obj$predicted,group = risk_group <- ifelse(p.obj$predicted > threshold, "High", "Low"),Status=mixed$Status,Time=mixed$Time/365)

rt=risk[order(risk$riskscore),]
riskClass=rt$group
lowLength=length(which(riskClass=="Low"))
highLength=length(which(riskClass=="High"))
lowMax=max(rt[which(rt$group=="Low"),]$riskscore)
line=rt$riskscore
line[line>1000]=1000
pdf(file=paste(name_fold,'/risk line.pdf',sep=''),width = 6,height = 6)
plot(line, type="p", pch=20,
     xlab="Patients (increasing risk socre)", ylab="Risk score",
     col=c(rep("lightblue",lowLength),rep("red",highLength)) )
abline(h=lowMax,v=lowLength,lty=2)
legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","lightblue"),cex=1.2)
dev.off()
color=as.vector(rt$Status)
color[color==1]="red"
color[color==0]="lightblue"
pdf(file=paste(name_fold,'/risk point.pdf',sep=''),width = 6,height = 6)
plot(rt$Time, pch=19,
     xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","lightblue"),cex=1.2)
abline(v=lowLength,lty=2)
dev.off()

library(timeROC)
predict_full <- p.obj$predicted
full_time <- mixed$Time/365
ROC_rt=timeROC(T=full_time,delta=mixed$Status,
               marker=predict_full,cause=1,
               weighting='aalen',
               times=c(1,2,3),ROC=TRUE)
pdf(file=paste(name_fold,'/roc.pdf',sep=''),width = 6,height = 6)
plot(ROC_rt,time=1,col='green3',title=FALSE,lwd=2)
plot(ROC_rt,time=2,col='blue4',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='darkred',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=c("green3",'blue4','darkred'),lwd=2,bty = 'n')
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
pdf(file=paste(name_fold,'/survival.pdf',sep=''),onefile = FALSE,width = 6,height =6)
print(surPlot)
dev.off()