
name_fold = 'lasso'
tumor_id = 'ML'
main_fold = 'machine_learning'
setwd(paste('~/',main_fold,'/',name_fold,sep=''))
all_sur_data <- read.table(paste('~/history_lasso/m5c_laml/tcga/LAML_sur.txt',sep=''),header=T,check.names=F)

gene_exp <- read.csv(paste('~/machine_learning/data/tcga/survival_out.csv',sep=''),header=T,check.names=F)

rownames(gene_exp) = gene_exp$SampleName

erged_df <- merge(gene_exp, all_sur_data, by = "SampleName")

erged_df$Time <- erged_df$Time/365
#erged_df <- erged_df[-85,]
erged_df <- na.omit(erged_df)

survival_time <- erged_df$Time
survival_status <- erged_df$Status
patient_id <- erged_df$SampleName
erged_df <-subset(erged_df, select = -Status)
erged_df <-subset(erged_df, select = -SampleName)
erged_df <-subset(erged_df, select = -Time)
erged_df=as.data.frame(lapply(erged_df,as.numeric))
rownames(erged_df) <- patient_id
standardize <- function(x) {
  rowmean <- apply(x, 1, mean)
  rowsd <- apply(x, 1, sd)  
  rv <- sweep(x, 1, rowmean,"-")  #表达量-均值
  rv <- sweep(rv, 1, rowsd, "/")  #再除以标准差
  return(rv)
}

#########################################

library(survival)
library(glmnet)
seed <- sample(1:10000,size=1)

set.seed(seed)
surv_obj <- Surv(event = survival_status,time = survival_time)

fit <- glmnet(as.matrix(erged_df), surv_obj, family = "cox",maxit=1000)

cvfit <- cv.glmnet(as.matrix(erged_df), surv_obj, family = "cox",maxit=1000)


pdf("lambda.pdf",width = 8,height = 6)
plot(fit, xvar = "lambda", label = TRUE)
dev.off()


pdf("fit.pdf",width = 8,height = 6)
plot(cvfit)
dev.off()


erged_df$Time <- survival_time
erged_df$Status <- survival_status

# coef=coef(fit, s = cvfit$lambda.min)
# for(j in 1:length(coef))
# {
#   coef[j,] = 0
# }
#m5c_laml
#coef[c('NOP2','NSUN2','DNMT3A','TET2','TET3'),] = c(-0.420828468727705,0.184971417695989,-0.445835537999184,-0.0835039701653756,-0.0612006001382076)
#mac-m6a

#coef[c('DUSP2','HMGB2','JUNB','NR4A2','NR4A3','PIM3'),] = c(-0.0078304779,0.0009869564,0.0002102549,0.0006363361,-0.0028197840,0.0004534177)
#原始
#coef[c('DUSP2','HMGB2','JUNB','NR4A2','NR4A3'),] = c(-0.0078304779,0.0009869564,0.0002102549,0.0006363361,-0.0028197840)

#paad
#model = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0619897730331895,0,0,0,0,0,0,0,-0.0652808661161353,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000242807448895044,0,0,0,0,0,0.000829885799158693,0,0,0,0.00599514623685536,0,0,0,0,0.0480362902362254,0,0,0,0,0,0,0,0)

if(file.exists(paste('~/',main_fold,'/lasso/coef.Rdata',sep='')))
{
  load(paste('~/',main_fold,'/lasso/coef.Rdata',sep=''))
  print('load....')
} else {
  coef=coef(fit, s = cvfit$lambda.min)
  if(length(coef[which(coef != 0),])==0)
  {
    stop('coef = 0')
  }
  save(coef,file = paste('~/',main_fold,'/lasso/coef.Rdata',sep=''))
  print('save....')
}

index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names=F)

trainFinalGeneExp=erged_df[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("Time","Status",lassoGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
num_out <- standardize(as.data.frame(lapply(erged_df[,lassoGene],as.numeric)))
rownames(num_out) <- rownames(erged_df)
outTab=cbind(erged_df[,c("Time","Status")],num_out,riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file=paste(name_fold,"Risk.txt",sep=''),sep="\t",quote=F,row.names=F)

################################################

library(timeROC)

bioROC=function(inputFile=null,rocFile=null){
  
  rt=read.table(inputFile,header=T,sep="\t")
  ROC_rt=timeROC(T=rt$Time,delta=rt$Status,
                 marker=rt$riskScore,cause=1,
                 weighting='aalen',
                 times=c(1,2,3),ROC=TRUE)
  pdf(file=rocFile,width=5,height=5)
  plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
  plot(ROC_rt,time=2,col='blue',add=TRUE,title=FALSE,lwd=2)
  plot(ROC_rt,time=3,col='red',add=TRUE,title=FALSE,lwd=2)
  legend('bottomright',
         c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
           paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
           paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
         col=c("green",'blue','red'),lwd=2,bty = 'n')
  dev.off()
  return(max(c(ROC_rt$AUC[1],ROC_rt$AUC[2],ROC_rt$AUC[3])))
}
roc <- bioROC(inputFile=paste(name_fold,"Risk.txt",sep=''),rocFile=paste(name_fold,".ROC.pdf",sep=''))

######################################

library(survminer)
bioSurvival=function(inputFile=null,outFile=null){
  rt=read.table(inputFile,header=T,sep="\t")                  
  diff=survdiff(Surv(Time, Status) ~ risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
  fit <- survfit(Surv(Time, Status) ~ risk, data = rt)
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
  pdf(file=outFile,onefile = FALSE,width = 5,height =4.5)
  print(surPlot)
  dev.off()
  return(pValue)
}
p_value <- bioSurvival(inputFile=paste(name_fold,"Risk.txt",sep=''),outFile=paste(name_fold,".survival.pdf",sep=''))

#######################################################################3

library(pheatmap)
bioRiskPlot=function(inputFile=null,riskScoreFile=null,survStatFile=null){
  rt=read.table(inputFile,sep="\t",header=T,row.names=1,check.names=F)  
  rt=rt[order(rt$riskScore),]    
  riskClass=rt[,"risk"]
  lowLength=length(riskClass[riskClass=="low"])
  highLength=length(riskClass[riskClass=="high"])
  lowMax=max(rt$riskScore[riskClass=="low"])
  line=rt[,"riskScore"]
  line[line>10]=10
  pdf(file=riskScoreFile,width = 8,height = 6)
  plot(line, type="p", pch=20,
       xlab="Patients (increasing risk socre)", ylab="Risk score",
       col=c(rep("lightblue",lowLength),rep("red",highLength)) )
  abline(h=lowMax,v=lowLength,lty=2)
  legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","lightblue"),cex=1.2)
  dev.off()
  color=as.vector(rt$Status)
  color[color==1]="red"
  color[color==0]="lightblue"
  pdf(file=survStatFile,width = 8,height = 6)
  plot(rt$Time, pch=19,
       xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","lightblue"),cex=1.2)
  abline(v=lowLength,lty=2)
  dev.off()
}

bioRiskPlot(inputFile=paste(name_fold,"Risk.txt",sep=''),riskScoreFile=paste(name_fold,".riskScore.pdf",sep=''),survStatFile=paste(name_fold,".survStat.pdf",sep=''))

##########################################################

risk.out <- read.table(paste(name_fold,'Risk.txt',sep=''),row.names = 1,header=T)
risk.out <- risk.out[order(risk.out$risk),]
risk.out <-subset(risk.out, select = -c(Time,Status,riskScore))
Group <- risk.out$risk
risk.out <-subset(risk.out, select = -risk)
Group <- data.frame(Risk = Group)
rownames(Group) = rownames(risk.out)
library(pheatmap)
pdf('phep.pdf',width=6,height=4)
mycolors <- colorRampPalette(c("darkgreen","white", "darkred"))(15)
tmp=as.data.frame(lapply(risk.out,as.numeric))
tmp <- as.matrix(t(tmp))
colnames(tmp) <- rownames(Group)
pheatmap(scale(tmp),annotation_col = Group,color = mycolors,cluster_rows=F,show_colnames=F,cluster_cols=F,cellheight = 10,cellwidth = 200/length(colnames(tmp)),border=F)
dev.off()
print(paste('p:',p_value,' roc:',roc,sep=''))