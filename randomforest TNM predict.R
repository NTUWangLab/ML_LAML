

library(dplyr)
library(ggplot2)
library(ggpubr)
library(pROC)
library(randomForestSRC)
library(boot)

setwd('~/machine_learning/RF_TNM')

name_fold = 'tcga'

all_sur_data <- read.table(paste('~/machine_learning/data/',name_fold,'/LAML_patient.txt',sep=''),header=T,check.names = F)
gene_exp <- read.csv(paste('~/machine_learning/data/',name_fold,'/survival_out.csv',sep=''),header=T,check.names = F)
mixed_sur <- merge(gene_exp, all_sur_data, by = "SampleName")

all_TNM_data <- read.table(paste('~/machine_learning/data/',name_fold,'/',name_fold,'fenqi.txt',sep=''),header=T,check.names = F)
mixed_TNM <- merge(gene_exp, all_TNM_data, by = "SampleName")

result <- data.frame()

mixed_sur<- subset(mixed_sur, select = -c(SampleName))



set.seed(123456)

p.obj <- rfsrc(Surv(Time,Status)~.,data = mixed_sur,
               mtry=3,
               nodesize=5,
               ntree=2000,
               importance = T,
               splitrule = 'logrank',
               proximity = T,
               forest = T,
               seed=123456
)

rs <- cbind(mixed_sur[,c(ncol(mixed_sur),ncol(mixed_sur)-1)],RS=as.numeric(predict(p.obj,newdata = mixed_sur)$predicted))
cc <- data.frame(p=as.numeric(summary(coxph(Surv(Time,Status)~RS,rs))$coefficients[,5]),C=as.numeric(summary(coxph(Surv(Time,Status)~RS,rs))$concordance[1]))
rownames(cc) <- c('Risk')
result <- rbind(result,cc)

library(timeROC)
for(i in colnames(all_TNM_data)[-c(1)])
{
  #compare_means(Age ~ RS, data = rs,method = "t.test")
  mixed_temp <- merge(gene_exp,mixed_TNM[which(colnames(mixed_TNM)%in%c('SampleName',i))],by='SampleName')
  mixed_tempp <- merge(all_sur_data,mixed_temp[which(colnames(mixed_TNM)%in%c('SampleName'))],by='SampleName')
  mixed_temp <- mixed_temp[which(mixed_temp$SampleName%in%c(mixed_tempp$SampleName)),]
  mixed_temp <- mixed_temp[,-c(1)]
  eval(parse(text = paste('t.obj <- rfsrc(',i,'~.,data = mixed_temp,mtry=3,nodesize=5,ntree=2000,importance = T,proximity = T,forest = T,seed=123456)',sep='')))
  rs <- cbind(subset(mixed_tempp,select=c('Time','Status')),RS=as.numeric(predict(t.obj,newdata = mixed_temp)$predicted))
  cc <- data.frame(p=as.numeric(summary(coxph(Surv(Time,Status)~RS,rs))$coefficients[,5]),C=as.numeric(summary(coxph(Surv(Time,Status)~RS,rs))$concordance[1]))
  rownames(cc) <- c(i)
  result <- rbind(result,cc)
}

result$Category <- rownames(result)
result$Category <- factor(result$Category, levels = c('Risk',result$Category[which(result$Category!='Risk')]))

colors <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","orange2","#D0006E","#B1F100","#3914AF","#009999")
custom_colors <- c()
for(i in result$Category)
{
  eval(parse(text= paste('custom_colors <- append(custom_colors,c(',i,' = colors[which(result$Category=="',i,'")]))',sep='')))
}

result <- result %>%
  mutate(p_label = case_when(
    p < 0.0001 ~ "****",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

pdf(paste(name_fold,'.pdf',sep=''),width=length(result$Category)/1.3,height=5)
ggplot(result, aes(x = Category, y = C, fill = Category)) +
  geom_bar(stat = "identity",width=0.5) +
  theme_classic() +
  labs(title = name_fold, x = "Category", y = "C-index") +
  geom_text(aes(label = p_label), vjust = -0.5)+
  scale_fill_manual(values = custom_colors)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()
