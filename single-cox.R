setwd('~/machine_learning/cox')


library(dplyr)
library(survival)
library(data.table)
store_data = list()

for(sample_name in c('tcga','GSE12417','GSE37642','GSE146173','GSE106291'))
{
####
dat <- read.csv(paste('~/machine_learning/data/',sample_name,'/survival_out.csv',sep=''),header=T,row.name=1,check.names = F)
####
  
####
# dat <- fread(paste('~/machine_learning/data/',sample_name,'/all_data.txt',sep=''), sep = "\t",header = T,stringsAsFactors = F,check.names = F,na.strings="NA",data.table = F)
# dat <- dat[!duplicated(dat[,c(1)]),]
# rownames(dat) <- dat[,c(1)]
# dat <- dat[,-c(1)]
# dat <- na.omit(dat)
# dat <- as.data.frame(t(dat))
####

all_sur_data <- read.table(paste('~/machine_learning/data/',sample_name,'/LAML_patient.txt',sep=''),header=T)

output <- data.frame()

for(gene_name in colnames(dat))
{
  dat1 <- dat[which(colnames(dat)==gene_name)]
  
  dat1$SampleName = rownames(dat1)

  erged_df <- merge(dat1, all_sur_data, by = "SampleName")
  
  erged_df <- erged_df[,-c(1)]
  
  cox_model <- coxph(Surv(Time,Status) ~ as.numeric(erged_df[,1]),data=erged_df)
  
  lower = round(summary(cox_model)$conf.int[,3],2)
  upper = round(summary(cox_model)$conf.int[,4],2)
  hr = round(summary(cox_model)$conf.int[,1],2)
  p = summary(cox_model)$coefficients[,5]
  if(is.na(p))
  {
    next
  }
  if(p>0.05||hr<1)
  {
    next
  }
  ci = paste(hr,'(',lower,'|',upper,')',sep='')
  temp_data <- data.frame(genename = c(gene_name),pvalue = c(p),HR=c(hr),Lower=c(lower),Upper=c(upper),CI=c(ci))
  output <- rbind(output,temp_data)
  eval(parse(text = paste('store_data$',sample_name,"<-append(store_data$",sample_name,",'",gene_name,"')",sep='')))
}
library(ggplot2)


data <- data.frame(
  group = output$genename,
  HR = output$HR,
  p_value = output$pvalue,
  lower=output$Lower,
  upper = output$Upper
)

data <- data[order(data$HR),]

data$point_size <- data$p_value


data$group <- factor(data$group, levels = data$group)
data$HR[which(data$HR>2)] = 2
data$upper[which(data$upper>2)] = 2


forest_plot <- ggplot(data, aes(x=HR,y=group)) +
  geom_point(aes(size = point_size), color = "red3") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 2)) +
  labs(title = toupper(sample_name),
       x = "HR",
       size = "P-Value") +
  theme(plot.title = element_text(aes(x=1)))+
  theme_light()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12))

pdf(paste(sample_name,'.pdf',sep=''),height=length(output$genename)/3)
print(forest_plot)
dev.off()
}