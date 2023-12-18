library(survival)
library(survivalsvm)


setwd('~/machine_learning')

all_sur_data <- read.table('~/machine_learning/data/tcga/LAML_patient.txt',header=T,check.names = F)

gene_exp <- read.csv('~/machine_learning/data/tcga/survival_out.csv',header=T,check.names = F)

mixed <- merge(gene_exp, all_sur_data, by = "SampleName")
rownames(mixed) <- mixed$SampleName

# for(i in 1:length(mixed$Status))
# {
#   if(mixed[i,length(mixed)] == 1)
#   {
#     mixed[i,length(mixed)-1] = -mixed[i,length(mixed)-1]
#   }
# }
mixed<- subset(mixed, select = -c(SampleName))

model <- survivalsvm(Surv(Time, Status) ~ ., mixed,  gamma.mu = 1,opt.meth = "ipop", kernel = "add_kernel")

