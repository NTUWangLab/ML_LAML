library(Boruta)
library(survival)

set.seed(123)
library(dplyr)
setwd('~/machine_learning/single_algorithm/Boruta')

all_sur_data <- read.table('~/machine_learning/data/tcga/LAML_patient.txt',header=T,check.names = F)

gene_exp <- read.csv('~/machine_learning/data/tcga/survival_out.csv',header=T,check.names = F)

mixed <- merge(gene_exp, all_sur_data, by = "SampleName")

rownames(mixed) <- mixed$SampleName

time <- mixed$Time
status <- mixed$Status

mixed <- select(mixed,-c(SampleName,Time,Status))

features <- mixed

surv_object <- Surv(time = time, event = status)

boruta_result <- Boruta(features, surv_object)

selected_features <- getSelectedAttributes(boruta_result,withTentative = TRUE)

Fun <- function(x){
  return ((x - min(x))/(max(x)-min(x)))
}

Boruta_imp <- attStats(boruta_result)[order(attStats(boruta_result)$meanImp),]
BORUTA_order <- data.frame(BORUTA=Fun(Boruta_imp$meanImp),gene=rownames(Boruta_imp))

save(BORUTA_order,file='~/machine_learning/single_algorithm/BORUTA_order.Rdata')

pdf('boruta_importance.pdf',width = 20,height=10)
par(oma=c(3,3,3,3)) 
plot(boruta_result,las=2,xlab='')
legend(x = 'topleft', 
       legend = c(paste('P-value:',boruta_result$pValue),sep=''),
       lty = 0,
       bty = 'n')
dev.off()

pdf('boruta_history.pdf',width = 20,height=10)
par(oma=c(3,3,3,3)) 
plot(plotImpHistory(boruta_result),las=2)
legend(x = 'topleft', 
       legend = c(paste('P-value:',boruta_result$pValue),sep=''),
       lty = 0,
       bty = 'n')
dev.off()