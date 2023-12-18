
setwd('~/machine_learning/linkET')

library(linkET)
library(ggplot2)
library(dplyr)
library(GSVA)
library(data.table)

Sel <- 'High'

rf_risk <- read.csv('~/diff/ml_laml/rfrisk.csv')
rf_risk <- rf_risk[,-c(ncol(rf_risk))]
rf_group <- rf_risk[which(rf_risk$risk==Sel),]$SampleName



riskscore_geneset <- c("ALDOC","KLF9","PPM1F","CSTB","CYYR1","DDIT4","DGAT2","IFITM3","PMM1","SASH1","SH2D3C","SUSD3","ANXA4","ARAP1","CYC1","DHRS9","ECHDC3","FBXO6","FKBP5","HIP1","HPCAL1","IL6R","LTB4R","PLCB2","RAB3D","RGL2","RPS6KA1","S100A13","TRIB1","ADGRE5","PELO","AVPI1","CYP19A1","HLX","PKN1","PPM1N","SHARPIN","SLA","DGAT1","GAS6","KRT5","AK1","CBLN3","CCND3","ECE1","EHBP1L1","GADD45A","GLTP","GNG11","ICAM3","IGF2R","ISG20","LSP1","MAD2L1BP","MAST3","MYL6","MZT2A","NEDD9","NFKBIL1","OPTN","OSGIN1","PCGF5","PF4","PIM1","PRDX5","RELB","RHOC","RPL3L","SEC14L1","SESN2","SH3BP5","SH3TC1","SIAH2","SLC10A3","SRXN1","TFE3","TMEM63C","TREML2","TUBA4A","TWIST2")


##################exp
CIBER <- read.csv('~/machine_learning/data/tcga/survival_out.csv',header=T,check.names = F)
##################


##################CIBER
# CIBER <- read.csv('~/diff/data/CIBER.csv',header=T)
# CIBER <- CIBER[grep('LAML',CIBER$CODE),]
# CIBER <- CIBER[grep('TCGA',CIBER$SampleName),]
# CIBER <- CIBER[which(CIBER$SampleName%in%rf_group),]
# CIBER <- subset(CIBER,select=-c(Group,CODE))
##################












###################ssGSVA
dat <- fread(paste('~/TCGA_DATA/counts/LAML.txt',sep=''), sep = "\t",header = T,stringsAsFactors = F,check.names = F,na.strings="NA",data.table = F)
dat <- dat[!duplicated(dat[,c(1)]),]
rownames(dat) <- dat[,c(1)]
dat <- dat[,-c(1)]
dat <- na.omit(dat)
dat <- dat[,which(colnames(dat)%in%rf_group)]


ssample <- list('risk'=riskscore_geneset)
ssgsea <- gsva(as.matrix(dat),ssample, method='ssgsea', kcdf='Poisson',abs.ranking=TRUE)

ssgsea <-as.data.frame(t(ssgsea))
ssgsea$SampleName <- rownames(ssgsea)

out <- merge(CIBER,ssgsea,by='SampleName')
ssgsea <- data.frame(risk = out$risk)
out <- subset(out,select=-c(risk,SampleName))

###################



mantel <- mantel_test(ssgsea, out,
                      spec_select = list(ssGSVA = 1:1)) %>% 
mutate(rd = cut(r, breaks = c(-1, 0.2, 0.4, 1),
                labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),#对相关系数进行分割，便于映射大小
       pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))#对P值进行分割，便于映射颜色


pdf(paste(Sel,'.pdf',sep=''),height=20,width=20)
qcorrplot(correlate(out), type = "lower", diag = FALSE) +#热图绘制
  geom_square() +#热图绘制
  geom_couple(aes(colour = pd, size = rd),data = mantel,curvature = nice_curvature()) +#aes里面是线条格式，data对应的是mantel test 计算结果，curvature控制线条曲率
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"),limits=c(-1,1)) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "r",##guides()函数调整标签样式
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "R", order = 3))
dev.off()