
library(ggplot2)
library(dplyr)
setwd('~/machine_learning/single_algorithm/draw_plot')

load('~/machine_learning/single_algorithm/BORUTA_order.Rdata')
load('~/machine_learning/single_algorithm/XGBOOST_order.Rdata') #少 c('KRT5','GNG11','ICAM3','IGF2R')
load('~/machine_learning/single_algorithm/ANN_order.Rdata')
load('~/machine_learning/single_algorithm/SVM_order.Rdata')
load('~/machine_learning/single_algorithm/RF_order.Rdata')
XGBOOST_order <- rbind(XGBOOST_order,data.frame(gene=c('KRT5','GNG11','ICAM3','IGF2R'),XGBOOST=c(0,0,0,0)))

out <- SVM_order %>% left_join(ANN_order,by='gene') %>% left_join(BORUTA_order,by='gene') %>% left_join(RF_order,by='gene') %>% left_join(XGBOOST_order,by='gene')
rownames(out) <- out$gene
out <- subset(out,select=-c(gene))
out <- cbind(out,.Mean = rowSums(out)/length(colnames(out)))
out <- out[order(-out$.Mean),]*3
value <- c()

for(i in 1:length(colnames(out)))
{
  eval(parse(text=paste('value<-append(value,',out[,i],')',sep='')))
}

color <- c()
COLORS <- c("#70f3ff","#44cef6","#3eede7","#1685a9","#177cb0","#065279","#003472","#4b5cc4","#a1afc9","#2e4e7e","#3b2e7e","#4a4266","#426666","#425066","#574266","#8d4bbb","#815463","#815476","#4c221b","#003371","#56004f","#801dae","#4c8dae","#b0a4e3","#cca4e3","#edd1d8","#e4c6d0","#ff461f","#ff2d51","#f36838","#ed5736","#ff4777","#f00056","#ffb3a7","#f47983","#db5a6b","#c93756","#f9906f","#f05654","#ff2121","#f20c00","#8c4356","#c83c23","#9d2933","#ff4c00","#ff4e20","#f35336","#dc3023","#ff3300","#cb3a56","#a98175","#b36d61","#ef7a82","#ff0097","#c32136","#be002f","#c91f37","#bf242a","#c3272b","#9d2933","#60281e","#622a1d","#bce672","#c9dd22","#bddd22","#afdd22","#a3d900","#9ed900","#9ed048","#96ce54","#00bc12","#0eb83a","#0eb83a","#0aa344","#16a951","#21a675","#057748","#0c8918","#00e500","#40de5a","#00e079","#00e09e","#3de1ad","#2add9c","#2edfa3","#7fecad","#a4e2c6","#7bcfa6","#1bd1a5","#48c0a3","#549688","#789262","#758a99","#50616d","#424c50","#41555d","#eaff56","#fff143","#faff72","#ffa631","#ffa400","#fa8c35","#ff8c31","#ff8936","#ff7500","#ffb61e","#ffc773","#ffc64b","#f2be45","#f0c239","#e9bb1d","#d9b611","#eacd76","#eedeb0","#d3b17d","#e29c45","#a78e44","#c89b40","#ae7000","#ca6924","#b25d25","#b35c44","#9b4400","#9c5333","#a88462","#896c39","#827100","#6e511e","#7c4b00","#955539","#845a33","#ffffff","#e9e7ef")


for(j in 1:length(colnames(out)))
{
  eval(parse(text=paste('color <- append(color,c(',colnames(out)[j],'="',sample(COLORS,size=1),'"))',sep='')))
}

my_data <- data.frame(
  Category = rep(c(colnames(out)[1], colnames(out)[2], colnames(out)[3], colnames(out)[4],colnames(out)[5],colnames(out)[6]), each = length(out[,1])),
  Sample = rep(1:80, times = length(colnames(out))),
  Value = value,
  Color = color,
  name = rownames(out)
)

color[6] = "#4c221b"

pdf('algorithm.pdf',width=20,height=10)
ggplot(my_data, aes(x = as.factor(Sample), y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  facet_grid(Category ~ ., scales = "free_y", space = "free") +
  scale_fill_manual(values = color) +
  labs(title = "plot", x = "Gene", y = "Important") +
  scale_x_discrete(labels = my_data$name) +  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
