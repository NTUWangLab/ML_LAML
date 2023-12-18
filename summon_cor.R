
setwd('~/machine_learning/KEGG_GSVA')


length_graph <- '50'

file_name <- 'c2.KEGG'

temp_var <- t_results_c2.KEGG

temp_var$group <- rep(1,length(temp_var))

temp_var <- temp_var[order(-temp_var$t_value),]

eval(parse(text=paste('temp_var <- as.data.frame(t(',file_name,'[which(rownames(',file_name,')%in%rownames(temp_var[1:',length_graph,',])),]))',sep='')))

temp_var$SampleName <- rownames(temp_var)
temp_var$Group <- rep('Tumor',length(rownames(temp_var)))
temp_var$CODE <- rep('LAML',length(rownames(temp_var)))

write.csv(temp_var,file='KEGG.csv',row.names = F)
