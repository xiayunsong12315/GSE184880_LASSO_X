rm(list = ls())

data_expression<-read.table("/data/zrx/WATER/GSE184880_LASSO/data_expression.txt",header = T,sep='\t',encoding='UTF-8')
data_expression = as.matrix(data_expression)#设置成矩阵数据
data_expression<-t(data_expression)

data_clinical<-read.table("/data/zrx/WATER/GSE184880_LASSO/data_clinical.txt",header = T,sep="\t",encoding='UTF-8')
OS = as.matrix(OS)#设置成矩阵数据
OS<-t(as.matrix(OS))

data_expression<-data_expression[-1,]
colnames(data_expression)<-data_expression[1,]
data_expression<-as.data.frame(data_expression)
data_expression<-data_expression[!duplicated(rownames(data_expression)), ]

rownames(data_expression)<-substring(rownames(data_expression),1,12)
rownames(data_expression)<-gsub("\\.","\\-",rownames(data_expression))
data_expression$ID<-rownames(data_expression)
train_data$ID<-rownames(train_data)
data_merge<-merge(data_expression,train_data,by="ID")
rownames(train_data_merge)<-train_data_merge[,1]
train_data_merge<-train_data_merge[,-1]

#train_data<-data_expression[,rownames(train_data)]
test_OS<-gene18[rownames(test_data),]
test_OS<-test_OS[,-(3:20)]
train_data_merge<-t(as.matrix(train_data_merge))
#GSE184880_momac_genelist_seurat_clusters<-as.matrix(GSE184880_momac_genelist_seurat_clusters)

#library(ConsensusClusterPlus)

#data_expression_macro<-data_expression[GSE184880_momac_genelist_seurat_clusters$GSE184880_momac_genelist_seurat_clusters,]
data_expression_macro<-data_expression[,GSE184880_momac_genelist_seurat_clusters[GSE184880_momac_genelist_seurat_clusters[,1]%in%colnames(data_expression),1]]
#gene18_macro<-data_expression_macro
#train_data_macro<-na.omit(train_data_macro)
data_expression_GSE_gene<-
data_expression_macro<-t(as.matrix(data_expression_macro))

colnames(data_expression_macro)%in%data_expression_GSE_gene
data_expression_GSE_gene%in%colnames(data_expression_macro)

write.table(data_expression_macro,"data_expression_macro.txt",sep = "\t",header = T)

Cluster <- ConsensusClusterPlus(d = data_expression_macro, # 分析矩阵
                                maxK = 7,  # 最大聚类数目
                                reps = 1000, # 重抽样的次数
                                pItem = 0.8, # 样品的重抽样比例
                                clusterAlg = "pam", # 使用的聚类算法，可以选择"hc"(hclust), "pam", "km"(k-means)
                                innerLinkage = "ward.D2", 
                                finalLinkage = "ward.D2",
                                distance = "euclidean",  # 计算距离的方法，可以选择pearson、spearman、euclidean、binary、maximum、canberra、minkowski
                                seed = 123456, # 设置随机种子，方便重复
                                plot = "pdf", # 结果图片的导出类型，可以选择"png"或者"pdf"
                                title = "(data_expression_GSE_gene聚类",writeTable=T)

######
# 重新绘制热图
annCol_macro2 <- data.frame(Cluster = paste0("Cluster",
                                             Cluster[[2]][["consensusClass"]]),
                            row.names = colnames(data_expression_macro))
head(annCol_macro2)
##  Cluster
#  Cluster
#Cluster
#TCGA-04-1331-01 Cluster1
#TCGA-04-1332-01 Cluster1
#TCGA-04-1337-01 Cluster2
#TCGA-04-1343-01 Cluster2
#TCGA-04-1348-01 Cluster2
#TCGA-04-1350-01 Cluster1                        
#install.packages("RColorBrewer")
mycol <- brewer.pal(2, "Set1")
annColors <- list(Cluster = c("Cluster1" = mycol[1],
                              "Cluster2" = mycol[2]))
heatdata <- Cluster[[2]][["consensusMatrix"]]
dimnames(heatdata) <- list(colnames(data_expression_macro), colnames(data_expression_macro))
heatdata[1:3, 1:3]
##TCGA-04-1331-01 TCGA-04-1332-01 TCGA-04-1337-01
##TCGA-04-1331-01       1.0000000       0.9890625       0.3011094
#TCGA-04-1332-01       0.9890625       1.0000000       0.2793210
#TCGA-04-1337-01       0.3011094       0.2793210       1.0000000
# 绘制热图
pheatmap(mat = heatdata,
         color = colorRampPalette((c("white","steelblue")))(100),
         border_color = NA,
         annotation_col = annCol_macro2,
         annotation_colors = annColors,
         show_colnames = F,
         show_rownames = F)

## png 
##   2

# 取出分型结果
Cluster_data_expression_macro_TCGA <- annCol_macro2 %>%
  as.data.frame() %>%
  rownames_to_column("patient_ID")
table(Cluster_gene18_macro_TCGA$Cluster)
write.csv(Cluster_train_macro_TCGA,'Cluster_train_macro_TCGA.csv')



gene18_OS<-gene18[,-(3:20)]
#library(tidyverse)
#library(dplyr)
###TCGA survival
data_clinical
identical(Cluster_data_expression_macro_TCGA$patient_ID,rownames(data_clinical))
Cluster_data_expression_macro_TCGA_1<- cbind(Cluster_data_expression_macro_TCGA,data_clinical)
head(Cluster_test_macro_TCGA_1)
test_OS<-test_data[,-(3:10)]

Cluster_data_expression_macro_TCGA_1$OS_MONTHS<-as.numeric(Cluster_data_expression_macro_TCGA_1$OS_MONTHS)
#library(survminer)
#library(survival)
fit <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ Cluster, 
               data =Cluster_data_expression_macro_TCGA_1 )
p<-ggsurvplot(fit,
              pval = TRUE,
              linetype = "solid",  
              palette = mycol,
              surv.median.line = "hv", 
              title = "Overall survival",
              ylab = "Cumulative survival (percentage)",
              xlab = " Time (Years)",
              legend.title = "Survival Plot",
              legend = c(0.85,0.70),
              risk.table = T,
              tables.height = 0.2,
              ggtheme = theme_bw(),
              tables.theme = theme_cleantable()
)
p
dev.off()
dev.copy2pdf(file = "TCGA_HNSC_survival.pdf",width = 9, height = 6) # 保存图片 <- na.omit(test_data)