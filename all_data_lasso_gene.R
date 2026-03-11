
data_exp<-read.table("/data/zrx/WATER/GSE184880_LASSO/data_expression.txt",sep = "\t",header = T)


allele <- rownames(GSE184880_momac_geneExp)
GSE184880_momac_geneExp_1<- cbind(allele,GSE184880_momac_geneExp)

GSE184880_momac_geneExp_1<-as.data.frame(GSE184880_momac_geneExp_1)
GSE184880_gene_name<-as.vector(GSE184880_momac_geneExp_1$allele)

data_expression_GSE_gene<-GSE184880_gene_name[GSE184880_gene_name%in%colnames(data_expression)]

data_clinical<-as.data.frame(data_clinical)
data_clinical<-na.omit(data_clinical)
data_expression<-data_expression[rownames(data_clinical),]


save(OS,file = "data_expression_GSE_gene_merge_clinical_OS.Rdata")

exprSet_hub <- t(data_expression[,data_expression_GSE_gene])
########https://www.jianshu.com/p/b4a19562ac5c
exprSet=exprSet_hub
meta=data_clinical
exprSet<-t(exprSet)
identical(rownames(exprSet),rownames(meta))
rownames(exprSet)%in%rownames(meta)
class(exprSet)
class(meta)
exprSet<-as.matrix(exprSet)

x=exprSet
class(meta)
x<-apply(x,2,as.numeric)
y=meta[,1]
y=meta$OS_STATUS

meta<-as.matrix(meta)
y<-as.numeric(y)
library(glmnet)
set.seed(123)
model_lasso <- glmnet(x, y,nlambda=10, alpha=1)
dev.off()
print(model_lasso)
plot(model_lasso, xvar = "lambda", label = T)
cv_fit <- cv.glmnet(x=x, y=y, nlambda = 1000,alpha = 1)
#这里选取了1000个，精确些
plot(cv_fit)
model_lasso_min <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
model_lasso_1se <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.1se)
choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
choose_gene_1se=rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]
length(choose_gene_min)  #70个
length(choose_gene_1se)  #40个
choose_gene_1se
lasso.prob <- predict(cv_fit, newx=x , s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
#如上得到根据模型预测每个样本的生存概率的单列矩阵
#载入数据，数据包括基本的临床信息和构建模型miRNA表达值
#高通量检测后会有很多差异miRNA，前期通过一系列分析进一步筛选出用于cox风险模型构建的miRNA
#这里默认已经完成了自变量的筛选，确定了8个miRNA用于cox风险比模型的构建
choose_gene_min
OS<-data_clinical
df_risk <- data_expression[,choose_gene_min]
df_risk <-df_risk[rownames(OS),]
meta<-apply(meta,2,as.numeric)
dat_risk <- cbind(df_risk, meta)
#构建多因素cox回归模型
#s=Surv(time, event) ~ CXCL10+JUN+CD2+CCR7+MAP3K13
s=Surv(OS_MONTHS, OS_STATUS) ~ RALA+IGKC+CX3CR1+RGS1+JUN+CD3E+CD3D+CCR7+LAMP3+TFPI2+NET1+TUBB2B+MARCKSL1+GADD45A+GSTP1+PPA1+CD40+NABP1
model <- coxph(s, data = dat_risk)
summary(model,data=dat)

#使用Survival程序包的Predict函数，计算出每位患者的风险评
RiskScore<-predict(model,type = "risk")
names(RiskScore) = rownames(dat_risk)


#开始绘制风险模型的生存点图
fp <- RiskScore
phe<-dat_risk
fp_dat=data.frame(patientid=1:length(fp),fp=as.numeric(sort(fp)))
#添加风险分组，以风险评分的中位值将患者分为c两组，大于中位值的 患者为高风险组，小于或等于中位值的患者为低风 险组
fp_dat$riskgroup= ifelse(fp_dat$fp>=median(fp_dat$fp),'high','low')

sur_dat=data.frame(patientid=1:length(fp),time=phe[names(sort(fp)),'OS_MONTHS'],event=phe[names(sort(fp)),'OS_STATUS']) 
sur_dat$event=ifelse(sur_dat$event==0,'alive','death')
sur_dat$event=factor(sur_dat$event,levels = c("death","alive"))
exp_dat=dat_risk[names(sort(fp)),1:(ncol(dat_risk)-2)]
#fp_dat用来绘制第一幅图
#sur_dat用来绘制第二幅图
#exp_dat用来绘制第三幅图
fp_dat$fp <- log(fp_dat$fp)
###第一个图
library(ggplot2)
p1=ggplot(fp_dat,aes(x=patientid,y=fp))+geom_point(aes(color=riskgroup))+
  scale_colour_manual(values = c("red","green"))+
  theme_bw()+labs(x="Patient ID(increasing risk score)",y="Risk score")+
  geom_hline(yintercept=median(fp_dat$fp),colour="black", linetype="dotted",size=0.8)+
  geom_vline(xintercept=sum(fp_dat$riskgroup=="low"),colour="black", linetype="dotted",size=0.8)
p1

#第二个图
p2=ggplot(sur_dat,aes(x=patientid,y=time))+geom_point(aes(col=event))+theme_bw()+
  scale_colour_manual(values = c("red","green"))+
  labs(x="Patient ID(increasing risk score)",y="Survival time(year)")+
  geom_vline(xintercept=sum(fp_dat$riskgroup=="low"),colour="black", linetype="dotted",size=0.8)
p2

#第三个图
fp_dat$riskgroup
dim(df_risk)
an = data.frame(group = fp_dat$riskgroup,
                row.names = rownames(df_risk))
library(pheatmap)
mycolors <- colorRampPalette(c("white", "green", "red"), bias = 1.2)(100)

tmp=t(apply(df_risk, 2, as.numeric))
colnames(tmp)<-rownames(df_risk)
#tmp[tmp > 1] = 1
#tmp[tmp < -1] = -1
an$group<-factor(an$group,levels = c("high","low"))

p3 <- pheatmap(tmp,col= mycolors,annotation_col = an,show_colnames = F,cluster_cols = F)
p3

#拼图实现三图联动
library(ggplotify)
plots = list(p1,p2,as.ggplot(as.grob(p3)))
library(gridExtra)
lay1 = rbind(c(rep(1,7)),c(rep(2,7)),c(rep(3,7))) #布局矩阵
grid.arrange(grobs = plots, layout_matrix = lay1, heigths = c(2,3,2),weights=c(10,10,10))

############
risk_sv <- as.data.frame(RiskScore)
risk_sv$barcode <- rownames(risk_sv)

identical(rownames(meta),rownames(risk_sv))
#####
TCGA_surv <- cbind(meta,RiskScore)
TCGA_surv$time <- as.numeric(TCGA_surv$time)
my.surv <-Surv(as.numeric(TCGA_surv$time), TCGA_surv$event)
TCGA_surv$group=ifelse(TCGA_surv$RiskScore > median(TCGA_surv$RiskScore),"high","low")
kmfit1 <- survfit(my.surv~group, data = TCGA_surv)
ggsurvplot(kmfit1, conf.int = F, pval = T, risk.table = T, ncesor.plot = T)

##############画表达箱式图
library(tinyarray)
group_list = fp_dat$riskgroup #定义组别
group_list = factor(group_list,levels = c("low","high")) # 转换为因子
table(group_list) #查看各自的数量
df_risk1<-as.numeric(df_risk)
group_list
df_risk=apply(df_risk,2,as.numeric)
draw_boxplot(t(df_risk),group_list)

alldata_tRNA_risk <- cbind(dat_risk, risk_sv)
save(alldata_tRNA_risk, file = 'alldata_OV_macro_riskscore.Rdata')

##########
meta$riskscore <- fp
xx <- rownames(meta)
pdata_TCGA_OV_342 <- TCGA_OV_clinical
pdata_TCGA_OV_342<- pdata_TCGA_OV_342[xx,]
cox_data <- cbind(pdata_TCGA_OV_342,meta)
cox_data$stage <- ifelse(cox_data$figo_stage %in%  c("Stage IIA","Stage IIB","Stage IIC"),"stageII","stageIII-IV")
cox_data$race1 <- ifelse(cox_data$race %in%  c("white"),"white","others")
cox_data$age <- ifelse(cox_data$age_at_diagnosis > 21550, ">60", "<60") 
saveRDS(cox_data, "cox_data_TCGA_OV_342.rds")
########
#单因素cox回归分析，这里看性别sex这个特征
colnames(cox_data)
res.cox <- coxph(Surv(time, event) ~ riskscore, data = cox_data)
res.cox

ggforest(res.cox, data = cox_data,  #数据集
         main = 'Hazard ratio',  #标题
         cpositions = c(0.05, 0.15, 0.35),  #前三列距离
         fontsize = 1, #字体大小
         noDigits = 3 #保留HR值以及95%CI的小数位数
)

colnames(cox_data)

#批量单因素cox回归分析
covariates <- c("stage","race1","age","riskscore")
#分别对每一个变量，构建生存分析的公式
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, event)~', x)))
univ_formulas
#循环对每一个特征做cox回归分析
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = cox_data)})
univ_models

#提取HR，95%置信区间和p值
univ_results1 <- lapply(univ_models,
                        function(x){ 
                          x <- summary(x)
                          #获取p值
                          p.value<-signif(x$wald["pvalue"], digits=2)
                          #获取HR
                          HR <-signif(x$coef[2], digits=2);
                          #获取95%置信区间
                          HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                          HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                          HR1 <- paste0(HR, " (", 
                                        HR.confint.lower, "-", HR.confint.upper, ")")
                          res<-c(p.value,HR1,HR.confint.lower,HR.confint.upper,HR)
                          names(res)<-c("p.value","HR1","low.HR","up.HR","HR")
                          return(res)
                        })
#转换成数据框，并转置
res1 <- t(as.data.frame(univ_results1, check.names = FALSE))
res1 <- as.data.frame(res1)
write.csv(res1, "基因单因素.csv")

# 森林图
#install.packages("forestplot")
library(forestplot)
out_multi <- read.csv("基因单因素.csv", header = T, row.names = 1)


tabletext <- cbind(c(NA,"Risk Factors",rownames(res1)),
                   c(NA,"P value",res1$p.value),
                   c(NA,"Hazard Ratio(95% CI)",res1$HR1))
tabletext

forestplot(labeltext=tabletext, 
           graph.pos=3,  #为Pvalue箱线图所在的位置
           col=fpColors(box="#D55E00", lines="#CC79A7", zero = "gray50"),
           mean=c(NA,NA,res1$HR),
           lower=c(NA,NA,res1$low.HR), #95%置信区间下限
           upper=c(NA,NA,res1$up.HR), #95%置信区间上限
           boxsize=0.3,lwd.ci=2,   #箱子大小，线的宽度
           ci.vertices.height = 0.2,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=1,      #zero线宽 基准线的位置
           colgap=unit(7,"mm"),    #列间隙
           xticks = c(-20,-10,0, 10,20), #横坐标刻度
           lwd.xaxis=1,            #X轴线宽
           lineheight = unit(1.5,"cm"), #固定行高
           graphwidth = unit(.3,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines=list("2" = gpar(lwd=4, col="black"),
                           "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
                           "7" = gpar(lwd=4, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           mar=unit(rep(0.5, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1.2),
                          ticks=gpar(cex=1.2),
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.2)),
           xlab="Hazard Ratio",
           title = "Univariate cox regression analysis")

####多因素cox回归分析
table(cox_data$figo_stage)
table(cox_data$race)

res.cox <- coxph(Surv(time, event) ~  age + race1 + stage + riskscore, data =  cox_data)
x <- summary(res.cox)
pvalue=signif(as.matrix(x$coefficients)[,5],2)
HR=signif(as.matrix(x$coefficients)[,2],2)
low=signif(x$conf.int[,3],2)
high=signif(x$conf.int[,4],2)
multi_res=data.frame(p.value=pvalue,
                     HR=paste(HR," (",low,"-",high,")",sep=""),
                     stringsAsFactors = F
)
multi_res
write.table(file="multivariate_cox_result.txt",multi_res,quote=F,sep="\t")

ggforest(res.cox, data = cox_data,  #数据集
         main = 'Hazard ratio',  #标题
         cpositions = c(0.05, 0.15, 0.35),  #前三列距离
         fontsize = 1, #字体大小
         noDigits = 3 #保留HR值以及95%CI的小数位数
)

####ROC分析
?timeROC
library(timeROC)
library(survival)
install.packages("timeROC")
new_dat <- meta
new_dat$Riskscore <- Riskscore

saveRDS(new_dat,"new_dat_for_ROC.rds")
result <-with(meta, timeROC(T=OS_MONTHS,
                            delta=OS_STATUS,
                            marker=RiskScore,
                            cause=1,
                            times=c(365,1095,1825),
                            iid = TRUE))

#identical(c(result$TP[,1],result$TP[,2],result$TP[,3]),as.numeric(result$TP))
dat = data.frame(fpr = as.numeric(result$FP),
                 tpr = as.numeric(result$TP),
                 time = rep(as.factor(c(365,1095,1825)),each = nrow(result$TP)))

library(ggplot2)
ggplot() + 
  geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 1) + 
  scale_color_manual(name = NULL,values = c("#92C5DE", "#F4A582", "#66C2A5"),
                     labels = paste0("AUC of ",c(1,3,5),"-y survival: ",
                                     format(round(result$AUC,2),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()

######nomogram
library("survival")
library("survminer")
library(rms)
###单因素Cox回归分析
cox_data1<-cox_data[,c("age", "race1", "stage")]
dat_risk_nomogram <- cox_data1
dat_risk_nomogram$time <- meta$time
dat_risk_nomogram$event <- meta$event
dat_risk_nomogram$riskscore <- fp

#构建多因素cox回归模型
library(rms)
dd <- datadist(dat_risk_nomogram) ###转换数据格式
options(datadist="dd")
colnames(dat_risk_nomogram)
s=Surv(time, event) ~ race1+stage+age+riskscore
model <- coxph(s, data = dat_risk_nomogram )
summary(model,data=dat)

##########
ggforest(model, data = dat_risk_nomogram,  #数据集
         main = 'Hazard ratio',  #标题
         cpositions = c(0.05, 0.15, 0.35),  #前三列距离
         fontsize = 1, #字体大小
         noDigits = 3 #保留HR值以及95%CI的小数位数
)

coxpbc<-cph(formula = Surv(time, event) ~  race1+stage+age+riskscore ,data=dat_risk_nomogram,x=T,y=T,surv = T,na.action=na.delete)  #,time.inc =2920

coxm1 <- cph(Surv(time,event)~race1+stage+age+riskscore,x=T,y=T,data=dat_risk_nomogram,surv=T)
surv <- Survival(coxm1)
surv1 <- function(x)surv(1*365,lp=x)
surv2 <- function(x)surv(1*1095,lp=x)
surv3 <- function(x)surv(1*1825,lp=x)

nom1<-nomogram(coxm1,fun=list(surv1,surv2,surv3),lp = F,
               funlabel=c('1-Year Survival probability',
                          '3-Years survival probability',
                          '5-Years survival probability'),
               maxscale=100,
               fun.at=c('0.9','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1'))
plot(nom1)

######Calibration curve https://www.jianshu.com/p/e1259488ece2
###https://www.jianshu.com/p/0512bc3e3be9
#构建回归模型
## 参数说明：
## 1、绘制校正曲线前需要在模型函数中添加参数x=T, y=T，详细参考帮助
## 2、u需要与之前模型中定义好的time.inc一致，即365或730；
## 3、m要根据样本量来确定，由于标准曲线一般将所有样本分为3组（在图中显示3个点）,而m代表每组的样本量数，因此m*3应该等于或近似等于样本量；
## 4、b代表最大再抽样的样本量
f1<-cph(Surv(time,event)~race1+stage+age+riskscore,
        x=T,y=T,data=dat_risk_nomogram,surv=T,na.action=na.delete,time.inc = 365) 

#参数m=50表示每组50个样本进行重复计算
cal1<-calibrate(f1, cmethod="KM", method="boot",u=365,m=112,B=338) 
plot(cal1,lwd=2,lty=1,
     errbar.col = c("#2166AC"),
     xlim = c(0,1),ylim= c(0,1),
     xlab="Nomogram-Predicted Probability of 1-Year OS",
     ylab="Actual 1-Year DFS (proportion)",
     col=c(rgb(192,98,83,maxColorValue=255)))

#############
f2<-cph(Surv(time,event)~race1+stage+age+riskscore,
        x=T,y=T,data=dat_risk_nomogram,surv=T,na.action=na.delete,time.inc = 1095) 

#参数m=50表示每组50个样本进行重复计算
cal2<-calibrate(f2, cmethod="KM", method="boot",u=1095,m=112,B=338) 
plot(cal2,lwd=2,lty=1,
     errbar.col = c("#2166AC"),
     xlim = c(0,1),ylim= c(0,1),
     xlab="Nomogram-Predicted Probability of 2-Year OS",
     ylab="Actual 2-Year DFS (proportion)",
     col=c(rgb(192,98,83,maxColorValue=255)))

#############
f3<-cph(Surv(time,event)~race1+stage+age+riskscore,
        x=T,y=T,data=dat_risk_nomogram,surv=T,na.action=na.delete,time.inc = 1825) 

#参数m=50表示每组50个样本进行重复计算
range(meta$time)
cal3<-calibrate(f3, cmethod="KM", method="boot",u=1825,m=112,B=338) 
plot(cal3,lwd=2,lty=1,
     errbar.col = c("#2166AC"),
     xlim = c(0,1),ylim= c(0,1),
     xlab="Nomogram-Predicted Probability of 2-Year OS",
     ylab="Actual 3-Year DFS (proportion)",
     col=c(rgb(192,98,83,maxColorValue=255)))

pdf("calibration_compare.pdf",width = 8,height = 8)
plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal2,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

plot(cal3,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","2-year","3-year"), #图例文字
       col =c("#2166AC","#B2182B","#224444"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
dev.off()
dev.new()
