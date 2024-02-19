
#系统报错改为英文
Sys.setenv(LANGUAGE = "en")
#禁止转化为因子
options(stringsAsFactors = FALSE)
#清空环境
rm(list=ls())
#设置工作目录
setwd("F:/1/5LUSC")

#加载R包
library(rjson)
library(limma)   
library(GEOquery) 
library(stringr)
library(survival)
library(glmnet)
library(survminer)
library(timeROC)
library(data.table)
library(dplyr)
library(patchwork)
library(Matrix)
library(readr)
library(tibble)
library(ggplot2)
library(tidyverse) 
library(future)
library(pheatmap)
library(msigdbr)
library(clusterProfiler)
library(ggpubr)
library(tidyverse)
library(rjson)
library(limma)   
library(GEOquery) 
library(stringr)
library(survival)
library(glmnet)
library(survminer)
library(timeROC)
library(Seurat)
library(data.table)
library(dplyr)
library(patchwork)
library(Matrix)
library(readr)
library(tibble)
library(ggplot2)
library(SingleR)
library(tidyverse) 
library(future)
library(harmony)
library(pheatmap)
library(DoubletFinder)
library(devtools)
library(copykat)
library(GSVA)
library(AUCell)
library(msigdbr)
library(clusterProfiler)
library(monocle)
library(CellChat)
library(SCENIC)
library(ggpubr)
library(RColorBrewer)

#"MediumSeaGreen","Firebrick3"
#"Orange2","DarkOrchid"
#"DarkOrchid","Orange2"
#"Firebrick2","DodgerBlue1"
#"DarkOrchid","Orange2"


####1.1TCGA-LUSC下载整理####

#设置路径
setwd("F:/1/5LUSC/1TCGA-LUSC")
library(tidyverse)#读入meta.data文件
json <- jsonlite::fromJSON("metadata.cart.2023-06-13.json")
#查看
#View(json)
#获取需要的数据
sample_id <- sapply(json$associated_entities,function(x){x[,1]})
file_sample <- data.frame(sample_id,file_name=json$file_name)  
#读入counts文件
count_file <- list.files('gdc_download_20230613_075345.926514/',
                         pattern = '*.tsv',recursive = TRUE)
#整理
count_file_name <- strsplit(count_file,split='/')
count_file_name <- sapply(count_file_name,function(x){x[2]})
#必要时修改下面的基因数
matrix = data.frame(matrix(nrow=60660,ncol=0))
#下面的修改样本例数
for (i in 1:553){
  path = paste0('gdc_download_20230613_075345.926514//',count_file[i])   #Counts文件夹名
  data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
  colnames(data)<-data[2,]
  data <-data[-c(1:6),]
  #数据类型，选择其中之一 3：unstranded；4：stranded_first；5：stranded_second；6：tpm_unstranded；7：fpkm_unstranded；8：fpkm_uq_unstranded
  #data <- data[3]
  data <- data[6]
  #data <- data[7]
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix <- cbind(matrix,data)
}
#选择3，即为counts
#write.csv(matrix,'1TCGA/LUAD_Counts.csv',row.names = TRUE)
#选择6，即为TPM
#write.csv(matrix,'TPM.csv',row.names = TRUE)
#选择7，即为FPKM
#write.csv(matrix,'1TCGA/LUAD_FPKM.csv',row.names = TRUE)

#转化为gene_symbol
path = paste0('gdc_download_20230613_075345.926514//',count_file[1])
data<- as.matrix(read.delim(path,fill = TRUE,header = FALSE,row.names = 1))
gene_name <-data[-c(1:6),1]
matrix0 <- cbind(gene_name,matrix)
gene_type <- data[-c(1:6),2]
matrix0 <- cbind(gene_type,matrix0)
#将gene_name列去除重复的基因，保留每个基因最大表达量结果
matrix0 <- aggregate( . ~ gene_name,data=matrix0, max)    
#保留mRNA
matrix0 <- subset(x = matrix0, gene_type == "protein_coding")
#将gene_name列设为行名
rownames(matrix0) <- matrix0[,1]
matrix0 <- matrix0[,-c(1,2)]
#write.csv(matrix0,'1TCGA/LUAD_Counts_GeneSymbol.csv',row.names = TRUE)
#write.csv(matrix0,'1TCGA/LUAD_TPM_GeneSymbol.txt',row.names = TRUE)
#write.csv(matrix0,'1TCGA/LUAD_FPKM_GeneSymbol.csv',row.names = TRUE)
#write.table(matrix0,'TCGA.txt', sep="\t", quote=F, row.names = TRUE)

matrix1 = data.frame(ID=rownames(matrix0),matrix0)
#colnames(matrix1)
name123 = gsub('[.]', '-', colnames(matrix1))
colnames(matrix1) = name123
#write.table(matrix1,'TCGA_TPM.txt', sep="\t", quote=F, row.names = F)
write.table(matrix1,'TCGA_LUSC.txt', sep="\t", quote=F, row.names = F)

####2.1TCGA-LUAD下载整理####

#设置路径
setwd("F:/1/5LUSC/2TCGA-LUAD")
library(tidyverse)#读入meta.data文件
json <- jsonlite::fromJSON("metadata.cart.2023-06-13.json")
#查看
#View(json)
#获取需要的数据
sample_id <- sapply(json$associated_entities,function(x){x[,1]})
file_sample <- data.frame(sample_id,file_name=json$file_name)  
#读入counts文件
count_file <- list.files('gdc_download_20230613_041808.500125/',
                         pattern = '*.tsv',recursive = TRUE)
#整理
count_file_name <- strsplit(count_file,split='/')
count_file_name <- sapply(count_file_name,function(x){x[2]})
#必要时修改下面的基因数
matrix = data.frame(matrix(nrow=60660,ncol=0))
#下面的修改样本例数
for (i in 1:557){
  path = paste0('gdc_download_20230613_041808.500125//',count_file[i])   #Counts文件夹名
  data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
  colnames(data)<-data[2,]
  data <-data[-c(1:6),]
  #数据类型，选择其中之一 3：unstranded；4：stranded_first；5：stranded_second；6：tpm_unstranded；7：fpkm_unstranded；8：fpkm_uq_unstranded
  #data <- data[3]
  data <- data[6]
  #data <- data[7]
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix <- cbind(matrix,data)
}
#选择3，即为counts
#write.csv(matrix,'1TCGA/LUAD_Counts.csv',row.names = TRUE)
#选择6，即为TPM
#write.csv(matrix,'TPM.csv',row.names = TRUE)
#选择7，即为FPKM
#write.csv(matrix,'1TCGA/LUAD_FPKM.csv',row.names = TRUE)

#转化为gene_symbol
path = paste0('gdc_download_20230613_041808.500125//',count_file[1])
data<- as.matrix(read.delim(path,fill = TRUE,header = FALSE,row.names = 1))
gene_name <-data[-c(1:6),1]
matrix0 <- cbind(gene_name,matrix)
gene_type <- data[-c(1:6),2]
matrix0 <- cbind(gene_type,matrix0)
#将gene_name列去除重复的基因，保留每个基因最大表达量结果
matrix0 <- aggregate( . ~ gene_name,data=matrix0, max)    
#保留mRNA
matrix0 <- subset(x = matrix0, gene_type == "protein_coding")
#将gene_name列设为行名
rownames(matrix0) <- matrix0[,1]
matrix0 <- matrix0[,-c(1,2)]
#write.csv(matrix0,'1TCGA/LUAD_Counts_GeneSymbol.csv',row.names = TRUE)
#write.csv(matrix0,'1TCGA/LUAD_TPM_GeneSymbol.txt',row.names = TRUE)
#write.csv(matrix0,'1TCGA/LUAD_FPKM_GeneSymbol.csv',row.names = TRUE)
#write.table(matrix0,'TCGA.txt', sep="\t", quote=F, row.names = TRUE)

matrix1 = data.frame(ID=rownames(matrix0),matrix0)
#colnames(matrix1)
name123 = gsub('[.]', '-', colnames(matrix1))
colnames(matrix1) = name123
#write.table(matrix1,'TCGA_TPM.txt', sep="\t", quote=F, row.names = F)
write.table(matrix1,'TCGA_LUAD.txt', sep="\t", quote=F, row.names = F)

####3.1GSE30219-GPL570####

setwd("F:/1/5LUSC/3GEO/GSE30219-GPL570")
#下载矩阵
gset <- getGEO("GSE30219",destdir = "F:/1/5LUSC/3GEO/GSE30219-GPL570",AnnotGPL = F,getGPL = F) 
gset
#提取gset中第一个元素（包含基因表达矩阵和注释信息），并赋值给变量a
a=gset[[1]]
#提取a中的表达矩阵并赋值给dat，这时的表达矩阵的基因名字为探针ID
dat=exprs(a)
#展示表达矩阵的前6行，看下数据是否经过log转换，一般数据在20以内，是经过log转换的，若有成百上千的数据，表示这个数据还需要进一步log处理。
head(dat)
#也可以输入上面这段代码，如果数据已经经log后，显示log2 transform not needed；如果数据没有尽行log，需要log程序按照代码尽行log转换，完成后显示log2 transform finished。
ex <- dat
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NA
dat <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
#查看a的临床信息，为后面选择用于分组的变量做准备
pd=pData(a)
#查看pd
#View(pd)
#输出临床信息
write.csv(pd,'clinicalGSE30219.csv',row.names = TRUE)

#另外一种获取GPL文件方法
#如果没有下载过这个包，就下载一下，代码如下
#BiocManager::install("hgu133plus2.db")
#提取信息
library(hgu133plus2.db)
ls("package:hgu133plus2.db") #大致查看R包中的信息，寻找我们需要的SYMBOL
ids <- toTable(hgu133plus2SYMBOL)  #用toTable（）函数提取
head(ids) #查看提取内容
library(dplyr)
colnames(ids) = c("probe_id" ,"symbol")
#exp <- mutate(ex,probe_id=rownames(probe_exp))#将行名变为列名为probe_id的一列
# exp是原来的表达矩阵
#exp2= merge(exp,ids,by.x="probe_id", by.y="probe_id") # 合并数据
#exp2=exp2[!duplicated(probe_exp$symbol),]# 按照symbol列去重
#rownames(exp2)=exp2$symbol # 数据框probe_exp的行名变成symbol
#gene_exp_matrix <- na.omit(exp2 ) # 去空值
#输出文件
#write.table(gene_exp_matrix2,file = "Gastric.cancer.geneid.exprs.txt",sep = "\t",row.names=T,col.names = T)

#另外一种方式
#下载注释文件
#gpl <- getGEO('GPL10558', destdir="")
#我们将平台文件转为list格式，并赋值给gpl1，将gpl1保存为R文件，方便以后调用。
#gpl1<-Table(gpl)
#查看平台文件的列名，我们看到有ID和gene symbol,记住gene symbol在第几列，我们这里在第11列。
#colnames(Table(gpl))
#再了解一下平台文件的数据，当然这里大家可以直接选择gene symbol字段也行，主要了解一下symbol中基因symbol值是否只有一个。
#View(gpl1)
#我们这里的gene symbol字段中的symbol，有的基因就不止一个名称，///后面有重名，我们需要第一个名字，所以需要用字符分割包再处理一下,原理同上面处理title一样。
#提取平台文件gpl1中的ID和gene symbol，并赋值给probe2symbol_df
#gpl1$`Gene Symbol`=str_split(gpl1[,11],'///',simplify = T)[,1]
#probe2symbol_df<-gpl1[,c(1,11)]

probe2symbol_df <- ids
#将列名改为probe_id和symbol
#这两步是因为我懒，不想再调整代码，为了方便后面代码运行，我改的，大家也可以不改。
colnames(probe2symbol_df)=c('probe_id','symbol')
#查看symbol为唯一值的个数
length(unique(probe2symbol_df$symbol))
#查看每个基因出现n次的个数,我们可以看到，symbol出现一次的基因有7050个，出现2次有1493个。。。
table(sort(table(probe2symbol_df$symbol)))
#去掉没有注释symbol的探针（其实这里没有注释的探针数量即为上面出现次数最多的基因440，也就是说有440个探针没有symbol）
ids=probe2symbol_df[probe2symbol_df$symbol != '',]
#%in%用于判断是否匹配，
#注意这里dat是gset表达矩阵的数据，这一步就是把平台文件的ID和矩阵中的ID匹配。
ids=probe2symbol_df[probe2symbol_df$probe_id %in%rownames(dat),]
#取表达矩阵中可以与探针名匹配的那些，去掉无法匹配的表达数据
dat=dat[ids$probe_id,]
#ids新建mean这一列，列名为mean，同时对dat这个矩阵按行操作，取每一行（即每个样本在这个探针上的）的均值，将结果给到mean这一列的每一行，这里也可以改为中位值，median.
ids$mean=apply(dat,1,mean)
#即先按symbol排序，相同的symbol再按照均值从大到小排列，方便后续保留第一个值。
ids=ids[order(ids$symbol,ids$mean,decreasing = T),]
#将symbol这一列取出重复项，'!'为否，即取出不重复的项，去除重复的gene
#取median的话，这样就保留每个基因最大表达量结果.最后得到n个基因。
ids=ids[!duplicated(ids$symbol),]
#新的ids取出probe_id这一列，将dat按照取出的这一列，用他们的每一行组成一个新的dat
dat=dat[ids$probe_id,]
#把ids的symbol这一列中的每一行给dat作为dat的行名
rownames(dat)=ids$symbol
View(dat)
#但是我们看到矩阵最后一行，是没有symbol名字的，我们把他去掉，数字自己更改。
dim(dat)
#dat<-dat[-13437,]
#rownames(dat)[20911] = "A1BG"
#最后我们把结果保存。
#write.csv(dat,file="2GEO/GSE68465/GSE68465.csv")
#write.table(dat,file="GEO.txt", sep="\t", quote=F, row.names = TRUE)
write.table(data.frame(ID=rownames(dat),dat),file="GSE30219.txt", sep="\t", quote=F, row.names = F,col.names = T)




####3.1GSE37745-GPL570####

setwd("F:/1/5LUSC/3GEO/GSE37745-GPL570")
#下载矩阵
gset <- getGEO("GSE37745",destdir = "F:/1/5LUSC/3GEO/GSE37745-GPL570",AnnotGPL = F,getGPL = F) 
gset
#提取gset中第一个元素（包含基因表达矩阵和注释信息），并赋值给变量a
a=gset[[1]]
#提取a中的表达矩阵并赋值给dat，这时的表达矩阵的基因名字为探针ID
dat=exprs(a)
#展示表达矩阵的前6行，看下数据是否经过log转换，一般数据在20以内，是经过log转换的，若有成百上千的数据，表示这个数据还需要进一步log处理。
head(dat)
#也可以输入上面这段代码，如果数据已经经log后，显示log2 transform not needed；如果数据没有尽行log，需要log程序按照代码尽行log转换，完成后显示log2 transform finished。
ex <- dat
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NA
dat <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
#查看a的临床信息，为后面选择用于分组的变量做准备
pd=pData(a)
#查看pd
#View(pd)
#输出临床信息
write.csv(pd,'clinicalGSE37745.csv',row.names = TRUE)

#另外一种获取GPL文件方法
#如果没有下载过这个包，就下载一下，代码如下
#BiocManager::install("hgu133plus2.db")
#提取信息
library(hgu133plus2.db)
ls("package:hgu133plus2.db") #大致查看R包中的信息，寻找我们需要的SYMBOL
ids <- toTable(hgu133plus2SYMBOL)  #用toTable（）函数提取
head(ids) #查看提取内容
library(dplyr)
colnames(ids) = c("probe_id" ,"symbol")
#exp <- mutate(ex,probe_id=rownames(probe_exp))#将行名变为列名为probe_id的一列
# exp是原来的表达矩阵
#exp2= merge(exp,ids,by.x="probe_id", by.y="probe_id") # 合并数据
#exp2=exp2[!duplicated(probe_exp$symbol),]# 按照symbol列去重
#rownames(exp2)=exp2$symbol # 数据框probe_exp的行名变成symbol
#gene_exp_matrix <- na.omit(exp2 ) # 去空值
#输出文件
#write.table(gene_exp_matrix2,file = "Gastric.cancer.geneid.exprs.txt",sep = "\t",row.names=T,col.names = T)

#另外一种方式
#下载注释文件
#gpl <- getGEO('GPL10558', destdir="")
#我们将平台文件转为list格式，并赋值给gpl1，将gpl1保存为R文件，方便以后调用。
#gpl1<-Table(gpl)
#查看平台文件的列名，我们看到有ID和gene symbol,记住gene symbol在第几列，我们这里在第11列。
#colnames(Table(gpl))
#再了解一下平台文件的数据，当然这里大家可以直接选择gene symbol字段也行，主要了解一下symbol中基因symbol值是否只有一个。
#View(gpl1)
#我们这里的gene symbol字段中的symbol，有的基因就不止一个名称，///后面有重名，我们需要第一个名字，所以需要用字符分割包再处理一下,原理同上面处理title一样。
#提取平台文件gpl1中的ID和gene symbol，并赋值给probe2symbol_df
#gpl1$`Gene Symbol`=str_split(gpl1[,11],'///',simplify = T)[,1]
#probe2symbol_df<-gpl1[,c(1,11)]

probe2symbol_df <- ids
#将列名改为probe_id和symbol
#这两步是因为我懒，不想再调整代码，为了方便后面代码运行，我改的，大家也可以不改。
colnames(probe2symbol_df)=c('probe_id','symbol')
#查看symbol为唯一值的个数
length(unique(probe2symbol_df$symbol))
#查看每个基因出现n次的个数,我们可以看到，symbol出现一次的基因有7050个，出现2次有1493个。。。
table(sort(table(probe2symbol_df$symbol)))
#去掉没有注释symbol的探针（其实这里没有注释的探针数量即为上面出现次数最多的基因440，也就是说有440个探针没有symbol）
ids=probe2symbol_df[probe2symbol_df$symbol != '',]
#%in%用于判断是否匹配，
#注意这里dat是gset表达矩阵的数据，这一步就是把平台文件的ID和矩阵中的ID匹配。
ids=probe2symbol_df[probe2symbol_df$probe_id %in%rownames(dat),]
#取表达矩阵中可以与探针名匹配的那些，去掉无法匹配的表达数据
dat=dat[ids$probe_id,]
#ids新建mean这一列，列名为mean，同时对dat这个矩阵按行操作，取每一行（即每个样本在这个探针上的）的均值，将结果给到mean这一列的每一行，这里也可以改为中位值，median.
ids$mean=apply(dat,1,mean)
#即先按symbol排序，相同的symbol再按照均值从大到小排列，方便后续保留第一个值。
ids=ids[order(ids$symbol,ids$mean,decreasing = T),]
#将symbol这一列取出重复项，'!'为否，即取出不重复的项，去除重复的gene
#取median的话，这样就保留每个基因最大表达量结果.最后得到n个基因。
ids=ids[!duplicated(ids$symbol),]
#新的ids取出probe_id这一列，将dat按照取出的这一列，用他们的每一行组成一个新的dat
dat=dat[ids$probe_id,]
#把ids的symbol这一列中的每一行给dat作为dat的行名
rownames(dat)=ids$symbol
#View(dat)
#但是我们看到矩阵最后一行，是没有symbol名字的，我们把他去掉，数字自己更改。
dim(dat)
#dat<-dat[-13437,]
#rownames(dat)[20911] = "A1BG"
#最后我们把结果保存。
#write.csv(dat,file="2GEO/GSE68465/GSE68465.csv")
#write.table(dat,file="GEO.txt", sep="\t", quote=F, row.names = TRUE)
write.table(data.frame(ID=rownames(dat),dat),file="GSE37745.txt", sep="\t", quote=F, row.names = F,col.names = T)




####3.1GSE41271-GPL6884####

setwd("F:/1/5LUSC/3GEO/GSE41271-GPL6884")
#下载矩阵
gset <- getGEO("GSE41271",destdir = "F:/1/5LUSC/3GEO/GSE41271-GPL6884",AnnotGPL = F,getGPL = F) 
gset
#提取gset中第一个元素（包含基因表达矩阵和注释信息），并赋值给变量a
a=gset[[1]]
#提取a中的表达矩阵并赋值给dat，这时的表达矩阵的基因名字为探针ID
dat=exprs(a)
#展示表达矩阵的前6行，看下数据是否经过log转换，一般数据在20以内，是经过log转换的，若有成百上千的数据，表示这个数据还需要进一步log处理。
head(dat)
#也可以输入上面这段代码，如果数据已经经log后，显示log2 transform not needed；如果数据没有尽行log，需要log程序按照代码尽行log转换，完成后显示log2 transform finished。
ex <- dat
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NA
dat <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
#查看a的临床信息，为后面选择用于分组的变量做准备
pd=pData(a)
#查看pd
#View(pd)
#输出临床信息
write.csv(pd,'clinicalGSE41271.csv',row.names = TRUE)

#另外一种获取GPL文件方法
#如果没有下载过这个包，就下载一下，代码如下
#BiocManager::install("illuminaHumanv4.db")
#提取信息
#library(illuminaHumanv4.db)
#ls("package:illuminaHumanv4.db") #大致查看R包中的信息，寻找我们需要的SYMBOL
#ids <- toTable(illuminaHumanv4SYMBOL)  #用toTable（）函数提取
#head(ids) #查看提取内容
#library(dplyr)
#colnames(ids) = c("probe_id" ,"symbol")
#exp <- mutate(ex,probe_id=rownames(probe_exp))#将行名变为列名为probe_id的一列
# exp是原来的表达矩阵
#exp2= merge(exp,ids,by.x="probe_id", by.y="probe_id") # 合并数据
#exp2=exp2[!duplicated(probe_exp$symbol),]# 按照symbol列去重
#rownames(exp2)=exp2$symbol # 数据框probe_exp的行名变成symbol
#gene_exp_matrix <- na.omit(exp2 ) # 去空值
#输出文件
#write.table(gene_exp_matrix2,file = "Gastric.cancer.geneid.exprs.txt",sep = "\t",row.names=T,col.names = T)


library(idmap3)
ls('package:idmap3')
idmap3::get_pipe_IDs('GPL6884')



ids <- idmap3::get_pipe_IDs('GPL6884')  #用toTable（）函数提取
head(ids) #查看提取内容
library(dplyr)
colnames(ids) = c("probe_id" ,"symbol")


#另外一种方式
#下载注释文件
#gpl <- getGEO('GPL10558', destdir="")
#我们将平台文件转为list格式，并赋值给gpl1，将gpl1保存为R文件，方便以后调用。
#gpl1<-Table(gpl)
#查看平台文件的列名，我们看到有ID和gene symbol,记住gene symbol在第几列，我们这里在第11列。
#colnames(Table(gpl))
#再了解一下平台文件的数据，当然这里大家可以直接选择gene symbol字段也行，主要了解一下symbol中基因symbol值是否只有一个。
#View(gpl1)
#我们这里的gene symbol字段中的symbol，有的基因就不止一个名称，///后面有重名，我们需要第一个名字，所以需要用字符分割包再处理一下,原理同上面处理title一样。
#提取平台文件gpl1中的ID和gene symbol，并赋值给probe2symbol_df
#gpl1$`Gene Symbol`=str_split(gpl1[,11],'///',simplify = T)[,1]
#probe2symbol_df<-gpl1[,c(1,11)]

probe2symbol_df <- ids
#将列名改为probe_id和symbol
#这两步是因为我懒，不想再调整代码，为了方便后面代码运行，我改的，大家也可以不改。
colnames(probe2symbol_df)=c('probe_id','symbol')
#查看symbol为唯一值的个数
length(unique(probe2symbol_df$symbol))
#查看每个基因出现n次的个数,我们可以看到，symbol出现一次的基因有7050个，出现2次有1493个。。。
table(sort(table(probe2symbol_df$symbol)))
#去掉没有注释symbol的探针（其实这里没有注释的探针数量即为上面出现次数最多的基因440，也就是说有440个探针没有symbol）
ids=probe2symbol_df[probe2symbol_df$symbol != '',]
#%in%用于判断是否匹配，
#注意这里dat是gset表达矩阵的数据，这一步就是把平台文件的ID和矩阵中的ID匹配。
ids=probe2symbol_df[probe2symbol_df$probe_id %in%rownames(dat),]
#取表达矩阵中可以与探针名匹配的那些，去掉无法匹配的表达数据
dat=dat[ids$probe_id,]
#ids新建mean这一列，列名为mean，同时对dat这个矩阵按行操作，取每一行（即每个样本在这个探针上的）的均值，将结果给到mean这一列的每一行，这里也可以改为中位值，median.
ids$mean=apply(dat,1,mean)
#即先按symbol排序，相同的symbol再按照均值从大到小排列，方便后续保留第一个值。
ids=ids[order(ids$symbol,ids$mean,decreasing = T),]
#将symbol这一列取出重复项，'!'为否，即取出不重复的项，去除重复的gene
#取median的话，这样就保留每个基因最大表达量结果.最后得到n个基因。
ids=ids[!duplicated(ids$symbol),]
#新的ids取出probe_id这一列，将dat按照取出的这一列，用他们的每一行组成一个新的dat
dat=dat[ids$probe_id,]
#把ids的symbol这一列中的每一行给dat作为dat的行名
rownames(dat)=ids$symbol
View(dat)
#但是我们看到矩阵最后一行，是没有symbol名字的，我们把他去掉，数字自己更改。
dim(dat)
#dat<-dat[-13437,]
#rownames(dat)[20911] = "A1BG"
#最后我们把结果保存。
#write.csv(dat,file="2GEO/GSE68465/GSE68465.csv")
#write.table(dat,file="GEO.txt", sep="\t", quote=F, row.names = TRUE)
write.table(data.frame(ID=rownames(dat),dat),file="GSE41271.txt", sep="\t", quote=F, row.names = F)




####3.1GSE68465-GPL96,只有1W####

setwd("F:/1/5LUSC/3GEO/GSE68465-GPL96")
#下载矩阵
gset <- getGEO("GSE68465",destdir = "F:/1/5LUSC/3GEO/GSE68465-GPL96",AnnotGPL = F,getGPL = F) 
gset
#提取gset中第一个元素（包含基因表达矩阵和注释信息），并赋值给变量a
a=gset[[1]]
#提取a中的表达矩阵并赋值给dat，这时的表达矩阵的基因名字为探针ID
dat=exprs(a)
#展示表达矩阵的前6行，看下数据是否经过log转换，一般数据在20以内，是经过log转换的，若有成百上千的数据，表示这个数据还需要进一步log处理。
head(dat)
#也可以输入上面这段代码，如果数据已经经log后，显示log2 transform not needed；如果数据没有尽行log，需要log程序按照代码尽行log转换，完成后显示log2 transform finished。
ex <- dat
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NA
dat <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
#查看a的临床信息，为后面选择用于分组的变量做准备
pd=pData(a)
#查看pd
#View(pd)
#输出临床信息
write.csv(pd,'clinicalGSE68465.csv',row.names = TRUE)

#另外一种获取GPL文件方法
#如果没有下载过这个包，就下载一下，代码如下
#BiocManager::install("hgu133a.db")
#提取信息
library(hgu133a.db)
ls("package:hgu133a.db") #大致查看R包中的信息，寻找我们需要的SYMBOL
ids <- toTable(hgu133aSYMBOL)  #用toTable（）函数提取
head(ids) #查看提取内容
library(dplyr)
colnames(ids) = c("probe_id" ,"symbol")
#exp <- mutate(ex,probe_id=rownames(probe_exp))#将行名变为列名为probe_id的一列
# exp是原来的表达矩阵
#exp2= merge(exp,ids,by.x="probe_id", by.y="probe_id") # 合并数据
#exp2=exp2[!duplicated(probe_exp$symbol),]# 按照symbol列去重
#rownames(exp2)=exp2$symbol # 数据框probe_exp的行名变成symbol
#gene_exp_matrix <- na.omit(exp2 ) # 去空值
#输出文件
#write.table(gene_exp_matrix2,file = "Gastric.cancer.geneid.exprs.txt",sep = "\t",row.names=T,col.names = T)

#另外一种方式
#下载注释文件
#gpl <- getGEO('GPL10558', destdir="")
#我们将平台文件转为list格式，并赋值给gpl1，将gpl1保存为R文件，方便以后调用。
#gpl1<-Table(gpl)
#查看平台文件的列名，我们看到有ID和gene symbol,记住gene symbol在第几列，我们这里在第11列。
#colnames(Table(gpl))
#再了解一下平台文件的数据，当然这里大家可以直接选择gene symbol字段也行，主要了解一下symbol中基因symbol值是否只有一个。
#View(gpl1)
#我们这里的gene symbol字段中的symbol，有的基因就不止一个名称，///后面有重名，我们需要第一个名字，所以需要用字符分割包再处理一下,原理同上面处理title一样。
#提取平台文件gpl1中的ID和gene symbol，并赋值给probe2symbol_df
#gpl1$`Gene Symbol`=str_split(gpl1[,11],'///',simplify = T)[,1]
#probe2symbol_df<-gpl1[,c(1,11)]

probe2symbol_df <- ids
#将列名改为probe_id和symbol
#这两步是因为我懒，不想再调整代码，为了方便后面代码运行，我改的，大家也可以不改。
colnames(probe2symbol_df)=c('probe_id','symbol')
#查看symbol为唯一值的个数
length(unique(probe2symbol_df$symbol))
#查看每个基因出现n次的个数,我们可以看到，symbol出现一次的基因有7050个，出现2次有1493个。。。
table(sort(table(probe2symbol_df$symbol)))
#去掉没有注释symbol的探针（其实这里没有注释的探针数量即为上面出现次数最多的基因440，也就是说有440个探针没有symbol）
ids=probe2symbol_df[probe2symbol_df$symbol != '',]
#%in%用于判断是否匹配，
#注意这里dat是gset表达矩阵的数据，这一步就是把平台文件的ID和矩阵中的ID匹配。
ids=probe2symbol_df[probe2symbol_df$probe_id %in%rownames(dat),]
#取表达矩阵中可以与探针名匹配的那些，去掉无法匹配的表达数据
dat=dat[ids$probe_id,]
#ids新建mean这一列，列名为mean，同时对dat这个矩阵按行操作，取每一行（即每个样本在这个探针上的）的均值，将结果给到mean这一列的每一行，这里也可以改为中位值，median.
ids$mean=apply(dat,1,mean)
#即先按symbol排序，相同的symbol再按照均值从大到小排列，方便后续保留第一个值。
ids=ids[order(ids$symbol,ids$mean,decreasing = T),]
#将symbol这一列取出重复项，'!'为否，即取出不重复的项，去除重复的gene
#取median的话，这样就保留每个基因最大表达量结果.最后得到n个基因。
ids=ids[!duplicated(ids$symbol),]
#新的ids取出probe_id这一列，将dat按照取出的这一列，用他们的每一行组成一个新的dat
dat=dat[ids$probe_id,]
#把ids的symbol这一列中的每一行给dat作为dat的行名
rownames(dat)=ids$symbol
#View(dat)
#但是我们看到矩阵最后一行，是没有symbol名字的，我们把他去掉，数字自己更改。
dim(dat)
#dat<-dat[-13437,]
#rownames(dat)[20911] = "A1BG"
#最后我们把结果保存。
#write.csv(dat,file="2GEO/GSE68465/GSE68465.csv")
#write.table(dat,file="GEO.txt", sep="\t", quote=F, row.names = TRUE)
write.table(data.frame(ID=rownames(dat),dat),file="GSE68465.txt", sep="\t", quote=F, row.names = F,col.names = T)




####3.1GSE73403-GPL6480####

setwd("F:/1/5LUSC/3GEO/GSE73403-GPL6480")
#下载矩阵
gset <- getGEO("GSE73403",destdir = "F:/1/5LUSC/3GEO/GSE73403-GPL6480",AnnotGPL = F,getGPL = F) 
gset
#提取gset中第一个元素（包含基因表达矩阵和注释信息），并赋值给变量a
a=gset[[1]]
#提取a中的表达矩阵并赋值给dat，这时的表达矩阵的基因名字为探针ID
dat=exprs(a)
#展示表达矩阵的前6行，看下数据是否经过log转换，一般数据在20以内，是经过log转换的，若有成百上千的数据，表示这个数据还需要进一步log处理。
head(dat)
#也可以输入上面这段代码，如果数据已经经log后，显示log2 transform not needed；如果数据没有尽行log，需要log程序按照代码尽行log转换，完成后显示log2 transform finished。
ex <- dat
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NA
dat <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
#查看a的临床信息，为后面选择用于分组的变量做准备
pd=pData(a)
#查看pd
#View(pd)
#输出临床信息
write.csv(pd,'clinicalGSE73403.csv',row.names = TRUE)

#另外一种获取GPL文件方法
#如果没有下载过这个包，就下载一下，代码如下
#BiocManager::install("illuminaHumanv4.db")
#提取信息
#library(illuminaHumanv4.db)
#ls("package:illuminaHumanv4.db") #大致查看R包中的信息，寻找我们需要的SYMBOL
#ids <- toTable(illuminaHumanv4SYMBOL)  #用toTable（）函数提取
#head(ids) #查看提取内容
#library(dplyr)
#colnames(ids) = c("probe_id" ,"symbol")
#exp <- mutate(ex,probe_id=rownames(probe_exp))#将行名变为列名为probe_id的一列
# exp是原来的表达矩阵
#exp2= merge(exp,ids,by.x="probe_id", by.y="probe_id") # 合并数据
#exp2=exp2[!duplicated(probe_exp$symbol),]# 按照symbol列去重
#rownames(exp2)=exp2$symbol # 数据框probe_exp的行名变成symbol
#gene_exp_matrix <- na.omit(exp2 ) # 去空值
#输出文件
#write.table(gene_exp_matrix2,file = "Gastric.cancer.geneid.exprs.txt",sep = "\t",row.names=T,col.names = T)


library(idmap3)
ls('package:idmap3')
idmap3::get_pipe_IDs('GPL6480')



ids <- idmap3::get_pipe_IDs('GPL6480')  #用toTable（）函数提取
head(ids) #查看提取内容
library(dplyr)
colnames(ids) = c("probe_id" ,"symbol")


#另外一种方式
#下载注释文件
#gpl <- getGEO('GPL10558', destdir="")
#我们将平台文件转为list格式，并赋值给gpl1，将gpl1保存为R文件，方便以后调用。
#gpl1<-Table(gpl)
#查看平台文件的列名，我们看到有ID和gene symbol,记住gene symbol在第几列，我们这里在第11列。
#colnames(Table(gpl))
#再了解一下平台文件的数据，当然这里大家可以直接选择gene symbol字段也行，主要了解一下symbol中基因symbol值是否只有一个。
#View(gpl1)
#我们这里的gene symbol字段中的symbol，有的基因就不止一个名称，///后面有重名，我们需要第一个名字，所以需要用字符分割包再处理一下,原理同上面处理title一样。
#提取平台文件gpl1中的ID和gene symbol，并赋值给probe2symbol_df
#gpl1$`Gene Symbol`=str_split(gpl1[,11],'///',simplify = T)[,1]
#probe2symbol_df<-gpl1[,c(1,11)]

probe2symbol_df <- ids
#将列名改为probe_id和symbol
#这两步是因为我懒，不想再调整代码，为了方便后面代码运行，我改的，大家也可以不改。
colnames(probe2symbol_df)=c('probe_id','symbol')
#查看symbol为唯一值的个数
length(unique(probe2symbol_df$symbol))
#查看每个基因出现n次的个数,我们可以看到，symbol出现一次的基因有7050个，出现2次有1493个。。。
table(sort(table(probe2symbol_df$symbol)))
#去掉没有注释symbol的探针（其实这里没有注释的探针数量即为上面出现次数最多的基因440，也就是说有440个探针没有symbol）
ids=probe2symbol_df[probe2symbol_df$symbol != '',]
#%in%用于判断是否匹配，
#注意这里dat是gset表达矩阵的数据，这一步就是把平台文件的ID和矩阵中的ID匹配。
ids=probe2symbol_df[probe2symbol_df$probe_id %in%rownames(dat),]
#取表达矩阵中可以与探针名匹配的那些，去掉无法匹配的表达数据
dat=dat[ids$probe_id,]
#ids新建mean这一列，列名为mean，同时对dat这个矩阵按行操作，取每一行（即每个样本在这个探针上的）的均值，将结果给到mean这一列的每一行，这里也可以改为中位值，median.
ids$mean=apply(dat,1,mean)
#即先按symbol排序，相同的symbol再按照均值从大到小排列，方便后续保留第一个值。
ids=ids[order(ids$symbol,ids$mean,decreasing = T),]
#将symbol这一列取出重复项，'!'为否，即取出不重复的项，去除重复的gene
#取median的话，这样就保留每个基因最大表达量结果.最后得到n个基因。
ids=ids[!duplicated(ids$symbol),]
#新的ids取出probe_id这一列，将dat按照取出的这一列，用他们的每一行组成一个新的dat
dat=dat[ids$probe_id,]
#把ids的symbol这一列中的每一行给dat作为dat的行名
rownames(dat)=ids$symbol
View(dat)
#但是我们看到矩阵最后一行，是没有symbol名字的，我们把他去掉，数字自己更改。
dim(dat)
#dat<-dat[-13437,]
#rownames(dat)[20911] = "A1BG"
#最后我们把结果保存。
#write.csv(dat,file="2GEO/GSE68465/GSE68465.csv")
#write.table(dat,file="GEO.txt", sep="\t", quote=F, row.names = TRUE)
write.table(data.frame(ID=rownames(dat),dat),file="GSE73403.txt", sep="\t", quote=F, row.names = F)


####3.1GSE74777-GPL17586####

setwd("F:/1/5LUSC/3GEO/GSE74777-GPL17586")
#下载矩阵,destdir填写后会自动识别
gset <- getGEO("GSE74777",destdir = "F:/1/5LUSC/3GEO/GSE74777-GPL17586",AnnotGPL = F,getGPL = F) 
gset
#提取gset中第一个元素（包含基因表达矩阵和注释信息），并赋值给变量a
a=gset[[1]]
#提取a中的表达矩阵并赋值给dat，这时的表达矩阵的基因名字为探针ID
dat=exprs(a)
#展示表达矩阵的前6行，看下数据是否经过log转换，一般数据在20以内，是经过log转换的，若有成百上千的数据，表示这个数据还需要进一步log处理。
head(dat)
#也可以输入上面这段代码，如果数据已经经log后，显示log2 transform not needed；如果数据没有尽行log，需要log程序按照代码尽行log转换，完成后显示log2 transform finished。
ex <- dat
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NA
dat <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
#查看a的临床信息，为后面选择用于分组的变量做准备
pd=pData(a)
#查看pd
#View(pd)
#输出临床信息
write.csv(pd,'clinical_GSE74777.csv',row.names = TRUE)


#另外一种方式
#下载注释文件
#gpl <- getGEO('GPL17586', destdir="F:/1/5LUSC/3GEO/GSE74777-GPL17586")
#df2=read.table("GPL25318_family.soft",sep = "\t")
#GPL=getGEO(filename = 'GPL25318_family.soft.gz') 
#gpl=read.table("GPL17586-45144.txt", header=T, sep="\t", check.names=F)


options(stringsAsFactors = F)
gpl=read.table("GPL17586-45144.txt",
               header = TRUE,fill = T,sep = "\t",
               comment.char = "#",
               stringsAsFactors = FALSE,
               quote = "")
head(gpl)
colnames(gpl)

ids=gpl[,c("ID","gene_assignment")]
head(ids)
colnames(ids)=c('probe_id','symbol')
head(ids)
probe2gene <- ids

library(stringr)  
probe2gene$symbol=trimws(str_split(probe2gene$symbol,'//',simplify = T)[,2])
plot(table(table(probe2gene$symbol)),xlim=c(1,50))
head(probe2gene)

#我们将平台文件转为list格式，并赋值给gpl1，将gpl1保存为R文件，方便以后调用。
#gpl1<-Table(df2)
#查看平台文件的列名，我们看到有ID和gene symbol,记住gene symbol在第几列，我们这里在第11列。
#colnames(Table(gpl))
#再了解一下平台文件的数据，当然这里大家可以直接选择gene symbol字段也行，主要了解一下symbol中基因symbol值是否只有一个。
#View(gpl1)
#我们这里的gene symbol字段中的symbol，有的基因就不止一个名称，
#后面有重名，我们需要第一个名字，所以需要用字符分割包再处理一下,原理同上面处理title一样。
#ids=gpl[,c(1,7)]
#提取平台文件gpl1中的ID和gene symbol，并赋值给probe2symbol_df
#gpl$`Gene Symbol`=str_split(gpl[,7],'///',simplify = T)[,1]
#probe2symbol_df<-gpl[,c(1,7)]
#ids <- probe2symbol_df
#probe2symbol_df <- ids
#将列名改为probe_id和symbol
#这两步是因为我懒，不想再调整代码，为了方便后面代码运行，我改的，大家也可以不改。
#colnames(probe2symbol_df)=c('probe_id','symbol')


probe2symbol_df=probe2gene
#查看symbol为唯一值的个数
length(unique(probe2symbol_df$symbol))
#查看每个基因出现n次的个数,我们可以看到，symbol出现一次的基因有7050个，出现2次有1493个。。。
table(sort(table(probe2symbol_df$symbol)))
#去掉没有注释symbol的探针（其实这里没有注释的探针数量即为上面出现次数最多的基因440，也就是说有440个探针没有symbol）
ids=probe2symbol_df[probe2symbol_df$symbol != '',]
#%in%用于判断是否匹配，
#注意这里dat是gset表达矩阵的数据，这一步就是把平台文件的ID和矩阵中的ID匹配。
ids=probe2symbol_df[probe2symbol_df$probe_id %in%rownames(dat),]
#取表达矩阵中可以与探针名匹配的那些，去掉无法匹配的表达数据
dat=dat[ids$probe_id,]
#ids新建mean这一列，列名为mean，同时对dat这个矩阵按行操作，取每一行（即每个样本在这个探针上的）的均值，将结果给到mean这一列的每一行，这里也可以改为中位值，median.
ids$mean=apply(dat,1,mean)
#即先按symbol排序，相同的symbol再按照均值从大到小排列，方便后续保留第一个值。
ids=ids[order(ids$symbol,ids$mean,decreasing = T),]
#将symbol这一列取出重复项，'!'为否，即取出不重复的项，去除重复的gene
#取median的话，这样就保留每个基因最大表达量结果.最后得到n个基因。
ids=ids[!duplicated(ids$symbol),]
#新的ids取出probe_id这一列，将dat按照取出的这一列，用他们的每一行组成一个新的dat
dat=dat[ids$probe_id,]
#把ids的symbol这一列中的每一行给dat作为dat的行名
rownames(dat)=ids$symbol
View(dat)
#但是我们看到矩阵最后一行，是没有symbol名字的，我们把他去掉，数字自己更改。
dim(dat)
dat<-dat[-30906,]
#最后我们把结果保存。
write.table(data.frame(ID=rownames(dat),dat),file="GSE74777.txt", sep="\t", quote=F, row.names = F)


####3.1GSE157010-GPL570####

setwd("F:/1/5LUSC/3GEO/GSE157010-GPL570")
#下载矩阵
gset <- getGEO("GSE157010",destdir = "F:/1/5LUSC/3GEO/GSE157010-GPL570",AnnotGPL = F,getGPL = F) 
gset
#提取gset中第一个元素（包含基因表达矩阵和注释信息），并赋值给变量a
a=gset[[1]]
#提取a中的表达矩阵并赋值给dat，这时的表达矩阵的基因名字为探针ID
dat=exprs(a)
#展示表达矩阵的前6行，看下数据是否经过log转换，一般数据在20以内，是经过log转换的，若有成百上千的数据，表示这个数据还需要进一步log处理。
head(dat)
#也可以输入上面这段代码，如果数据已经经log后，显示log2 transform not needed；如果数据没有尽行log，需要log程序按照代码尽行log转换，完成后显示log2 transform finished。
ex <- dat
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NA
dat <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
#查看a的临床信息，为后面选择用于分组的变量做准备
pd=pData(a)
#查看pd
#View(pd)
#输出临床信息
write.csv(pd,'clinicalGSE157010.csv',row.names = TRUE)

#另外一种获取GPL文件方法
#如果没有下载过这个包，就下载一下，代码如下
#BiocManager::install("hgu133plus2.db")
#提取信息
library(hgu133plus2.db)
ls("package:hgu133plus2.db") #大致查看R包中的信息，寻找我们需要的SYMBOL
ids <- toTable(hgu133plus2SYMBOL)  #用toTable（）函数提取
head(ids) #查看提取内容
library(dplyr)
colnames(ids) = c("probe_id" ,"symbol")
#exp <- mutate(ex,probe_id=rownames(probe_exp))#将行名变为列名为probe_id的一列
# exp是原来的表达矩阵
#exp2= merge(exp,ids,by.x="probe_id", by.y="probe_id") # 合并数据
#exp2=exp2[!duplicated(probe_exp$symbol),]# 按照symbol列去重
#rownames(exp2)=exp2$symbol # 数据框probe_exp的行名变成symbol
#gene_exp_matrix <- na.omit(exp2 ) # 去空值
#输出文件
#write.table(gene_exp_matrix2,file = "Gastric.cancer.geneid.exprs.txt",sep = "\t",row.names=T,col.names = T)

#另外一种方式
#下载注释文件
#gpl <- getGEO('GPL10558', destdir="")
#我们将平台文件转为list格式，并赋值给gpl1，将gpl1保存为R文件，方便以后调用。
#gpl1<-Table(gpl)
#查看平台文件的列名，我们看到有ID和gene symbol,记住gene symbol在第几列，我们这里在第11列。
#colnames(Table(gpl))
#再了解一下平台文件的数据，当然这里大家可以直接选择gene symbol字段也行，主要了解一下symbol中基因symbol值是否只有一个。
#View(gpl1)
#我们这里的gene symbol字段中的symbol，有的基因就不止一个名称，///后面有重名，我们需要第一个名字，所以需要用字符分割包再处理一下,原理同上面处理title一样。
#提取平台文件gpl1中的ID和gene symbol，并赋值给probe2symbol_df
#gpl1$`Gene Symbol`=str_split(gpl1[,11],'///',simplify = T)[,1]
#probe2symbol_df<-gpl1[,c(1,11)]

probe2symbol_df <- ids
#将列名改为probe_id和symbol
#这两步是因为我懒，不想再调整代码，为了方便后面代码运行，我改的，大家也可以不改。
colnames(probe2symbol_df)=c('probe_id','symbol')
#查看symbol为唯一值的个数
length(unique(probe2symbol_df$symbol))
#查看每个基因出现n次的个数,我们可以看到，symbol出现一次的基因有7050个，出现2次有1493个。。。
table(sort(table(probe2symbol_df$symbol)))
#去掉没有注释symbol的探针（其实这里没有注释的探针数量即为上面出现次数最多的基因440，也就是说有440个探针没有symbol）
ids=probe2symbol_df[probe2symbol_df$symbol != '',]
#%in%用于判断是否匹配，
#注意这里dat是gset表达矩阵的数据，这一步就是把平台文件的ID和矩阵中的ID匹配。
ids=probe2symbol_df[probe2symbol_df$probe_id %in%rownames(dat),]
#取表达矩阵中可以与探针名匹配的那些，去掉无法匹配的表达数据
dat=dat[ids$probe_id,]
#ids新建mean这一列，列名为mean，同时对dat这个矩阵按行操作，取每一行（即每个样本在这个探针上的）的均值，将结果给到mean这一列的每一行，这里也可以改为中位值，median.
ids$mean=apply(dat,1,mean)
#即先按symbol排序，相同的symbol再按照均值从大到小排列，方便后续保留第一个值。
ids=ids[order(ids$symbol,ids$mean,decreasing = T),]
#将symbol这一列取出重复项，'!'为否，即取出不重复的项，去除重复的gene
#取median的话，这样就保留每个基因最大表达量结果.最后得到n个基因。
ids=ids[!duplicated(ids$symbol),]
#新的ids取出probe_id这一列，将dat按照取出的这一列，用他们的每一行组成一个新的dat
dat=dat[ids$probe_id,]
#把ids的symbol这一列中的每一行给dat作为dat的行名
rownames(dat)=ids$symbol
#View(dat)
#但是我们看到矩阵最后一行，是没有symbol名字的，我们把他去掉，数字自己更改。
dim(dat)
#dat<-dat[-13437,]
#rownames(dat)[20911] = "A1BG"
#最后我们把结果保存。
#write.csv(dat,file="2GEO/GSE68465/GSE68465.csv")
#write.table(dat,file="GEO.txt", sep="\t", quote=F, row.names = TRUE)
write.table(data.frame(ID=rownames(dat),dat),file="GSE157010.txt", sep="\t", quote=F, row.names = F,col.names = T)






####4.1scRNA-seq####

setwd("F:/1/5LUSC/4scRNA")
#加载功能包
source("sc_function.R")

pbmc <-fread("GSE153935_TLDS_AllCells.txt", sep="\t")
pbmc[1:10,1:10]
pbmc <- pbmc %>% column_to_rownames("V1")
pbmc[1:10,1:10]
dim(pbmc)



#dim(pbmc1)
#dim(pbmc2)
#合并
#pbmc=rbind(pbmc1,pbmc2,pbmc3,pbmc4,pbmc5,pbmc6)
#dim(pbmc)
#去除重复
#pbmc=pbmc[!duplicated(pbmc$barcode),]
#dim(pbmc)
#pbmc <- pbmc %>% column_to_rownames("barcode")
#dim(pbmc)

#创建seurat对象
pbmc <- CreateSeuratObject(pbmc, min.cells = 3, min.features = 300)
#pbmc <- CreateSeuratObject(pbmc)
dim(pbmc)

#检测双细胞
pbmc <- findDoublets(data = pbmc, SingleRref = "ref_Human_all.RData")
table(pbmc@meta.data$DF.classifications)
#创建QC文件夹
dir.create("1ResData")
#保存seurat
saveRDS(pbmc, file = "1ResData/1pbmc_Seurat.rds")

####4.2细胞质控####

setwd("F:/1/5LUSC/4scRNA")
#读取seurat
pbmc <- readRDS("1ResData/1pbmc_Seurat.rds")
#创建QC文件夹
dir.create("2QC")
#计算质控指标
#计算细胞中线粒体基因比例
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#计算细胞中核糖体基因比例
pbmc[["percent.rb"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[LS]")
#计算红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(pbmc))
pbmc[["percent.HB"]]<-PercentageFeatureSet(pbmc, features=HB.genes) 

#查看质控指标
#设置绘图元素
theme.set2 = theme(axis.title.x=element_blank())
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB")
group = "orig.ident"
#质控前小提琴图
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(pbmc, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)    
ggsave("2QC/1vlnplot_before_qc.pdf", plot = violin, width = 14, height = 8) 

#设置质控指标
quantile(pbmc$nFeature_RNA, seq(0.01, 0.1, 0.01))
quantile(pbmc$nFeature_RNA, seq(0.9, 1, 0.01))
plots[[1]] + geom_hline(yintercept = 500) + geom_hline(yintercept = 4500)
quantile(pbmc$nCount_RNA, seq(0.01, 0.1, 0.01))
plots[[2]] + geom_hline(yintercept = 22000)
quantile(pbmc$percent.mt, seq(0.9, 1, 0.01))
plots[[3]] + geom_hline(yintercept = 20)
quantile(pbmc$percent.HB, seq(0.9, 1, 0.01))
plots[[4]] + geom_hline(yintercept = 1)

#质控
#设置质控标准
minGene=300
maxGene=10000
minUMI=600
pctMT=20
pctHB=1
dim(pbmc)

#数据质控并绘制小提琴图
pbmc <- subset(pbmc, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & nCount_RNA > minUMI &
                 percent.mt < pctMT & percent.HB < pctHB)
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(pbmc, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=1)    
ggsave("2QC/2vlnplot_after_qc.pdf", plot = violin, width = 20, height = 5) 
dim(pbmc)

#细胞周期评分
pbmc <- NormalizeData(pbmc)  #解决每个细胞测序深度不同的问题
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(pbmc))
s_genes <- cc.genes$s.genes    
s_genes <- CaseMatch(search=s_genes, match=rownames(pbmc))
pbmc <- CellCycleScoring(pbmc, g2m.features=g2m_genes, s.features=s_genes)
colnames(pbmc@meta.data)
table(pbmc$Phase)

#保存
saveRDS(pbmc, file = "1ResData/2pbmc_qc.rds")

####4.3降维聚类####

setwd("F:/1/5LUSC/4scRNA")
#读取seurat
pbmc <- readRDS("1ResData/2pbmc_qc.rds")
#标准流程
#SCT
pbmc <- SCTransform(pbmc, vars.to.regress = c("S.Score", "G2M.Score"))
#pbmc <- SCTransform(pbmc)
#PCA
pbmc <- RunPCA(pbmc, verbose = F)

#保存
saveRDS(pbmc, file = "1ResData/3pbmc_PCA.rds")
#读取seurat
pbmc <- readRDS("1ResData/3pbmc_PCA.rds")

#
ElbowPlot(pbmc, ndims = 50)
pcs = 1:30
table(pbmc@meta.data$orig.ident)


#去除批次
dir.create("3UMAP")
#harmony
pbmc <- RunHarmony(pbmc, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20)
#降维聚类
#pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = pcs) %>% FindClusters(resolution = 1)
#pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = pcs) %>% RunTSNE(reduction = "harmony", dims = pcs)
#确定合适的分辨率
#devtools::install_github("PaulingLiu/ROGUE")
library(ROGUE)

seq = seq(0.5,2,by=0.1)
pbmc <- FindNeighbors(pbmc,  dims = pcs) 
for (res in seq){
  pbmc = FindClusters(pbmc, resolution = res)
}

library(clustree)

p1 = clustree(pbmc,prefix = "SCT_snn_res.")+coord_flip()
p = p1+plot_layout(widths = c(3,1))
ggsave("3UMAP/SCT_sun_res.png", p, width = 30, height = 14)

#降维聚类
pbmc <- FindNeighbors(pbmc,  dims = pcs) %>% FindClusters(resolution = 1)
pbmc <- RunUMAP(pbmc,  dims = pcs) %>% RunTSNE(dims = pcs)
#画图
#p <- DimPlot(pbmc, reduction = "umap", label = T)
#ggsave("3Harmony/1Cluster_reduction.png", p, width = 7, height = 6)
#p = DimPlot(pbmc,reduction = "umap",label = F,group.by = "orig.ident")
#ggsave("3Harmony/2Cluster_reduction_sample.png", p, width = 7, height = 6)
#p = DimPlot(pbmc,reduction = "umap",label = F,group.by = "site")
#ggsave("3Harmony/3Cluster_reduction_site.png", p, width = 7, height = 6)

Biocols = c('#AB3282', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
            '#E95C59', '#E59CC4','#E5D2DD' , '#23452F', '#BD956A', '#8C549C', '#585658',
            '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
            '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
            '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
            '#968175',"DarkOrchid")

p <- DimPlot(pbmc, reduction = "umap", label = T)#,cols = Biocols)
ggsave("3UMAP/1Cluster_reduction.png", p, width = 7, height = 6)
p <- DimPlot(pbmc, reduction = "tsne", label = T)#,cols = Biocols)
ggsave("3UMAP/1Cluster_reductionTSNE.png", p, width = 7, height = 6)
p = DimPlot(pbmc,reduction = "umap",label = F,group.by = "orig.ident")
ggsave("3UMAP/2Cluster_reduction_sample.png", p, width = 7, height = 6)
#保存
saveRDS(pbmc, file = "1ResData/4pbmc_TSNE.rds")

####4.4细胞注释####

setwd("F:/1/5LUSC/4scRNA")
#读取
pbmc <- readRDS("1ResData/4pbmc_TSNE.rds")
dir.create("4Celltype")
#Marker基因标记

markers <- c( "ACTA2","FAP","PDGFRB","NOTCH3", #Fibroblast
             "PECAM1","VWF","ENG", #Endothelia
             "CMA1","TPSAB1","TPSB2", #Mast
             "AR","KRT19","KRT18","KRT8","TP63","KRT14","KRT5","EPCAM",#Epithelial
             "LYZ","CSF1R","CD14","HAVCR2",#Monolytic
             "CD163","CD68","FCGR3A", #macrophage
             "MS4A1","MZB1","IGHG1","MKI67","TOP2A","CD27", #B/plasma
             "PTPRC","CD7","CD2","CD3G","CD3E","CD3D","IL7R") #T

#featureplot
#将表达矩阵SCT转为RNA
DefaultAssay(pbmc) <- "RNA"
p <- FeaturePlot(pbmc, features = markers, ncol = 3)
ggsave("4Celltype/1Markers_featureplot.pdf", p, width = 15, height = 32)
#dotplot
p <- DotPlot(pbmc, features = markers) + RotatedAxis()
ggsave("4Celltype/2Markers_dotplot.pdf", p, width = 14, height = 6)
#vlnplot
p <- VlnPlot(pbmc, features = markers, stack = T, flip = T) + NoLegend()
ggsave("4Celltype/3Markers_vlnplot.pdf", p, width = 14, height = 6)

#Fibroblast  Endothelia  Mast  Epithelial  Monolytic Macrophage B/Plasma T
pbmc$celltype.main <- recode(pbmc@meta.data$seurat_clusters,
                             "0" = "T",
                             "1" = "Monolytic/Macrophage",
                             "2" = "T",
                             "3" = "T",
                             "4" = "T",
                             "5" = "B/Plasma",
                             "6" = "T",
                             "7" = "T",
                             "8" = "Fibroblast",
                             "9" = "Fibroblast",
                             "10" = "Mast",
                             "11" = "Endothelia",
                             "12" = "Epithelial",
                             "13" = "Epithelial",
                             "14" = "Monolytic/Macrophage",
                             "15" = "Fibroblast",
                             "16" = "Monolytic/Macrophage",
                             "17" = "B/Plasma",
                             "18" = "Epithelial",
                             "19" = "Fibroblast",
                             "20" = "Fibroblast",
                             "21" = "Fibroblast",
                             "22" = "B/Plasma",
                             "23" = "B/Plasma",
                             "24" = "Epithelial",
                             "25" = "Epithelial")

#pbmc = subset(pbmc,idents = c(15,16,17),invert = T)

#table(pbmc@active.ident)
#table(pbmc1@active.ident)

#查看各个细胞群的亚型
table(pbmc@meta.data$celltype.main)
Biocols = c('#AB3282', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
            '#E95C59', '#E59CC4','#E5D2DD' , '#23452F', '#BD956A', '#8C549C', '#585658',
            '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
            '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
            '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
            '#968175',"DarkOrchid")

#画图
DimPlot(pbmc, reduction = "tsne", label = T,group.by = "celltype.main",cols = Biocols)
DimPlot(pbmc, reduction = "umap", label = T,group.by = "DF.classifications",cols = Biocols)
#画图
p = DimPlot(pbmc, reduction = "umap", label = T,group.by = "celltype.main")#,cols = Biocols)
ggsave("4Celltype/11Cluster_umap.png", p, width = 7, height = 6)
p = DimPlot(pbmc, reduction = "tsne", label = T,group.by = "celltype.main")#,cols = Biocols)
ggsave("4Celltype/12Cluster_tsne.png", p, width = 7, height = 6)
#保存

markers <- c( "ACTA2","FAP","PDGFRB", #Fibroblast
              "VWF","ENG", #Endothelia
              "TPSAB1", #Mast
              "AR","KRT19","KRT18","KRT8","EPCAM",#Epithelial
              "LYZ","CSF1R","CD14","HAVCR2",#Monolytic
              "CD163","CD68","FCGR3A", #macrophage
              "MS4A1","MZB1","IGHG1", #B/plasma
              "PTPRC","CD7","CD2","CD3G","CD3E","CD3D","IL7R") #T

DotPlot(pbmc, features = markers,group.by = "celltype.main") + RotatedAxis()

#featureplot
#将表达矩阵SCT转为RNA
DefaultAssay(pbmc) <- "RNA"
p <- FeaturePlot(pbmc, features = markers, ncol = 3)#,cols = c('#E39A35','#AB3282'))
ggsave("4Celltype/20Markers_featureplot.pdf", p, width = 15, height = 32)
#dotplot
p <- DotPlot(pbmc, features = markers,group.by = "celltype.main") + RotatedAxis()
             #,cols = c('#E39A35','#AB3282'))
ggsave("4Celltype/21Markers_dotplot.pdf", p, width = 14, height = 5)
#vlnplot
p <- VlnPlot(pbmc, features = markers, stack = T, flip = T) + NoLegend()
             #,cols = Biocols)
ggsave("4Celltype/22Markers_vlnplot.pdf", p, width = 14, height = 6)
DefaultAssay(pbmc) <- "SCT"

saveRDS(pbmc, file = "1ResData/5pbmc_cluster.rds")

pbmc <- readRDS("1ResData/5pbmc_cluster.rds")

Idents(pbmc)<-pbmc$celltype.main
DefaultAssay(pbmc) <- "RNA"
markers<-FindAllMarkers(pbmc,only.pos = T,min.pct = 0.3,logfc.threshold = 0.3,min.diff.pct = 0.2)
DefaultAssay(pbmc) <- "SCT"
write.csv(markers,file="markers.csv")


#提取Fibroblast
pbmc_Fib <- subset(pbmc,celltype.main %in% "Fibroblast")
table(pbmc@meta.data$celltype.main)
dim(pbmc_Fib)
#保存
saveRDS(pbmc_Fib, file = "1ResData/6pbmc_Fib.rds")

####4.5Fibroblast再分群####

setwd("F:/1/5LUSC/4scRNA")
#读取
pbmc_Fib <- readRDS("1ResData/6pbmc_Fib.rds")
#标准流程
#SCT
pbmc_Fib <- SCTransform(pbmc_Fib, vars.to.regress = c("S.Score", "G2M.Score"))
#PCA
pbmc_Fib <- RunPCA(pbmc_Fib, verbose = F)

#
ElbowPlot(pbmc_Fib, ndims = 50)
pcs = 1:30
table(pbmc_Fib@meta.data$orig.ident)


#去除批次
dir.create("5Fib")
#harmony
#pbmc_Fib <- RunHarmony(pbmc_Fib, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20)
#降维聚类
#pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = pcs) %>% FindClusters(resolution = 1)
#pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = pcs) %>% RunTSNE(reduction = "harmony", dims = pcs)
#确定合适的分辨率
library("ROGUE")

seq = seq(0.1,2,by=0.1)
pbmc_Fib <- FindNeighbors(pbmc_Fib,  dims = pcs) 
for (res in seq){
  pbmc_Fib = FindClusters(pbmc_Fib, resolution = res)
}

library(clustree)

p1 = clustree(pbmc_Fib,prefix = "SCT_snn_res.")+coord_flip()
p = p1+plot_layout(widths = c(3,1))
ggsave("5Fib/SCT_sun_res.png", p, width = 30, height = 10)

#library(Seurat)
#library(clustree)
#pbmc_Fib <- FindClusters(pbmc_Fib, resolution = seq(0.1,1,by=0.1))
#clustree(pbmc_Fib)

#install.packages("ggraph", version = "2.1.0")

pbmc_Fib <- RunHarmony(pbmc_Fib, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20)

pbmc_Fib <- FindNeighbors(pbmc_Fib,  dims = pcs) %>% FindClusters(resolution = 0.6)
pbmc_Fib <- RunUMAP(pbmc_Fib,  dims = pcs) %>% RunTSNE(dims = pcs)

#画图
Biocols = c('#AB3282', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
            '#E95C59', '#E59CC4','#E5D2DD' , '#23452F', '#BD956A', '#8C549C', '#585658',
            '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
            '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
            '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
            '#968175',"DarkOrchid")


p <- DimPlot(pbmc_Fib, reduction = "tsne", label = T)#,cols = Biocols)
ggsave("5Fib/1FibCluster_reduction.png", p, width = 7, height = 6)

p <- DimPlot(pbmc_Fib, reduction = "umap", label = T)#,cols = Biocols)
ggsave("5Fib/1FibCluster_reduction_umap.png", p, width = 7, height = 6)
#保存
#保存

#去除批次
dir.create("6FibCelltype")

write.csv(rownames(pbmc_Fib),file="6FibCelltype/rownames(pbmc_Fib).csv")

#Marker基因标记
markers <- c("COL1A1","COL10A1","COL4A1","MMP3","IL13","TGFB1", #myCAF "IL4",
             "IL6","IL11","IL8","LIF","CXCL1","CXCL2",
             "CXCL10","CCL2","CCL8","C1QA","C1QB","C1QC","CXCL12","CFD",#iCAF "CXCL8","GM-CSF",
             "CXCL14")#apCAF "CXCL12","CXCL14","CFD"

#将表达矩阵SCT转为RNA
DefaultAssay(pbmc_Fib) <- "RNA"
#p <- FeaturePlot(pbmc_Fib, reduction = "tsne", features = markers, ncol = 4)
#ggsave("6FibCelltype/1Fib_Markers_featureplot.pdf", p, width = 14, height = 8)
#dotplot
p <- DotPlot(pbmc_Fib, features = markers) + RotatedAxis()
ggsave("6FibCelltype/2Fib_Markers_dotplot.pdf", p, width = 14, height = 6)
#vlnplot
p <- VlnPlot(pbmc_Fib, features = markers, stack = T, flip = T) + NoLegend()
ggsave("6FibCelltype/3Fib_Markers_vlnplot.pdf", p, width = 14, height = 6)
DefaultAssay(pbmc_Fib) <- "SCT"


p1 <- DimPlot(pbmc_Fib, reduction = "tsne", label = T)
p2 <- DimPlot(pbmc_Fib, reduction = "umap", label = T)
p = p1+p2
ggsave("6FibCelltype/4FibCluster_reduction.png", p, width = 10, height = 6)




#保存
saveRDS(pbmc_Fib, file = "1ResData/7pbmc_Fib.rds")

#myCAF,apCAF,iCAF
pbmc_Fib$fib.celltype.main <- recode(pbmc_Fib@meta.data$seurat_clusters,
                                  "0" = "myCAF",
                                  "1" = "apCAF",
                                  "2" = "myCAF",
                                  "3" = "iCAF",
                                  "4" = "iCAF",
                                  "5" = "myCAF",
                                  "6" = "myCAF",
                                  "7" = "myCAF")


p1 <- DimPlot(pbmc_Fib, reduction = "tsne", label = T,group.by = "fib.celltype.main")
p2 <- DimPlot(pbmc_Fib, reduction = "umap", label = T,group.by = "fib.celltype.main")
p = p1+p2
ggsave("6FibCelltype/5FibCluster_reduction.png", p, width = 7, height = 6)


p1 <- DimPlot(pbmc_Fib, reduction = "umap", label = T,group.by = "fib.celltype.main")
ggsave("6FibCelltype/5FibCluster_reduction_umap.png", p1, width = 7, height = 6)



#Marker基因标记
markers <- c("COL10A1","COL4A1","MMP3","IL13","TGFB1", #myCAF "IL4",
             "IL6","IL8","CXCL1","CXCL10","CCL2","CXCL14","CFD",#iCAF "CXCL8","GM-CSF",
             "CXCL12")#apCAF "CXCL12","CXCL14","CFD"

#将表达矩阵SCT转为RNA
DefaultAssay(pbmc_Fib) <- "RNA"
#p <- FeaturePlot(pbmc_Fib, reduction = "tsne", features = markers, ncol = 4)
#ggsave("6FibCelltype/1Fib_Markers_featureplot.pdf", p, width = 14, height = 8)
#dotplot
p <- DotPlot(pbmc_Fib, features = markers,group.by = "fib.celltype.main") + RotatedAxis()
ggsave("6FibCelltype/12Fib_Markers_dotplot.pdf", p, width = 14, height = 6)
#vlnplot
p <- VlnPlot(pbmc_Fib, features = markers, stack = T, flip = T,group.by = "fib.celltype.main") + NoLegend()
ggsave("6FibCelltype/13Fib_Markers_vlnplot.pdf", p, width = 14, height = 6)
DefaultAssay(pbmc_Fib) <- "SCT"


p1 <- DimPlot(pbmc_Fib, reduction = "tsne", label = T,group.by = "fib.celltype.main")
p2 <- DimPlot(pbmc_Fib, reduction = "umap", label = T,group.by = "fib.celltype.main")
p = p1+p2
ggsave("6FibCelltype/14FibCluster_reduction.png", p, width = 10, height = 6)



#保存
saveRDS(pbmc_Fib, file = "1ResData/8pbmc_Fib.rds")

###寻找每个群高变基因，

Idents(pbmc_Fib)<-pbmc_Fib$fib.celltype.main
DefaultAssay(pbmc_Fib) <- "RNA"
markers<-FindAllMarkers(pbmc_Fib,only.pos = T,min.pct = 0.1,logfc.threshold = 0.3,min.diff.pct = 0.1)
DefaultAssay(pbmc_Fib) <- "SCT"
write.csv(markers,file="6FibCelltype/markers0.1.0.3.0.1.csv")


####4.6GSVA未做####

setwd("F:/1/5LUSC/4scRNA")

dir.create("8GSVA")
rm(list = ls())
### 平均表达矩阵
scRNA <- readRDS("1ResData/7pbmc_Fib.rds")
DefaultAssay(scRNA) <- "RNA"
scRNA <- NormalizeData(scRNA)

#scRNA$group <- paste0(scRNA$subtype, ":", scRNA$celltype_major)
#table(scRNA$group)
table(scRNA@meta.data[["seurat_clusters"]])

expr <- AverageExpression(scRNA, assays = "RNA", slot = "data", group.by = "seurat_clusters")[[1]] #选取平均表达值
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)
print(n = 100,msigdbr_collections())

#读入感兴趣的条目名称
CAF1 = read.table(file ="8GSVA/1CAF2.txt", header=T, sep="\t", check.names=F)
CAF2 = read.table(file ="8GSVA/2CAF2.txt", header=T, sep="\t", check.names=F)
CAF3 = read.table(file ="8GSVA/3CAF2.txt", header=T, sep="\t", check.names=F)
CAFALL = rbind(CAF1,CAF2,CAF3)

#获取GO所有条目
genesetsBP <- msigdbr(species = "Homo sapiens", category = "C5") 
genesetsBP <- subset(genesetsBP, gs_subcat=="GO:BP", select = c("gs_name", "gene_symbol")) %>% as.data.frame()
genesetsBP <- split(genesetsBP$gene_symbol, genesetsBP$gs_name)

genesetsCC <- msigdbr(species = "Homo sapiens", category = "C5") 
genesetsCC <- subset(genesetsCC, gs_subcat=="GO:CC", select = c("gs_name", "gene_symbol")) %>% as.data.frame()
genesetsCC <- split(genesetsCC$gene_symbol, genesetsCC$gs_name)

genesetsMF <- msigdbr(species = "Homo sapiens", category = "C5") 
genesetsMF <- subset(genesetsMF, gs_subcat=="GO:MF", select = c("gs_name", "gene_symbol")) %>% as.data.frame()
genesetsMF <- split(genesetsMF$gene_symbol, genesetsMF$gs_name)

genesetsALL = append(genesetsBP,genesetsCC)
genesetsALL = append(genesetsALL,genesetsMF)

#提取感兴趣的条目
genesetsCAF = genesetsALL[CAFALL[,1]]

# gsva默认开启全部线程计算
gsva.res <- gsva(expr, genesetsCAF, method="ssgsea", kcdf="Gaussian", parallel.sz=1) 
saveRDS(gsva.res, "8GSVA/ssGSEA_genesetsCAF.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "8GSVA/ssGSEA_genesetsCAF.csv", row.names = F)
# 热图展示
pheatmap::pheatmap(gsva.res, scale = "row", filename = "8GSVA/ssGSEA_genesetsCAF.pdf", width = 15, height = 10)






#提取感兴趣的条目
genesetsCAF1 = genesetsALL[CAF1[,1]]
genesetsCAF2 = genesetsALL[CAF2[,1]]
genesetsCAF3 = genesetsALL[CAF3[,1]]

# gsva默认开启全部线程计算
gsva.res <- gsva(expr, genesetsCAF1, method="ssgsea", kcdf="Gaussian", parallel.sz=1) 
saveRDS(gsva.res, "8GSVA/ssGSEA_genesetsCAF1.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "8GSVA/ssGSEA_genesetsCAF1.csv", row.names = F)
# 热图展示
pheatmap::pheatmap(gsva.res, scale = "row", filename = "8GSVA/ssGSEA_genesetsCAF1.pdf", width = 20, height = 20)

# gsva默认开启全部线程计算
gsva.res <- gsva(expr, genesetsCAF2, method="ssgsea", kcdf="Gaussian", parallel.sz=1) 
saveRDS(gsva.res, "8GSVA/ssGSEA_genesetsCAF2.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "8GSVA/ssGSEA_genesetsCAF2.csv", row.names = F)
# 热图展示
pheatmap::pheatmap(gsva.res, scale = "row", filename = "8GSVA/ssGSEA_genesetsCAF2.pdf", width = 20, height = 20)

# gsva默认开启全部线程计算
gsva.res <- gsva(expr, genesetsCAF3, method="ssgsea", kcdf="Gaussian", parallel.sz=1) 
saveRDS(gsva.res, "8GSVA/ssGSEA_genesetsCAF3.rds")
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "8GSVA/ssGSEA_genesetsCAF3.csv", row.names = F)
# 热图展示
pheatmap::pheatmap(gsva.res, scale = "row", filename = "8GSVA/ssGSEA_genesetsCAF3.pdf", width = 20, height = 20)



table(scRNA@meta.data$seurat_clusters)

dim(scRNA)


####4.7拟时序分析#####
setwd("F:/1/5LUSC/4scRNA")
dir.create("9Monocle")

#创建monocle分析对象
pbmc_Fib <- readRDS("1ResData/8pbmc_Fib.rds")
pbmc_Fib <- NormalizeData(pbmc_Fib)
DimPlot(pbmc_Fib, group.by = "fib.celltype.main", label = T)
DimPlot(pbmc_Fib, group.by = "seurat_clusters", label = T)


#pbmc_Fib$celltype <- pbmc_Fib$copykat.pred
colnames(pbmc_Fib@meta.data)

#monocle不推荐使用slot="data"
data <- GetAssayData(pbmc_Fib, assay = "RNA", slot = "counts")
#可以选取特定的meta.data
#pd <- new('AnnotatedDataFrame', data = pbmc_Fib@meta.data[,c(1,35,36)])
pd <- new('AnnotatedDataFrame', data = pbmc_Fib@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit = 0.5,
                        expressionFamily = negbinomial.size())

#数据预处理
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=8)

#选择排序基因
disp_table <- dispersionTable(mycds)
order.genes <- subset(disp_table, mean_expression >= 0.005 & dispersion_empirical >= 
                        1 * dispersion_fit) %>% pull(gene_id) %>% as.character()
#增加或剔除特定基因，默认不运行
if(F){
  order.genes.adj <- order.genes
  #增加基因
  add.genes <- c('CCR7','LEF1','TCF7','SELL')
  order.genes.adj <- unique(c(order.genes.adj, add.genes))
  #减少基因
  del.genes <- unique(c(grep("^MT-", rownames(pbmc_Fib), v=T), "NEAT1","TMSB4X","TMSB10"))
  order.genes.adj <- setdiff(order.genes.adj, del.genes)
  #改变高变基因
  order.genes <- order.genes.adj
}
mycds <- setOrderingFilter(mycds, order.genes)
#设置排序基因
p <- plot_ordering_genes(mycds)
ggsave("9Monocle/OrderGenes.pdf", p, width = 8, height = 6)

#降维排序
mycds <- reduceDimension(mycds, max_components = 2, reduction_method = 'DDRTree')
mycds <- orderCells(mycds)

display.brewer.all()
cols_1<-brewer.pal(8, 'Dark2')

Biocols = c('#AB3282', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
            '#E95C59', '#E59CC4','#E5D2DD' , '#23452F', '#BD956A', '#8C549C', '#585658',
            '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
            '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
            '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
            '#968175',"DarkOrchid")


p1 <- plot_cell_trajectory(mycds, color_by = "State")+  
  scale_colour_manual(
    values = Biocols
    # aesthetics = c("colour", "fill")
  )


#结果可视化
#State轨迹分布图
p1 <- plot_cell_trajectory(mycds, color_by = "State")+  
  scale_colour_manual(
    values = Biocols
    # aesthetics = c("colour", "fill")
  )
ggsave("9Monocle/Trajectory_State.pdf", plot = p1, width = 10, height = 6.5)
#Pseudotime轨迹图
p2 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
ggsave("9Monocle/Trajectory_Pseudotime.pdf", plot = p2, width = 10, height = 6.5)
#Celltype轨迹分布图
p3 <- plot_cell_trajectory(mycds, color_by = "fib.celltype.main")+  
  scale_colour_manual(
    values = Biocols
    # aesthetics = c("colour", "fill")
  )
ggsave("9Monocle/Trajectory_Celltype.pdf", plot = p3, width = 10, height = 6.5)
#Sample轨迹分布图
p4 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")+  
  scale_colour_manual(
    values = Biocols
    # aesthetics = c("colour", "fill")
  )
ggsave("9Monocle/Trajectory_seurat_clusters.pdf", plot = p4, width = 10, height = 6.5)

#调整排序重新画图（看情况而定，不是每次都要调整排序）
if(F){
  mycds <- orderCells(mycds, root_state = 5)
  p1 <- plot_cell_trajectory(mycds, color_by = "State")
  ggsave("Trajectory_State5.pdf", plot = p1, width = 10, height = 6.5)
  #Pseudotime轨迹图
  p2 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
  ggsave("Trajectory_Pseudotime5.pdf", plot = p2, width = 10, height = 6.5)
  #Celltype轨迹分布图
  p3 <- plot_cell_trajectory(mycds, color_by = "celltype")
  ggsave("Trajectory_Celltype5.pdf", plot = p3, width = 10, height = 6.5)
  #Sample轨迹分布图
  p4 <- plot_cell_trajectory(mycds, color_by = "orig.ident")
  ggsave("Trajectory_Sample5.pdf", plot = p4, width = 10, height = 6.5)
}

#提取拟时信息给seurat对象
pdata <- Biobase::pData(mycds)
pbmc_Fib <- AddMetaData(pbmc_Fib, metadata = pdata[,c("Pseudotime","State")])

#以下未做
#寻找拟时差异基因一
#拟时差异基因分析
Time_diff <- differentialGeneTest(mycds, cores = 10, fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.csv(Time_diff, "9Monocle/Time_diff_all.csv", row.names = F)
#显著差异基因的作图
Time_genes <- Time_diff[order(Time_diff$qval), "gene_short_name"][1:100] #提取qval最小的100个基因
p = plot_pseudotime_heatmap(mycds[Time_genes,], num_clusters=3, show_rownames=T, return_heatmap=T)
ggsave("9Monocle/Time_heatmap.pdf", p, width = 5, height = 10)
#显著差异基因按热图结果排序并保存
hp.genes <- p$tree_row$labels[p$tree_row$order]
Time_diff_sig <- Time_diff[hp.genes, c("gene_short_name", "pval", "qval")]
write.csv(Time_diff_sig, "9Monocle/Time_diff_sig.csv", row.names = F)

#寻找拟时差异基因二
Idents(pbmc_Fib) <- "State"
deg <- FindAllMarkers(pbmc_Fib, logfc.threshold = 0.5, only.pos = T)
State_genes <- group_by(deg, cluster) %>% top_n(25, avg_log2FC) %>% pull(gene) %>% unique()
p = plot_pseudotime_heatmap(mycds[State_genes,], num_clusters=3, show_rownames=T, return_heatmap=T)
ggsave("9Monocle/State_heatmap.pdf", p, width = 5, height = 10)
#显著差异基因按热图结果排序并保存
hp.genes <- p$tree_row$labels[p$tree_row$order]
State_diff_sig <- deg[hp.genes, ]
write.csv(State_diff_sig, "9Monocle/State_diff_sig.csv", row.names = F)

#指定基因的可视化
s.genes <- c("IFNG","FOXP3","CCR7")
p1 <- plot_genes_jitter(mycds[s.genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(mycds[s.genes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "State")
plotc <- p1|p2|p3
ggsave("Genes_Jitterplot.pdf", plot = plotc, width = 16, height = 8)
#模型基因的可视化
s.genes <- c("S100A9","CFH","PPP1R16A","POP4","PRDX6","GAGE2A",
             "UGT2B15","LDHA","TMEM106C","ALAS1","RTN3","DNAJB4")
plot_genes_in_pseudotime(mycds[s.genes,], color_by = "Pseudotime")
ggsave("Model_gene5.pdf", plot = plotc, width = 8, height = 16)
#保存
saveRDS(pbmc_Fib, file = "1ResData/9pbmc.pseudotime.rds")

####4.8转录因子分析####

setwd("F:/1/5LUSC/4scRNA")
#创建PySCENIC分析对象
dir.create("10PySCENIC")
pbmc_Fib <- readRDS("1ResData/8pbmc_Fib.rds")
write.csv(t(as.matrix(pbmc_Fib@assays$RNA@counts)),file = "10PySCENIC/pbmc_Fib_exp.csv")

#以下内容在linux系统中运行
#转化为loom文件，Linux下的python脚本
#编辑脚本
vim trans.py
#输入以下内容
import os, sys
os.getcwd()
os.listdir(os.getcwd())
import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("pbmc_Fib_exp.csv");#R中导出的表达矩阵
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("pbmc_Fib.loom",x.X.transpose(),row_attrs,col_attrs)

#保存并退出
#运行trans.py
python trans.py
ls
#这样在文件夹中会出现pbmc_Fib.loom文件，就是接下来输入pyscenic的文件。

pyscenic grn --num_workers 10 \
--sparse \
--method grnboost2 \
--output pbmc_Fib.adj.csv \
pbmc_Fib.loom \
hs_hgnc_tfs.txt
#这一步的目的
#推断转录因子与提供的表达矩阵基因的共表达模块，基于grnboost2，R中时GENIE3

pyscenic ctx --num_workers 10 \
--output pbmc_Fib.regulons.csv \
--expression_mtx_fname pbmc_Fib.loom \
--all_modules \
--mask_dropouts \
--mode "dask_multiprocessing" \
--min_genes 10 \
--annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
pbmc_Fib.adj.csv \
hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
#这一步的目的
#进行TF-motif富集分析，识别直接靶标
#得到转录因子(TF)与其对应的直接作用的靶点,称为regulon(每一个regulon是1个TF和其调控的靶基因)


pyscenic aucell --num_workers 3 \
--output pbmc_Fib_SCENIC.loom \
pbmc_Fib.loom \
pbmc_Fib.regulons.csv
#这一步的目的
#使用AUCell对每个细胞的每个regulon活性进行评分。


#以下在R中运行
setwd("F:/1/5LUSC/4scRNA/10PySCENIC")
#加载分析包
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
#可视化相关包，多加载点没毛病
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)
library(pheatmap)

sce_SCENIC <- open_loom("pbmc_Fib_SCENIC.loom")
#exprMat <- get_dgem(sce_SCENIC)#从sce_SCENIC文件提取表达矩阵
# exprMat_log <- log2(exprMat+1) # log处理
regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
class(regulons)

regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC)


#修改成自己的
human_data <- readRDS("F:/1/5LUSC/4scRNA/1ResData/8pbmc_Fib.rds")
colnames(human_data@meta.data)
cellinfo <- human_data@meta.data[,c('fib.celltype.main','orig.ident',"nFeature_RNA","nCount_RNA")]#细胞meta信息
colnames(cellinfo)=c('celltype', 'orig.ident','nGene' ,'nUMI')
######计算细胞特异性TF
cellTypes <-  as.data.frame(subset(cellinfo,select = 'celltype'))
selectedResolution <- "celltype"
sub_regulonAUC <- regulonAUC

rss <- calcRSS(AUC=getAUC(sub_regulonAUC),
               cellAnnotation=cellTypes[colnames(sub_regulonAUC),
                                        selectedResolution])

rss=na.omit(rss)
rssPlot <- 
  plotRSS(
    rss,
    zThreshold = 1,
    cluster_columns = FALSE,
    order_rows = TRUE,
    thr=0.01,
    varName = "cellType",
    col.low = '#330066',
    col.mid = '#66CC66',
    col.high = '#FFCC33')
rssPlot


rss_data <- rssPlot$plot$data
#devtools::install_github("XiaoLuo-boy/ggheatmap")
library(ggheatmap)
library(reshape2)
rss_data<-dcast(rss_data, 
                Topic~rss_data$cellType,
                value.var = 'Z')
rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]
colnames(rss_data)
col_ann <- data.frame(group= c(rep("apCAF",1),
                               rep("myCAF",1),
                               rep("iCAF",1)))#列注释
rownames(col_ann) <- colnames(rss_data)
groupcol <- c("#D9534F", "#96CEB4", "#CBE86B", "#EDE574", "#0099CC")
names(groupcol) <- c("apCAF","myCAF","iCAF")
col <- list(group=groupcol)

text_columns <- sample(colnames(rss_data),0)#不显示列名

p <- ggheatmap(rss_data,color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
               cluster_rows = T,cluster_cols = F,scale = "row",
               annotation_cols = col_ann,
               annotation_color = col,
               legendName="Relative value",
               text_show_cols = text_columns)
p


next_regulonAUC <- regulonAUC[,match(colnames(human_data),colnames(regulonAUC))]
dim(next_regulonAUC)

regulon_AUC <- regulonAUC@NAMES
human_data@meta.data = cbind(human_data@meta.data ,t(assay(next_regulonAUC[regulon_AUC,])))

#自己选定感兴趣的或者比较重要的转录因子，这里我是随机的
TF_plot <- c("BRCA1(+)")

DotPlot(human_data, features = TF_plot)+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(hjust =1,vjust=1, angle = 45))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))


DotPlot(human_data, features = TF_plot, group.by = 'fib.celltype.main')+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(hjust =1,vjust=1, angle = 45))+
  theme(legend.direction = "horizontal", 
        legend.position = "bottom")+
  labs(x=NULL,y=NULL)

colnames(human_data@meta.data)



cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution])
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 


regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)


hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byGroup_Scaled, name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6),
                                   show_row_names = F)) 
hm  #可视化所有的TF


#提取转录因子靶基因
source("F:/1/5LUSC/4scRNA/sc_function.R")
TFtargets <- Regulons2Genesets("regulons.csv")
saveRDS(TFtargets, file = "TFtargets.rds")

regulonAUC = assay(regulonAUC)

#AUC矩阵二进制转化
#读取regulons活性矩阵
#regulonAUC <- read.csv("auc_mtx.csv", row.names = 1,check.names = F)

#rownames(regulonAUC) = regulonAUC[,1]

#regulonAUC <-open_loom("pbmc_Fib_SCENIC.loom")
regulonAUC <- t(regulonAUC)
saveRDS(regulonAUC, file = "regulonAUC.rds")
#计算活性阈值，运行时间比较长
bin.T <- AUCell_exploreThresholds(regulonAUC,
                                  smallestPopPercent=0.25,
                                  assignCells=TRUE, 
                                  plotHist=FALSE,
                                  verbose=FALSE)
saveRDS(bin.T, "binThresholds.rds")
#二进制转化
regulonBin <- lapply(rownames(regulonAUC), function(reg){
  as.numeric(colnames(regulonAUC) %in% bin.T[[reg]][["assignment"]])
})
regulonBin <- do.call("rbind", regulonBin)
dimnames(regulonBin) <- list(rownames(regulonAUC), colnames(regulonAUC))
saveRDS(regulonBin, "regulonBin.rds")

#Regulon特异性分析
#RSS
sco <- readRDS("F:/1/5LUSC/4scRNA/1ResData/8pbmc_Fib.rds")


group = substr(colnames(regulonAUC), nchar(colnames(regulonAUC))-1, nchar(colnames(regulonAUC))-1)
dim(regulonAUC)
data1 = regulonAUC[,group == "+"]
dim(data1)
data0 = regulonAUC[,group == "-"]
dim(data0)
data = cbind(data1,data0)
dim(data)
colnames(data1) = substring(colnames(data1),1, nchar(colnames(data1))-3)
regulonAUC = t(data1)

rss <- calcRSS(regulonAUC, sco$fib.celltype.main)



#rownames(rss) = sub('...$','',rownames(rss))
rss.plot <- rssplot(rss = rss)
p <- wrap_plots(rss.plot, ncol = 3)
ggsave("RSS.pdf", p, width = 9, height = 8, limitsize = F)

#差异分析
sco[["Regulon"]] <- CreateAssayObject(counts = regulonAUC)
sco[["binRegulon"]] <- CreateAssayObject(counts = regulonBin)
DefaultAssay(sco) <- "Regulon"
sco <- ScaleData(sco, features = rownames(sco))
Idents(sco) <- "fib.celltype.main"
deg <- FindAllMarkers(sco, only.pos = T, logfc.threshold = 0)

#rownames(deg)= sub('...$','',rownames(deg))

top <- group_by(deg, cluster) %>% top_n(10, avg_log2FC) %>% pull(gene) %>% unique()


p <- DoHeatmap(sco, features = top, label = F)
ggsave("DEG.png", p, width = 12, height = 6.5, limitsize = F)


#Regulon概览
rn  = "TFE3(-)"
tf = "TFE3"
gp = "fib.celltype.main"
p1 <- DimPlot(sco, reduction = "umap", group.by = "fib.celltype.main", label = T) + NoLegend()
DefaultAssay(sco) <- "binRegulon"
p2 <- FeaturePlot(sco, reduction = "umap", features = rn) + ggtitle(paste0(rn,"_binRAS"))
DefaultAssay(sco) <- "Regulon"
p3 <- FeaturePlot(sco, reduction = "umap", features = rn) + ggtitle(paste0(rn,"_RAS"))
DefaultAssay(sco) <- "RNA"
p4 <- FeaturePlot(sco, reduction = "umap", features = tf) + ggtitle(paste0(tf,"_Expression"))
df1 <- data.frame(cluster=sco@meta.data[,gp,drop=T], auc=sco@assays$Regulon@counts[rn,], row.names = colnames(sco))
p5 <- ggboxplot(df1, x='cluster', y='auc', fill='cluster', bxp.errorbar = T, outlier.shape = NA) +
  ggtitle(paste0(rn,"_RAS")) + NoLegend() + theme(plot.title = element_text(hjust = 0.5)) + RotatedAxis() 
df2 <- data.frame(cluster=sco@meta.data[,gp,drop=T], expr=sco@assays$RNA@data[tf,], row.names = colnames(sco))
p6 <- ggboxplot(df2, x='cluster', y='expr', fill='cluster', bxp.errorbar = T, outlier.shape = NA) +
  ggtitle(paste0(tf,"_Expression")) + NoLegend() + theme(plot.title = element_text(hjust = 0.5)) + RotatedAxis() 
p <- (p1|p3|p4)/(p2|p5|p6)
ggsave(paste0(tf,"_overview.pdf"), p, width = 16, height = 9)  


####4.9细胞通讯分析####

setwd("F:/1/5LUSC/4scRNA")

options(stringsAsFactors = FALSE)
dir.create("11CellChat")
#创建CellChat对象
sco.brca.9s <- readRDS("1ResData/5pbmc_cluster.rds")
#s.cells <- c('CAFs MSC iCAF-like','CAFs myCAF-like','Macrophage','Monocyte','NK cells','NKT cells','T cells CD4+','T cells CD8+')
#sco <- subset(sco.brca.9s, subtype == "ER+"&celltype_minor %in% s.cells)
#sco <- subset(sco.brca.9s, subtype == "TNBC"&celltype_minor %in% s.cells)

table(sco.brca.9s@meta.data[["celltype.main"]])
#table(sco.brca.9s@meta.data[["copykat.pred"]])
#导出meta.data信息
write.table(sco.brca.9s@meta.data,"11CellChat/sco.brca.9s.metadata_change.xlsx",sep="\t",quote=F,col.names=T)

pbmc_Fib <- readRDS("1ResData/8pbmc_Fib.rds")
write.table(pbmc_Fib@meta.data,"11CellChat/pbmc_Fib.metadata_change.xlsx",sep="\t",quote=F,col.names=T)

#自行修改
#增加meta.data的信息
metadata_celltype_finnal <-fread("11CellChat/metadata_change2.xlsx", sep="\t",header = T,check.names = F)
#将第一列转变为行名
metadata_celltype_finnal <- metadata_celltype_finnal %>% column_to_rownames("V1")
dim(metadata_celltype_finnal)
#添加
sco.brca.9s <- AddMetaData(sco.brca.9s, metadata = metadata_celltype_finnal)
table(sco.brca.9s@meta.data$celltype.final)
#转换
sco <- sco.brca.9s
sco <- Seurat::NormalizeData(sco)
sco$celltype <- sco$celltype.final
cellchat <- createCellChat(sco@assays$RNA@data, meta = sco@meta.data, group.by = "celltype")
table(sco$celltype,sco$orig.ident)
#分析初始化
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human    #小鼠用CellChatDB.mouse
showDatabaseCategory(CellChatDB)
#可以选择数据库子集用于分析
#CellChatDB <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB

#数据预处理
#提取细胞通讯信号基因
cellchat <- subsetData(cellchat)
#识别过表达的配体、受体基因
future::plan("multisession", workers = 40)
cellchat <- identifyOverExpressedGenes(cellchat)
#识别过表达的配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
future::plan("sequential", workers = 6)
#数据校正（可选）
#cellchat <- projectData(cellchat, PPI.human)

#计算互作概率
#配体受体层面计算互作概率
cellchat <- computeCommunProb(cellchat, raw.use = T) #如果用非校正数据，raw.use = T
cellchat <- filterCommunication(cellchat, min.cells = 10)
#信号通路层面计算互作概率
cellchat <- computeCommunProbPathway(cellchat)
#提取细胞互作结果
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "11CellChat/CommunProb.csv", row.names = F)
df.netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netP, "11CellChat/CommunProbPathway.csv", row.names = F)
#聚合通讯概率或数量（细胞间会多个配体-受体互作关系，这一步是统计数量或聚合概率）
cellchat <- aggregateNet(cellchat)

#细胞互作概览
pdf("11CellChat/fig1a_netVisual_overview_all.pdf", width = 8, height = 6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")
dev.off()

pdf("11CellChat/fig1b_netVisual_overview_split.pdf", width = 6, height = 5)
mat <- cellchat@net$weight
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

#按信号通路探索细胞互作
mypathways <- cellchat@netP$pathways
#mypathways <- mypathways[1:10] #用10个信号通路演示，实际分析时不运行这行代码！

#细胞间互作概率展示
#网络图展示
pdf("11CellChat/fig3a_netVisual_pathways_circle.pdf", width = 6, height = 5)
for(pathways.show in mypathways){
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
}
dev.off()
#和弦图展示互作概览
pdf("11CellChat/fig3b_netVisual_pathways_chord.pdf", width = 10, height = 8)
for(pathways.show in mypathways){
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
}
dev.off()
#热图展示互作概览
pdf("11CellChat/fig3c_netVisual_pathways_heatmap.pdf", width = 10, height = 6.5)
for(pathways.show in mypathways){
  p <- netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  plot(p)
}
dev.off()

#信号通路内的配体-受体
if(!file.exists("11CellChat/pathways_detail")) dir.create("11CellChat/pathways_detail")
for(pathways.show in mypathways){
  pdf(paste0("11CellChat/pathways_detail/", pathways.show, ".pdf"), width = 8, height = 6.5)
  #配体-受体贡献度展示
  netAnalysis_contribution(cellchat, signaling = pathways.show)
  pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)$interaction_name
  for(LR.show in pairLR){
    #网络图展示细胞间的配体-受体互作
    netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
  }
  for(LR.show in pairLR){
    #和弦图展示细胞间的配体-受体互作
    netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
  }
  dev.off()
}

#按配体-受体探索细胞互作
levels(cellchat@idents)

#显示所有的细胞间配体-受体互作
p <- netVisual_bubble(cellchat, sources.use = 1:2, targets.use = 1:2, remove.isolate = FALSE)
ggsave("11CellChat/fig4a_CCI_all.pdf", p, width = 10, height = 15, limitsize = F)

p <- netVisual_bubble(cellchat, sources.use = c(1,5,8), targets.use = c(2,3,4,6,7), remove.isolate = FALSE,
                      font.size = 17,font.size.title = 17,line.size = 0.3,angle.x = 45)
ggsave("11CellChat/fig4a_CCI_1.pdf", p, width = 10, height = 20, limitsize = F)

p <- netVisual_bubble(cellchat, sources.use = c(2,3,4,6,7), targets.use = c(1,5,8), remove.isolate = FALSE,
                      font.size = 17,font.size.title = 17,line.size = 0.3,angle.x = 45)
ggsave("11CellChat/fig4a_CCI_2.pdf", p, width = 10, height = 20, limitsize = F)



#指定细胞间的配体-受体互作

pairLR.use <- c("FGF1_FGFR1","FGF1_FGFR2","FGF1_FGFR3","FGF1_FGFR4","FGF10_FGFR1","FGF10_FGFR2","FGF2_FGFR1","FGF2_FGFR2","FGF2_FGFR3","FGF2_FGFR4","FGF4_FGFR1","FGF4_FGFR2",
                "FGF4_FGFR3","FGF4_FGFR4","FGF5_FGFR1","FGF5_FGFR2","FGF5_FGFR3","FGF5_FGFR4","FGF7_FGFR1","FGF7_FGFR2","FGF9_FGFR1","FGF9_FGFR2","FGF9_FGFR3","FGF9_FGFR4",
                "FGF21_FGFR1","FGF21_FGFR3","AREG_(EGFR+ERBB2)","AREG_EGFR","EGF_(EGFR+ERBB2)","EGF_EGFR","HBEGF_(EGFR+ERBB2)","HBEGF_EGFR","PGF_VEGFR1","TGFA_(EGFR+ERBB2)",
                "TGFA_EGFR","VEGFA_VEGFR1","VEGFA_VEGFR1R2","VEGFA_VEGFR2","VEGFB_VEGFR1","VEGFC_VEGFR2","VEGFC_VEGFR2R3","VEGFC_VEGFR3","HGF_MET")
pairLR.use <- data.frame(interaction_name = pairLR.use)
p <- netVisual_bubble(cellchat, sources.use = 1:6, targets.use = 7, pairLR.use =  pairLR.use, remove.isolate = FALSE)
ggsave("11CellChat/Tumor/targets.pdf", p, width = 6, height = 10, limitsize = F)
p <- netVisual_bubble(cellchat, sources.use = 7, targets.use = 1:6, pairLR.use =  pairLR.use, remove.isolate = FALSE)
ggsave("11CellChat/Tumor/sources.pdf", p, width = 6, height = 10, limitsize = F)


#指定细胞和信号通路的配体-受体互作
p <- netVisual_bubble(cellchat, sources.use = 1:4, targets.use = 7:8, signaling = c("CCL","CXCL"), remove.isolate = FALSE)
ggsave("11CellChat/fig4c_CCI_subcell_subpathway.pdf", p, width = 6, height = 6, limitsize = F)

#指定细胞和配体-受体
#pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CXCL"))
pairLR.use <- c('CXCL2_CXCR2','CXCL12_CXCR4','CCL3_CCR5','CCL4_CCR5' )
pairLR.use <- data.frame(interaction_name = pairLR.use)
p <- netVisual_bubble(cellchat, sources.use = 1:4, targets.use = 7, pairLR.use =  pairLR.use, remove.isolate = FALSE)
ggsave("11CellChat/fig4d_CCI_subcell_subLR.pdf", p, width = 6, height = 4, limitsize = F)

#配体、受体基因表达
p <- plotGeneExpression(cellchat, signaling = "CXCL")
ggsave("11CellChat/fig5a_GeneExpression_violin_sig.pdf", p, width = 10, height = 9, limitsize = F)
p <- plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = FALSE)
ggsave("11CellChat/fig5b_GeneExpression_violin_all.pdf", p, width = 10, height = 9, limitsize = F)

#信号网络生态分析
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pdf("11CellChat/fig6_signalingRole.pdf", width = 6, height = 4.5)
for(pathways.show in mypathways){
  netAnalysis_signalingRole_network(cellchat, signaling=pathways.show, width=8, height=2.5, font.size=10)
}
dev.off()

#保存结果
saveRDS(cellchat, file = "1ResData/8pbmc_CellChat.rds")

#读取
cellchat <- readRDS("1ResData/8pbmc_CellChat.rds")



####5.1单因素cox,未运行####


#引用包
library(limma)
library(survival)	
library(survminer)

expFile="TCGA_LUSC.txt"      #表达数据文件
cliFile="time_LUSC.txt"            #生存数据文件
setwd("F:/1/5LUSC/5Survival")    #设置工作目录

#读取表达文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data2=data
data=t(data)

#读取生存数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli$futime=cli$futime/365


markers <- read_csv("markers0.1.0.3.0.1.csv")
apCAF = subset(markers,cluster == "apCAF")

samegeneapCAF = intersect(colnames(data),as.list(apCAF[,1])[["...1"]])
data = data[,samegeneapCAF]

#rownames(data) = gsub('GSE116918_', '', rownames(data))
#rownames(data) = gsub('TCGA_', '', rownames(data))

#数据合并
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",row.names(data))
sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
rt=cbind(cli, data)

#对基因进行循环，查找预后相关的基因
outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  #cox分析
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<0.05){
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}

#输出单因素的结果
write.table(outTab, file="uniCoxapCAF.txt", sep="\t", row.names=F, quote=F)

#输出单因素显著基因的表达量
sigGenes=as.vector(outTab[,1])
uniSigExp=data2[sigGenes,]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
#write.table(uniSigExp,file="uniSigExpapCAF.txt",sep="\t",row.names=F,quote=F)

#输出单因素显著基因表达数据和生存数据合并的文件
uniSigExpTime=rt[,c("futime","fustat",sigGenes)]
uniSigExpTime=cbind(id=row.names(uniSigExpTime),uniSigExpTime)
#write.table(uniSigExpTime, file="uniSigExpTimeapCAF.txt", sep="\t", row.names=F, quote=F)

############定义森林图函数
bioForest=function(coxFile=null,forestFile=null,forestCol=null){
  #读取输入文件
  rt <- read.table(coxFile,header=T,sep="\t",check.names=F,row.names=1)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #输出图形
  height=nrow(rt)/13+6
  pdf(file=forestFile, width = 7,height = height)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #绘制森林图左边的基因信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  
  #绘制森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol[1], forestCol[2])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)
  axis(1)
  dev.off()
}

#调用函数, 对单因素的结果进行可视化, 绘制森林图
bioForest(coxFile="uniCoxapCAF.txt", forestFile="forestapCAF.pdf", forestCol=c("Firebrick3","MediumSeaGreen"))








#读取表达文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data2=data
data=t(data)

#读取生存数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli$futime=cli$futime/365


myCAF = subset(markers,cluster == "myCAF")
samegenemyCAF = intersect(colnames(data),as.list(myCAF[,1])[["...1"]])
data = data[,samegenemyCAF]

#rownames(data) = gsub('GSE116918_', '', rownames(data))
#rownames(data) = gsub('TCGA_', '', rownames(data))

#数据合并
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",row.names(data))
sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
rt=cbind(cli, data)

#对基因进行循环，查找预后相关的基因
outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  #cox分析
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<0.05){
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}

#输出单因素的结果
write.table(outTab, file="uniCoxmyCAF.txt", sep="\t", row.names=F, quote=F)

#调用函数, 对单因素的结果进行可视化, 绘制森林图
bioForest(coxFile="uniCoxmyCAF.txt", forestFile="forestmyCAF.pdf", forestCol=c("Firebrick3","MediumSeaGreen"))






#读取表达文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data2=data
data=t(data)

#读取生存数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli$futime=cli$futime/365

iCAF = subset(markers,cluster == "iCAF")
samegeneiCAF = intersect(colnames(data),as.list(iCAF[,1])[["...1"]])
data = data[,samegeneiCAF]

#rownames(data) = gsub('GSE116918_', '', rownames(data))
#rownames(data) = gsub('TCGA_', '', rownames(data))

#数据合并
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",row.names(data))
sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
rt=cbind(cli, data)

#对基因进行循环，查找预后相关的基因
outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  #cox分析
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<0.05){
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}

#输出单因素的结果
write.table(outTab, file="uniCoxiCAF.txt", sep="\t", row.names=F, quote=F)

#调用函数, 对单因素的结果进行可视化, 绘制森林图
bioForest(coxFile="uniCoxiCAF.txt", forestFile="forestiCAF.pdf", forestCol=c("Firebrick3","MediumSeaGreen"))







####5.2GSVA计算每个样本的评分####

#引用包
library(readr)
library(limma)
library(sva)
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)
setwd("F:/1/5LUSC/5Survival")
tcgaExpFile="TCGA_LUSC.txt"         #TCGA表达数据文件

#gene1 = read.table("uniCoxapCAF.txt", header=T, sep="\t", check.names=F)
#gene2 = read.table("uniCoxmyCAF.txt", header=T, sep="\t", check.names=F)
#gene3 = read.table("uniCoxiCAF.txt", header=T, sep="\t", check.names=F)
markers <- read_csv("markers0.1.0.3.0.1.csv")
gene1 = subset(markers,cluster == "apCAF")
gene2 = subset(markers,cluster == "myCAF")
gene3 = subset(markers,cluster == "iCAF")

#gene1 = as.list(gene1[,1])[["...1"]]
#gene2 = as.list(gene2[,1])[["...1"]]
#gene3 = as.list(gene3[,1])[["...1"]]
#读取TCGA基因表达文件,并对数据进行处理
rt = read.table(tcgaExpFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
tcga=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
tcga=avereps(tcga)

#samegeneapCAF = intersect(rownames(tcga),as.list(apCAF[,1])[["...1"]])
#expapCAF = tcga[samegeneapCAF,]

#samegenemyCAF = intersect(rownames(tcga),as.list(myCAF[,1])[["...1"]])
#expmyCAF = tcga[samegenemyCAF,]

#samegeneiCAF = intersect(rownames(tcga),as.list(iCAF[,1])[["...1"]])
#expiCAF = tcga[samegeneiCAF,]

geneSets = list()
#geneSets[1] = as.list(apCAF[,1])
#names(geneSets)[1] <- "apCAF"
#geneSets[2] = as.list(myCAF[,1])
#names(geneSets)[2] <- "myCAF"
#geneSets[3] = as.list(iCAF[,1])
#names(geneSets)[3] <- "iCAF"

#geneSets[1] = list(gene1[,1])
#names(geneSets)[1] <- "apCAF"

#geneSets[2] = list(gene2[,1])
#names(geneSets)[2] <- "myCAF"

#geneSets[3] = list(gene3[,1])
#names(geneSets)[3] <- "iCAF"

gene1 = intersect(rownames(tcga),as.vector(unlist(gene1[,1])))
geneSets[1] = list(gene1)
names(geneSets)[1] <- "apCAF"

gene2 = intersect(rownames(tcga),as.vector(unlist(gene2[,1])))
geneSets[2] = list(gene2)
names(geneSets)[2] <- "myCAF"

gene3 = intersect(rownames(tcga),as.vector(unlist(gene3[,1])))
geneSets[3] = list(gene3)
names(geneSets)[3] <- "iCAF"

#ssGSEA分析
ssgseaScore=gsva(tcga, geneSets, method='gsva', kcdf='Gaussian', abs.ranking=T)


#对ssGSEA打分进行矫正
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
#输出ssGSEA打分结果
ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="ssGSEA.result.txt",sep="\t",quote=F,col.names=F)

riskapCAF=as.vector(ifelse(ssgseaScore[1,]>median(ssgseaScore[1,]),"high","low"))
OutapCAF= cbind(id=colnames(ssgseaScore),ssgseaScore[1,],riskapCAF)
colnames(OutapCAF)[2]="riskScore"
colnames(OutapCAF)[3]="risk"
write.table(OutapCAF,file="OutapCAF.txt",sep="\t",quote=F,col.names=T,row.names = F)

riskmyCAF=as.vector(ifelse(ssgseaScore[2,]>median(ssgseaScore[2,]),"high","low"))
OutmyCAF= cbind(id=colnames(ssgseaScore),ssgseaScore[2,],riskmyCAF)
colnames(OutmyCAF)[2]="riskScore"
colnames(OutmyCAF)[3]="risk"
write.table(OutmyCAF,file="OutmyCAF.txt",sep="\t",quote=F,col.names=T,row.names = F)

riskiCAF=as.vector(ifelse(ssgseaScore[3,]>median(ssgseaScore[3,]),"high","low"))
OutiCAF= cbind(id=colnames(ssgseaScore),ssgseaScore[3,],riskiCAF)
colnames(OutiCAF)[2]="riskScore"
colnames(OutiCAF)[3]="risk"
write.table(OutiCAF,file="OutiCAF.txt",sep="\t",quote=F,col.names=T,row.names = F)


####5.3合并时间####

library(limma)                       #引用包
expFile_tcga="OutmyCAF.txt"      #表达数据文件
cliFile_tcga="time_LUSC.txt"                   #临床数据

setwd("F:/1/5LUSC/5Survival")   #设置工作目录

#读取表达文件，并对输入文件整理
rt=read.table(expFile_tcga,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,1:ncol(rt)]

#读取生存数据
cli=read.table(cliFile_tcga,sep="\t",check.names=F,header=T,row.names=1)     #读取临床文件
cli$futime[cli$futime<=0]=1
cli$futime=cli$futime/365
#rownames(exp) = gsub('GSE116918_', '', rownames(exp))
#rownames(exp) = gsub('TCGA_', '', rownames(exp))
rownames(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(exp))
sample = intersect(rownames(exp),rownames(cli))
exp = exp[sample,]
cli = cli[sample,]
out=cbind(cli,exp[,2:3])
#colnames(out)[3] = "risk"
write.table(data.frame(ID=rownames(out),out),file="expTime.myCAF.txt",sep="\t",row.names=F,quote=F)



expFile_tcga="OutiCAF.txt"      #表达数据文件

#读取表达文件，并对输入文件整理
rt=read.table(expFile_tcga,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,1:ncol(rt)]

#读取生存数据
cli=read.table(cliFile_tcga,sep="\t",check.names=F,header=T,row.names=1)     #读取临床文件
cli$futime[cli$futime<=0]=1
cli$futime=cli$futime/365
#rownames(exp) = gsub('GSE116918_', '', rownames(exp))
#rownames(exp) = gsub('TCGA_', '', rownames(exp))
rownames(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(exp))
sample = intersect(rownames(exp),rownames(cli))
exp = exp[sample,]
cli = cli[sample,]
out=cbind(cli,exp[,2:3])
#colnames(out)[3] = "risk"
write.table(data.frame(ID=rownames(out),out),file="expTime.iCAF.txt",sep="\t",row.names=F,quote=F)





expFile_tcga="OutapCAF.txt"      #表达数据文件

#读取表达文件，并对输入文件整理
rt=read.table(expFile_tcga,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,1:ncol(rt)]

#读取生存数据
cli=read.table(cliFile_tcga,sep="\t",check.names=F,header=T,row.names=1)     #读取临床文件
cli$futime[cli$futime<=0]=1
cli$futime=cli$futime/365
#rownames(exp) = gsub('GSE116918_', '', rownames(exp))
#rownames(exp) = gsub('TCGA_', '', rownames(exp))
rownames(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(exp))
sample = intersect(rownames(exp),rownames(cli))
exp = exp[sample,]
cli = cli[sample,]
out=cbind(cli,exp[,2:3])
#colnames(out)[3] = "risk"
write.table(data.frame(ID=rownames(out),out),file="expTime.apCAF.txt",sep="\t",row.names=F,quote=F)





####5.4生存分析####

#引用包
library(survival)
library(survminer)
setwd("F:/1/5LUSC/5Survival")      #设置工作目录

#绘制生存曲线函数
bioSurvival=function(inputFile=null,outFile=null){
  #读取输入文件
  rt=read.table(inputFile,header=T,sep="\t")
  #比较高低风险组生存差异，得到显著性p值
  diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  
  #绘制生存曲线
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=T,
                     pval=pValue,
                     pval.size=6,
                     risk.table=TRUE,
                     legend.labs=c("High risk", "Low risk"),
                     legend.title="Risk",
                     xlab="Time(years)",
                     break.time.by = 1,
                     risk.table.title="",
                     palette=c("DarkOrchid","Orange2"),
                     risk.table.height=.25,
                     title=paste0("myCAF"))
  pdf(file=outFile,onefile = FALSE,width = 6.5,height =6)
  print(surPlot)
  dev.off()
}
bioSurvival(inputFile="expTime.myCAF.txt",outFile="survival.myCAF.pdf")



bioSurvival=function(inputFile=null,outFile=null){
  #读取输入文件
  rt=read.table(inputFile,header=T,sep="\t")
  #比较高低风险组生存差异，得到显著性p值
  diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  
  #绘制生存曲线
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=T,
                     pval=pValue,
                     pval.size=6,
                     risk.table=TRUE,
                     legend.labs=c("High risk", "Low risk"),
                     legend.title="Risk",
                     xlab="Time(years)",
                     break.time.by = 1,
                     risk.table.title="",
                     palette=c("DarkOrchid","Orange2"),
                     risk.table.height=.25,
                     title=paste0("apCAF"))
  pdf(file=outFile,onefile = FALSE,width = 6.5,height =6)
  print(surPlot)
  dev.off()
}
bioSurvival(inputFile="expTime.apCAF.txt",outFile="survival.apCAF.pdf")



bioSurvival=function(inputFile=null,outFile=null){
  #读取输入文件
  rt=read.table(inputFile,header=T,sep="\t")
  #比较高低风险组生存差异，得到显著性p值
  diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  
  #绘制生存曲线
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=T,
                     pval=pValue,
                     pval.size=6,
                     risk.table=TRUE,
                     legend.labs=c("High risk", "Low risk"),
                     legend.title="Risk",
                     xlab="Time(years)",
                     break.time.by = 1,
                     risk.table.title="",
                     palette=c("DarkOrchid","Orange2"),
                     risk.table.height=.25,
                     title=paste0("iCAF"))
  pdf(file=outFile,onefile = FALSE,width = 6.5,height =6)
  print(surPlot)
  dev.off()
}
bioSurvival(inputFile="expTime.iCAF.txt",outFile="survival.iCAF.pdf")



#引用包
library(survival)
library(survminer)

svdata1=read.table("expTime.myCAF.txt",header=T,sep="\t",row.names = 1)
svdata = svdata1[,c(1,2,length(colnames(svdata1))-1)]
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #默认组内sample不能低于30%
res.cat <- surv_categorize(res.cut)
##统计作图
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl<-list()
group <- res.cat[,"riskScore"] 
survival_dat <- data.frame(group = group)
fit <- survfit(my.surv ~ group)
group <- factor(group, levels = c("low", "high"))
data.survdiff <- survdiff(my.surv ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
if (p.val>0.05) next
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
svsort <- svdata[order(svdata[,"riskScore"]),]
pl[["riskScore"]]<-ggsurvplot(fit, data = survival_dat ,
                              ggtheme = theme_bw(), #想要网格就运行这行
                              conf.int = T, #不画置信区间，想画置信区间就把F改成T
                              #conf.int.style = "step",#置信区间的类型，还可改为ribbon
                              censor = F, #不显示观察值所在的位置
                              palette = c("#DA4A35","#44B1C9"), #线的颜色对应高、低
                              legend.title = "myCAF-riskScore",#基因名写在图例题目的位置
                              font.legend = 11,#图例的字体大小
                              #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                              #在图例上标出高低分界点的表达量，和组内sample数量
                              #在左下角标出pvalue、HR、95% CI
                              #太小的p value标为p < 0.001
                              pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                         paste("p = ",round(p.val,3), sep = "")),
                                           HR, CI, sep = "\n"))
ggsave("survival.myCAF2.pdf",width = 5,height = 5)



#引用包
library(survival)
library(survminer)

svdata1=read.table("expTime.apCAF.txt",header=T,sep="\t",row.names = 1)
svdata = svdata1[,c(1,2,length(colnames(svdata1))-1)]
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #默认组内sample不能低于30%
res.cat <- surv_categorize(res.cut)
##统计作图
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl<-list()
group <- res.cat[,"riskScore"] 
survival_dat <- data.frame(group = group)
fit <- survfit(my.surv ~ group)
group <- factor(group, levels = c("low", "high"))
data.survdiff <- survdiff(my.surv ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
if (p.val>0.05) next
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
svsort <- svdata[order(svdata[,"riskScore"]),]
pl[["riskScore"]]<-ggsurvplot(fit, data = survival_dat ,
                              ggtheme = theme_bw(), #想要网格就运行这行
                              conf.int = T, #不画置信区间，想画置信区间就把F改成T
                              #conf.int.style = "step",#置信区间的类型，还可改为ribbon
                              censor = F, #不显示观察值所在的位置
                              palette = c("#DA4A35","#44B1C9"), #线的颜色对应高、低
                              legend.title = "apCAF-riskScore",#基因名写在图例题目的位置
                              font.legend = 11,#图例的字体大小
                              #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                              #在图例上标出高低分界点的表达量，和组内sample数量
                              #在左下角标出pvalue、HR、95% CI
                              #太小的p value标为p < 0.001
                              pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                         paste("p = ",round(p.val,3), sep = "")),
                                           HR, CI, sep = "\n"))
ggsave("survival.apCAF2.pdf",width = 5,height = 5)



#引用包
library(survival)
library(survminer)

svdata1=read.table("expTime.iCAF.txt",header=T,sep="\t",row.names = 1)
svdata = svdata1[,c(1,2,length(colnames(svdata1))-1)]
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.3) #默认组内sample不能低于30%
res.cat <- surv_categorize(res.cut)
##统计作图
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl<-list()
group <- res.cat[,"riskScore"] 
survival_dat <- data.frame(group = group)
fit <- survfit(my.surv ~ group)
group <- factor(group, levels = c("low", "high"))
data.survdiff <- survdiff(my.surv ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
if (p.val>0.05) next
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
svsort <- svdata[order(svdata[,"riskScore"]),]
pl[["riskScore"]]<-ggsurvplot(fit, data = survival_dat ,
                              ggtheme = theme_bw(), #想要网格就运行这行
                              conf.int = T, #不画置信区间，想画置信区间就把F改成T
                              #conf.int.style = "step",#置信区间的类型，还可改为ribbon
                              censor = F, #不显示观察值所在的位置
                              palette = c("#DA4A35","#44B1C9"), #线的颜色对应高、低
                              legend.title = "iCAF-riskScore",#基因名写在图例题目的位置
                              font.legend = 11,#图例的字体大小
                              #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                              #在图例上标出高低分界点的表达量，和组内sample数量
                              #在左下角标出pvalue、HR、95% CI
                              #太小的p value标为p < 0.001
                              pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                         paste("p = ",round(p.val,3), sep = "")),
                                           HR, CI, sep = "\n"))
ggsave("survival.iCAF2.pdf",width = 5,height = 5)


####5.5TCGA,GEO,相关基因####

#引用包
library(limma)
library(sva)

tcgaExpFile1="TCGA_LUSC.txt" #TCGA表达数据文件
tcgaExpFile2="TCGA_LUAD.txt" #TCGA表达数据文件

geoExpFileGSE30219="GSE30219.txt"
geoExpFileGSE37745="GSE37745.txt"
geoExpFileGSE41271="GSE41271.txt"
geoExpFileGSE68465="GSE68465.txt"
geoExpFileGSE73403="GSE73403.txt"
geoExpFileGSE74777="GSE74777.txt"
geoExpFileGSE157010="GSE157010.txt"

setwd("F:/1/5LUSC/5Survival")  #设置工作目录

#读取TCGA基因表达文件,并对数据进行处理
rt = read.table(tcgaExpFile1, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
tcga1=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
tcga1=avereps(tcga1)
tcga1=log2(tcga1+1)

#读取TCGA基因表达文件,并对数据进行处理
rt = read.table(tcgaExpFile2, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
tcga2=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
tcga2=avereps(tcga2)
tcga2=log2(tcga2+1)

#读取geoGSE41271基因表达文件,并对数据进行处理
rtGSE41271 = read.table(geoExpFileGSE41271,header=T,sep="\t",check.names=F)
rtGSE41271=as.matrix(rtGSE41271)
rownames(rtGSE41271)=rtGSE41271[,1]
expGSE41271=rtGSE41271[,2:ncol(rtGSE41271)]
dimnames=list(rownames(expGSE41271),colnames(expGSE41271))
geoGSE41271=matrix(as.numeric(as.matrix(expGSE41271)),nrow=nrow(expGSE41271),dimnames=dimnames)
geoGSE41271=avereps(geoGSE41271)
#geoGSE41271=log2(geoGSE41271+1)      #需要修改，如果数值很大，去掉前面的#，如果数值很小，保留#
geoGSE41271=normalizeBetweenArrays(geoGSE41271)

#读取geoGSE68465基因表达文件,并对数据进行处理
rtGSE68465 = read.table(geoExpFileGSE68465,header=T,sep="\t",check.names=F)
rtGSE68465=as.matrix(rtGSE68465)
rownames(rtGSE68465)=rtGSE68465[,1]
expGSE68465=rtGSE68465[,2:ncol(rtGSE68465)]
dimnames=list(rownames(expGSE68465),colnames(expGSE68465))
geoGSE68465=matrix(as.numeric(as.matrix(expGSE68465)),nrow=nrow(expGSE68465),dimnames=dimnames)
geoGSE68465=avereps(geoGSE68465)
geoGSE68465=normalizeBetweenArrays(geoGSE68465)
#geoGSE68465=log2(geoGSE68465+1)      #需要修改，如果数值很大，去掉前面的#，如果数值很小，保留#

#读取geoGSE30219基因表达文件,并对数据进行处理
rtGSE30219 = read.table(geoExpFileGSE30219,header=T,sep="\t",check.names=F)
rtGSE30219=as.matrix(rtGSE30219)
rownames(rtGSE30219)=rtGSE30219[,1]
expGSE30219=rtGSE30219[,2:ncol(rtGSE30219)]
dimnames=list(rownames(expGSE30219),colnames(expGSE30219))
geoGSE30219=matrix(as.numeric(as.matrix(expGSE30219)),nrow=nrow(expGSE30219),dimnames=dimnames)
geoGSE30219=avereps(geoGSE30219)
#geoGSE30219=log2(geoGSE30219+1)      #需要修改，如果数值很大，去掉前面的#，如果数值很小，保留#
geoGSE30219=normalizeBetweenArrays(geoGSE30219)

#读取geoGSE73403基因表达文件,并对数据进行处理
rtGSE73403 = read.table(geoExpFileGSE73403,header=T,sep="\t",check.names=F)
rtGSE73403=as.matrix(rtGSE73403)
rownames(rtGSE73403)=rtGSE73403[,1]
expGSE73403=rtGSE73403[,2:ncol(rtGSE73403)]
dimnames=list(rownames(expGSE73403),colnames(expGSE73403))
geoGSE73403=matrix(as.numeric(as.matrix(expGSE73403)),nrow=nrow(expGSE73403),dimnames=dimnames)
geoGSE73403=avereps(geoGSE73403)
#geoGSE73403=log2(geoGSE73403+1)      #需要修改，如果数值很大，去掉前面的#，如果数值很小，保留#
geoGSE73403=normalizeBetweenArrays(geoGSE73403)



#读取geoGSE37745基因表达文件,并对数据进行处理
rtGSE37745 = read.table(geoExpFileGSE37745,header=T,sep="\t",check.names=F)
rtGSE37745=as.matrix(rtGSE37745)
rownames(rtGSE37745)=rtGSE37745[,1]
expGSE37745=rtGSE37745[,2:ncol(rtGSE37745)]
dimnames=list(rownames(expGSE37745),colnames(expGSE37745))
geoGSE37745=matrix(as.numeric(as.matrix(expGSE37745)),nrow=nrow(expGSE37745),dimnames=dimnames)
geoGSE37745=avereps(geoGSE37745)
#geoGSE37745=log2(geoGSE37745+1)      #需要修改，如果数值很大，去掉前面的#，如果数值很小，保留#
geoGSE37745=normalizeBetweenArrays(geoGSE37745)

#读取geoGSE74777基因表达文件,并对数据进行处理
rtGSE74777 = read.table(geoExpFileGSE74777,header=T,sep="\t",check.names=F)
rtGSE74777=as.matrix(rtGSE74777)
rownames(rtGSE74777)=rtGSE74777[,1]
expGSE74777=rtGSE74777[,2:ncol(rtGSE74777)]
dimnames=list(rownames(expGSE74777),colnames(expGSE74777))
geoGSE74777=matrix(as.numeric(as.matrix(expGSE74777)),nrow=nrow(expGSE74777),dimnames=dimnames)
geoGSE74777=avereps(geoGSE74777)
#geoGSE74777=log2(geoGSE74777+1)      #需要修改，如果数值很大，去掉前面的#，如果数值很小，保留#
geoGSE74777=normalizeBetweenArrays(geoGSE74777)


#读取geoGSE157010基因表达文件,并对数据进行处理
rtGSE157010 = read.table(geoExpFileGSE157010,header=T,sep="\t",check.names=F)
rtGSE157010=as.matrix(rtGSE157010)
rownames(rtGSE157010)=rtGSE157010[,1]
expGSE157010=rtGSE157010[,2:ncol(rtGSE157010)]
dimnames=list(rownames(expGSE157010),colnames(expGSE157010))
geoGSE157010=matrix(as.numeric(as.matrix(expGSE157010)),nrow=nrow(expGSE157010),dimnames=dimnames)
geoGSE157010=avereps(geoGSE157010)
#geoGSE157010=log2(geoGSE157010+1)      #需要修改，如果数值很大，去掉前面的#，如果数值很小，保留#
geoGSE157010=normalizeBetweenArrays(geoGSE157010)

#读取基因列表文件, 提取失巢凋亡基因的表达量
gene=read.table("CAFgene.txt", header=F, sep="\t", check.names=F)
#gene=gene[,1][-1]

sameGene1=intersect(gene[,1], row.names(tcga1))
#sameGene1=intersect(row.names(tcga1), row.names(tcga2))
sameGene2=intersect(row.names(tcga2), row.names(geoGSE41271))
sameGene3=intersect(row.names(geoGSE68465), row.names(geoGSE30219))
sameGene4=intersect(row.names(geoGSE73403),row.names(geoGSE37745) )
sameGene5=intersect(row.names(geoGSE74777), row.names(geoGSE157010))

sameGeneA=intersect(sameGene1, sameGene2)
sameGeneB=intersect(sameGene3, sameGene4)
sameGeneA=intersect(sameGeneA, sameGene5)
sameGeneC=intersect(sameGeneA, sameGeneB)

#输出表达量
tcga1=tcga1[sameGeneC,]
tcga2=tcga2[sameGeneC,]
geoGSE41271=geoGSE41271[sameGeneC,]
geoGSE68465=geoGSE68465[sameGeneC,]
geoGSE30219=geoGSE30219[sameGeneC,]
geoGSE73403=geoGSE73403[sameGeneC,]
geoGSE37745=geoGSE37745[sameGeneC,]
geoGSE74777=geoGSE74777[sameGeneC,]
geoGSE157010=geoGSE157010[sameGeneC,]

write.table(data.frame(ID=rownames(tcga1),tcga1,check.names = F),file="diffgene.tcga.LUSC.txt",sep="\t",quote=F,col.names=T,row.names = F)
write.table(data.frame(ID=rownames(tcga2),tcga2,check.names = F),file="diffgene.tcga.LUAD.txt",sep="\t",quote=F,col.names=T,row.names = F)
write.table(data.frame(ID=rownames(geoGSE41271),geoGSE41271),file="diffgene.geoGSE41271.txt",sep="\t",quote=F,col.names=T,row.names = F)
write.table(data.frame(ID=rownames(geoGSE68465),geoGSE68465),file="diffgene.geoGSE68465.txt",sep="\t",quote=F,col.names=T,row.names = F)
write.table(data.frame(ID=rownames(geoGSE30219),geoGSE30219),file="diffgene.geoGSE30219.txt",sep="\t",quote=F,col.names=T,row.names = F)
write.table(data.frame(ID=rownames(geoGSE73403),geoGSE73403),file="diffgene.geoGSE73403.txt",sep="\t",quote=F,col.names=T,row.names = F)
write.table(data.frame(ID=rownames(geoGSE37745),geoGSE37745),file="diffgene.geoGSE37745.txt",sep="\t",quote=F,col.names=T,row.names = F)
write.table(data.frame(ID=rownames(geoGSE74777),geoGSE74777),file="diffgene.geoGSE74777.txt",sep="\t",quote=F,col.names=T,row.names = F)
write.table(data.frame(ID=rownames(geoGSE157010),geoGSE157010),file="diffgene.geoGSE157010.txt",sep="\t",quote=F,col.names=T,row.names = F)






####5.6单因素cox####


#引用包
library(limma)
library(survival)	
library(survminer)

expFile="diffgene.tcga.LUSC.txt"      #表达数据文件
cliFile="time_LUSC.txt"            #生存数据文件
setwd("F:/1/5LUSC/5Survival")    #设置工作目录

#读取表达文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data2=data
data=data[rowMeans(data)>2,]
data=t(data)
dim(data)

#读取生存数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli$futime=cli$futime/365


#markers <- read_csv("markers0.1.0.3.0.1.csv")
#apCAF = subset(markers,cluster == "apCAF")

#samegeneapCAF = intersect(colnames(data),as.list(apCAF[,1])[["...1"]])
#data = data[,samegeneapCAF]

#rownames(data) = gsub('GSE116918_', '', rownames(data))
#rownames(data) = gsub('TCGA_', '', rownames(data))

#数据合并
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",row.names(data))
sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
rt=cbind(cli, data)

#对基因进行循环，查找预后相关的基因
outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  #cox分析
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<0.05){
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}

#输出单因素的结果
write.table(outTab, file="uniCoxCAF.txt", sep="\t", row.names=F, quote=F)

#输出单因素显著基因的表达量
sigGenes=as.vector(outTab[,1])
uniSigExp=data2[sigGenes,]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExpCAF.txt",sep="\t",row.names=F,quote=F)

#输出单因素显著基因表达数据和生存数据合并的文件
uniSigExpTime=rt[,c("futime","fustat",sigGenes)]
uniSigExpTime=cbind(id=row.names(uniSigExpTime),uniSigExpTime)
write.table(uniSigExpTime, file="uniSigExpTimeCAF.txt", sep="\t", row.names=F, quote=F)

############定义森林图函数
bioForest=function(coxFile=null,forestFile=null,forestCol=null){
  #读取输入文件
  rt <- read.table(coxFile,header=T,sep="\t",check.names=F,row.names=1)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #输出图形
  height=nrow(rt)/13+7
  pdf(file=forestFile, width = 7,height = height)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #绘制森林图左边的基因信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  
  #绘制森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol[1], forestCol[2])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)
  axis(1)
  dev.off()
}

#调用函数, 对单因素的结果进行可视化, 绘制森林图
bioForest(coxFile="uniCoxCAF.txt", forestFile="forestCAF.pdf", forestCol=c("#DA4A35","#44B1C9"))



####6.1CNV####
####6.2CNVfreq####


inputFile="cnvMatrix.txt"     #输入文件
setwd("F:/1/5LUSC/6Mutation")   #设置工作目录

rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    #读取输入文件
GAIN=rowSums(rt> 0)         #获取拷贝数增加的样品数目
LOSS=rowSums(rt< 0)         #获取拷贝数缺失的样品数目
GAIN=GAIN/ncol(rt)*100      #计算拷贝数增加的百分率
LOSS=LOSS/ncol(rt)*100      #计算拷贝数缺失的百分率
data=cbind(GAIN, LOSS)
data=data[order(data[,"GAIN"],decreasing = T),]

#绘制图形
data.max = apply(data, 1, max)
pdf(file="CNVfreq.pdf", width=9, height=5)
cex=1.3
par(cex.lab=cex, cex.axis=cex, font.axis=2, las=1, xpd=T)
bar=barplot(data.max, col="grey80", border=NA,
            xlab="", ylab="CNV.frequency(%)", space=1.5,
            xaxt="n", ylim=c(0,1.2*max(data.max)))
points(bar,data[,"GAIN"], pch=20, col="DarkOrchid", cex=3)
points(bar,data[,"LOSS"], pch=20, col="Orange2", cex=3)
legend("top", legend=c('GAIN','LOSS'), col=c("DarkOrchid","Orange2"), pch=20, bty="n", cex=2, ncol=2)
par(srt=45)
text(bar, par('usr')[3]-0.2, rownames(data), adj=1, cex=0.8)
dev.off()


####6.3preRcircos####
####6.4Rcircos####


library("RCircos")       #引用包
setwd("F:/1/5LUSC/6Mutation")   #设置工作目录

#初始化圈图
cytoBandIdeogram=read.table("refer.txt", header=T, sep="\t")
chr.exclude <- NULL
cyto.info <- cytoBandIdeogram
tracks.inside <- 5
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

#设置圈图参数
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$text.size=0.6
rcircos.params$point.size=5
RCircos.Reset.Plot.Parameters(rcircos.params)

#输出文件
pdf(file="RCircos.pdf", width=7, height=7)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

#读取拷贝数的文件，绘制散点图
RCircos.Scatter.Data=read.table("Rcircos.scatter.txt", header=T, sep="\t", check.names=F)
data.col <- 4
track.num <- 1
side <- "in"
RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col, track.num, side, by.fold=0.1)

#读取基因注释文件，标注基因的名称
RCircos.Gene.Label.Data=read.table("Rcircos.geneLabel.txt", header=T, sep="\t", check.names=F)
name.col <- 4
side <- "in"
track.num <- 2
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side)
track.num <- 3
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col, track.num, side)
dev.off()


####6.5network####


#引用包
library(igraph)
library(psych)
library(reshape2)
library(RColorBrewer)

GeneExpfile="diffgene.tcga.LUSC.txt"      #表达数据文件
Coxfile="uniCox.txt"              #生存分析的结果文件
setwd("F:/1/5LUSC/6Mutation")   #设置工作目录

#读取输入文件
gene.exp <- read.table(GeneExpfile,header=T,sep="\t",row.names=1)
gene.cox <- read.table(Coxfile,header=T,sep="\t")

#基因取交集
gene.group=data.frame(id=gene.cox[,1], group="CAFRG")
genelist <- intersect(gene.group$id, gene.cox$id)
genelist <- intersect(genelist, rownames(gene.exp))
gene.group <- gene.group[match(genelist,gene.group$id),]
gene.group <- gene.group[order(gene.group$group),]
gene.exp <- gene.exp[match(gene.group$id,rownames(gene.exp)),]
gene.cox <- gene.cox[match(gene.group$id,gene.cox$id),]

#准备网络关系文件
gene.cor <- corr.test(t(gene.exp))
gene.cor.cor <- gene.cor$r
gene.cor.pvalue <- gene.cor$p
gene.cor.cor[upper.tri(gene.cor.cor)] = NA
gene.cor.pvalue[upper.tri(gene.cor.pvalue)] = NA
gene.cor.cor.melt <- melt(gene.cor.cor)   #gene1 \t gene2 \t cor
gene.cor.pvalue.melt <- melt(gene.cor.pvalue)
gene.melt <- data.frame(from = gene.cor.cor.melt$Var2,to=gene.cor.cor.melt$Var1,cor=gene.cor.cor.melt$value,pvalue=gene.cor.pvalue.melt$value)
gene.melt <- gene.melt[gene.melt$from!=gene.melt$to&!is.na(gene.melt$pvalue),,drop=F]
gene.edge <- gene.melt[gene.melt$pvalue<0.0001,,drop=F]
gene.edge$color <- ifelse(gene.edge$cor>0,'pink','#6495ED')
gene.edge$weight <- abs(gene.edge$cor)*4

#准备节点属性文件
gene.node <- gene.group
group.color <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(gene.node$group)))
gene.node$color <- group.color[as.numeric(as.factor(gene.node$group))]
gene.node$shape <- "circle"
gene.node$frame <- ifelse(gene.cox$HR>1,"DarkOrchid","Orange2")
gene.node$pvalue <- gene.cox$pvalue
#节点的大小
pvalue.breaks <- c(0,0.0001,0.001,0.01,0.05,1)
pvalue.size <- c(16,14,12,10,8)
cutpvalue <- cut(gene.node$pvalue,breaks=pvalue.breaks)
gene.node$size <- pvalue.size[as.numeric(cutpvalue)]
nodefile <- "network.node.txt"
edgefile <- "network.edge.txt"
write.table(gene.node, nodefile, sep="\t", col.names=T, row.names=F, quote=F)
write.table(gene.edge, edgefile, sep="\t", col.names=T, row.names=F, quote=F)


#绘制网络图
node = read.table(nodefile, header=T, sep="\t", comment.char="")
edge = read.table(edgefile, header=T, sep="\t", comment.char="")

g = graph.data.frame(edge,directed = FALSE)
node = node[match(names(components(g)$membership),node$id),]

if(!is.na(match('color',colnames(node)))) V(g)$color = node$color
if(!is.na(match('size',colnames(node)))) V(g)$size = node$size
if(!is.na(match('shape',colnames(node)))) V(g)$shape = node$shape
if(!is.na(match('frame',colnames(node)))) V(g)$frame = node$frame

#输出图形
pdf(file="network.pdf", width=9.5, height=7.5)
par(mar=c(0,0,0,0))
layout(matrix(c(1,1,4,2,3,4),nc=2),height=c(4,4,2),width=c(8,3))

#节点坐标 
coord = layout_in_circle(g)
degree.x = acos(coord[,1])
degree.y = asin(coord[,2])
degree.alpha = c()
for(i in 1:length(degree.x)){
  if(degree.y[i]<0) degree.alpha=c(degree.alpha,2*pi-degree.x[i]) else degree.alpha=c(degree.alpha,degree.x[i])
}
degree.cut.group = (0:8)/4*pi
degree.cut.group[1] = -0.0001
degree.cut = cut(degree.alpha,degree.cut.group)
#degree.degree = c(-pi/4,-pi/4,-pi/2,-pi/2,pi/2,pi/2,pi/2,pi/4)
degree.degree = c(-pi/6,-pi/4,-pi/2,-pi/2,pi/2,pi/2,pi/2,pi/6)*1.45
degree = degree.degree[as.numeric(degree.cut)]

#定义饼图,左半圆颜色代表基因的属性,右半圆代表基因的风险,哪些基因是高风险基因,哪些基因是低风险基因
values <- lapply(node$id,function(x)c(1,1))
V(g)$pie.color = lapply(1:nrow(node),function(x)c(node$color[x],node$frame[x]))
V(g)$frame = NA 

#绘制图形的主体部分
plot(g,layout=layout_in_circle,vertex.shape="pie",vertex.pie=values,
     vertex.label.cex=V(g)$lable.cex,edge.width = E(g)$weight,edge.arrow.size=0,
     vertex.label.color=V(g)$color,vertex.frame.color=V(g)$frame,edge.color=E(g)$color,
     vertex.label.cex=2,vertex.label.font=2,vertex.size=V(g)$size,edge.curved=0.5,
     vertex.color=V(g)$color,vertex.label.dist=1.35,vertex.label.degree=degree)
# label.degree : zero means to the right; and pi means to the left; up is -pi/2 and down is pi/2;  The default value is -pi/4
# label.dist If it is 0 then the label is centered on the vertex; If it is 1 then the label is displayed beside the vertex.

#绘制节点属性图例(基因的属性)
par(mar=c(0,0,0,0))
plot(1,type="n",xlab="",ylab="",axes=F)
groupinfo = unique(data.frame(group=node$group,color=node$color))
legend("left",legend=groupinfo$group,col=groupinfo$color,pch=16,bty="n",cex=3)
#绘制基因风险的图例(哪些基因是高风险的基因,哪些基因是低风险的基因)
par(mar=c(0,0,0,0))
plot(1,type="n",xlab="",ylab="",axes=F)
legend("left",legend=c('Risk factors','Favorable factors'),col=c("DarkOrchid","Orange2"),pch=16,bty="n",cex=2.5)
#绘制预后pvalue的图例
par(mar=c(0,0,0,0))
plot(1,type="n",xlab="",axes=F,ylab="")
legend("top",legend=c('Postive correlation with P<0.0001','Negative correlation with P<0.0001'),lty=1,lwd=4,col=c('pink','#6495ED'),bty="n",cex=2.2)
legend('bottom',legend=c(0.0001,0.001,0.01,0.05,1),pch=16,pt.cex=c(1.6,1.4,1.2,1,0.8)*6,bty="n",ncol=5,cex=2.2,col="black",title="Cox test, pvalue")
dev.off()


####7.1cluster####


#引用包
library(limma)
library(survival)
library(ConsensusClusterPlus)

expFile="uniSigExpTimeCAF.txt"      #表达数据文件
workDir="F:/1/5LUSC/7Cluster"     #工作目录
setwd(workDir)       #设置工作目录

#读取输入文件
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)
data = t(data[,-c(1,2)])


#对样品进行分型
#pam,km
#pearson,spearman,euclidean
maxK=9      #最大的k值(最多可以将样品分成几个组)
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=100,
                             pItem=0.8,
                             pFeature=1,
                             title=workDir,
                             clusterAlg="km",
                             distance="euclidean",
                             seed=123,
                             plot="png")

#输出分型结果
clusterNum=3      #分成几个亚型
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("CAFcluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$CAFcluster))
cluster$CAFcluster=letter[match(cluster$CAFcluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="CAFcluster.txt", sep="\t", quote=F, col.names=F)




####7.2Sur####

#引用包
library(survival)
library(survminer)

clusterFile="CAFcluster.txt"     #分型的结果文件
cliFile="time_LUSC.txt"               #生存数据文件
setwd("F:/1/5LUSC/7Cluster")     #设置工作目录

#读取分型的结果文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
rownames(cluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(cluster))
#读取生存数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365

#数据合并
sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])

#生存差异分析,观察不同分型之间病人的生存是否具有差异
length=length(levels(factor(rt$CAFcluster)))
diff=survdiff(Surv(futime, fustat) ~ CAFcluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ CAFcluster, data = rt)
#print(surv_median(fit))

#绘制生存曲线
bioCol=c("#DA4A35","#44B1C9","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="CAFcluster",
                   legend.labs=levels(factor(rt[,"CAFcluster"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 2,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)

#输出图形
pdf(file="survival.pdf", width=6.5, height=6.25, onefile=FALSE)
print(surPlot)
dev.off()







####7.3PCA####


#引用包
library(limma)
library(Rtsne)
library(umap)
library(ggplot2)

expFile="uniSigExpTimeCAF.txt"           #表达数据文件
clusterFile="CAFcluster.txt"      #分型的结果文件
setwd("F:/1/5LUSC/7Cluster")     #设置工作目录

#读取输入文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)

rt = rt[,-c(2,3)]

rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#data=t(data)

############PCA分析
data.pca=prcomp(data, scale. = TRUE)
pcaPredict=predict(data.pca)

#读取分型的结果文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
CAFcluster=as.vector(cluster[,1])

#设置分型的颜色
bioCol=c("#DA4A35","#44B1C9","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
clusterCol=bioCol[1:length(levels(factor(CAFcluster)))]

#绘制图形
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], CAFcluster=CAFcluster)
PCA.mean=aggregate(PCA[,1:2], list(CAFcluster=PCA$CAFcluster), mean)
pdf(file="PCA.pdf", width=5.5, height=4.75)
p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color=CAFcluster, shape=CAFcluster)) +
  scale_colour_manual(name="CAFcluster", values =clusterCol)+
  theme_bw()+
  labs(title ="PCA")+
  theme(plot.margin=unit(rep(1.5,4),'lines'), plot.title = element_text(hjust=0.5))+
  annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$CAFcluster, cex=7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()


############UMAP分析
umapOut=umap(data)
umap=data.frame(UMAP1=umapOut$layout[,1], UMAP2=umapOut$layout[,2], CAFcluster=CAFcluster)
umap.mean=aggregate(umap[,1:2], list(CAFcluster=umap$CAFcluster), mean)	
#绘制分型的UAMP图
pdf(file="UMAP.pdf", width=5.5, height=4.75)       #保存输入出文件
p=ggplot(data=umap, aes(UMAP1, UMAP2)) + geom_point(aes(color=CAFcluster, shape=CAFcluster)) +
  scale_colour_manual(name="CAFcluster",  values =clusterCol)+
  theme_bw()+
  labs(title ="UAMP")+
  theme(plot.margin=unit(rep(1.5,4),'lines'), plot.title = element_text(hjust=0.5))+
  annotate("text",x=umap.mean$UMAP1, y=umap.mean$UMAP2, label=umap.mean$CAFcluster, cex=7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()


############tSNE分析
tsneOut=Rtsne(data, dims=2, perplexity=10, verbose=F, max_iter=500,check_duplicates=F)
tsne=data.frame(tSNE1=tsneOut$Y[,1], tSNE2=tsneOut$Y[,2], CAFcluster=CAFcluster)
tSNE.mean=aggregate(tsne[,1:2], list(CAFcluster=tsne$CAFcluster), mean)	
#绘制分型的tSNE图
pdf(file="tSNE.pdf", width=5.5, height=4.75)       #保存输入出文件
p=ggplot(data = tsne, aes(tSNE1, tSNE2)) + geom_point(aes(color=CAFcluster, shape=CAFcluster)) +
  scale_colour_manual(name="CAFcluster",  values =clusterCol)+
  theme_bw()+
  labs(title ="tSNE")+
  theme(plot.margin=unit(rep(1.5,4),'lines'), plot.title = element_text(hjust=0.5))+
  annotate("text",x=tSNE.mean$tSNE1, y=tSNE.mean$tSNE2, label=tSNE.mean$CAFcluster, cex=7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()




####7.4clusterDiff####

#引用包
library(limma)
library(reshape2)
library(ggpubr)

expFile="uniSigExpCAF.txt"      #表达数据文件
cluFile="CAFcluster.txt"     #分型的结果文件
setwd("F:/1/5LUSC/7Cluster")     #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
colnames(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(exp))
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)

#读取分型的结果文件
cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
sameSample=intersect(row.names(data), row.names(cluster))
expClu=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])

#提取差异显著的基因
sigGene=c()
for(i in colnames(expClu)[1:(ncol(expClu)-1)]){
  if(sd(expClu[,i])<0.001){next}
  if(length(levels(factor(expClu[,"CAFcluster"])))>2){
    test=kruskal.test(expClu[,i] ~ expClu[,"CAFcluster"])
  }else{
    test=wilcox.test(expClu[,i] ~ expClu[,"CAFcluster"])
  }
  pvalue=test$p.value
  if(pvalue<0.05){
    sigGene=c(sigGene, i)
  }
}
sigGene=c(sigGene, "CAFcluster")
expClu=expClu[,sigGene]

#把数据转换成ggplot2输入文件
data=melt(expClu, id.vars=c("CAFcluster"))
colnames(data)=c("CAFcluster", "Gene", "Expression")

#设置图形颜色
bioCol=c("#DA4A35","#44B1C9","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"CAFcluster"])))]

#绘制箱线图
p=ggboxplot(data, x="Gene", y="Expression", color="CAFcluster",
            xlab="",
            ylab="Gene expression",
            legend.title="CAFcluster",
            palette = bioCol,
            width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=CAFcluster),
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

#输出箱线图
pdf(file="geneboxplot.diff.pdf", width=10, height=4.5)
print(p1)
dev.off()






####7.5heatmap####

library(pheatmap)       #引用包
expFile="uniSigExpCAF.txt"          #表达数据文件
clusterFile="CAFcluster.txt"     #分型的结果文件
cliFile="clinical5.txt"           #临床数据文件
setwd("F:/1/5LUSC/7Cluster")  #设置工作目录

#读取表达数据文件
exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(exp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(exp))
exp=t(exp)
#读取分型的结果文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#合并表达和分型数据
sameSample=intersect(row.names(exp), row.names(cluster))
exp=exp[sameSample, , drop=F]
cluster=cluster[sameSample, , drop=F]
expCluster=cbind(exp, cluster)
#Project=gsub("(.*?)\\_.*", "\\1", rownames(expCluster))
rownames(expCluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(expCluster))
#expCluster=cbind(expCluster, Project)

#合并临床数据
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))
sameSample=intersect(row.names(expCluster), row.names(cli))
expCluster=expCluster[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
data=cbind(expCluster, cli)

#提取热图的数据
data=data[order(data$CAFcluster),]
Type=data[,((ncol(exp)+1):ncol(data))]
data=t(data[,1:ncol(exp)])

#定义热图注释的颜色
bioCol=c("#DA4A35","#44B1C9","#7CC767","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
prgCluCol=bioCol[1:length(levels(factor(Type$CAFcluster)))]
names(prgCluCol)=levels(factor(Type$CAFcluster))
ann_colors[["CAFcluster"]]=prgCluCol

#热图可视化
pdf("heatmap.pdf", width=8, height=6)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("#44B1C9",5), "white", rep("#DA4A35",5)))(100),
         cluster_cols =F,
         cluster_rows =F,
         show_colnames=F,
         scale="row",
         fontsize=6,
         fontsize_row=6,
         fontsize_col=6)
dev.off()




####7.6ssGSEA####

#引用包
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)

expFile="TCGA_LUSC.txt"               #表达数据文件
clusterFile="CAFcluster.txt"      #分型的结果文件
gmtFile="immune.gmt"              #免疫基因集文件
setwd("F:/1/5LUSC/7Cluster")     #设置工作目录

#读取表达数据文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#读取基因集文件
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#ssGSEA分析
ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#对ssGSEA打分进行矫正
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
#输出ssGSEA打分结果
ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="ssGSEA.result.txt",sep="\t",quote=F,col.names=F)

#读取分型的结果文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并(将分型的结果与免疫细胞的打分进行合并)
ssgseaScore=t(ssgseaScore)
row.names(ssgseaScore)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",row.names(ssgseaScore))
sameSample=intersect(row.names(ssgseaScore), row.names(cluster))
ssgseaScore=ssgseaScore[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
scoreCluster=cbind(ssgseaScore, cluster)

#把数据转换成ggplot2输入文件
data=melt(scoreCluster, id.vars=c("CAFcluster"))
colnames(data)=c("CAFcluster", "Immune", "Fraction")

#绘制箱线图
bioCol=c("#DA4A35","#44B1C9","#7CC767","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"CAFcluster"])))]
p=ggboxplot(data, x="Immune", y="Fraction", color="CAFcluster",
            xlab="",
            ylab="Immune infiltration",
            legend.title="CAFcluster",
            palette=bioCol)
p=p+rotate_x_text(50)

#保存图形
pdf(file="boxplot.pdf", width=11, height=5.2)
p+stat_compare_means(aes(group=CAFcluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")
dev.off()





####7.7GSVA####

#引用包
library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)

expFile="TCGA_LUSC.txt"                 #表达数据文件
clusterFile="CAFcluster.txt"        #分型的结果文件
gmtFile="c2.cp.kegg.symbols.gmt"    #基因集文件
setwd("F:/1/5LUSC/7Cluster")     #设置工作目录

#读取表达输入文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#GSVA分析
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())
gsvaResult=gsva(data, 
                geneSets, 
                min.sz=10, 
                max.sz=500, 
                verbose=TRUE,
                parallel.sz=1)
gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)
write.table(gsvaOut, file="gsvaOut.txt", sep="\t", quote=F, col.names=F)

#读取分型的结果
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并(将分型和GSVA的结果进行合并)
gsvaResult=t(gsvaResult)
row.names(gsvaResult)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",row.names(gsvaResult))
sameSample=intersect(row.names(gsvaResult), row.names(cluster))
gsvaResult=gsvaResult[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
gsvaCluster=cbind(gsvaResult, cluster)
#Project=gsub("(.*?)\\_.*", "\\1", rownames(gsvaCluster))
#gsvaCluster=cbind(gsvaCluster, Project)

#设置比较组
adj.P.Val.Filter=0.05
allType=as.vector(gsvaCluster$CAFcluster)
comp=combn(levels(factor(allType)), 2)

#对比较组进行循环, 观察哪些通路在不同的分型之间是具有差异的
for(i in 1:ncol(comp)){
  #样品分组
  treat=gsvaCluster[gsvaCluster$CAFcluster==comp[2,i],]
  con=gsvaCluster[gsvaCluster$CAFcluster==comp[1,i],]
  data=rbind(con, treat)
  #对通路进行差异分析
  Type=as.vector(data$CAFcluster)
  ann=data[,c(ncol(data), (ncol(data)-1))]
  data=t(data[,-c((ncol(data)-1), ncol(data))])
  design=model.matrix(~0+factor(Type))
  colnames(design)=levels(factor(Type))
  fit=lmFit(data, design)
  contrast=paste0(comp[2,i], "-", comp[1,i])
  cont.matrix=makeContrasts(contrast, levels=design)
  fit2=contrasts.fit(fit, cont.matrix)
  fit2=eBayes(fit2)
  
  #输出所有通路的差异情况
  allDiff=topTable(fit2,adjust='fdr',number=200000)
  allDiffOut=rbind(id=colnames(allDiff),allDiff)
  write.table(allDiffOut, file=paste0(contrast, ".all.txt"), sep="\t", quote=F, col.names=F)
  
  #输出差异显著的通路
  diffSig=allDiff[with(allDiff, (abs(logFC)>0.1 & adj.P.Val < adj.P.Val.Filter )), ]
  diffSigOut=rbind(id=colnames(diffSig),diffSig)
  write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)
  
  #设置热图注释的颜色
  bioCol=c("#DA4A35","#44B1C9","#7CC767","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
  ann_colors=list()
  m6aCluCol=bioCol[1:length(levels(factor(allType)))]
  names(m6aCluCol)=levels(factor(allType))
  ann_colors[["CAFcluster"]]=m6aCluCol[c(comp[1,i], comp[2,i])]
  
  #绘制差异通路热图
  termNum=30     #设置显示通路的数目
  diffTermName=as.vector(rownames(diffSig))
  diffLength=length(diffTermName)
  if(diffLength<termNum){termNum=diffLength}
  hmGene=diffTermName[1:termNum]
  hmExp=data[hmGene,]
  pdf(file=paste0(contrast,".heatmap.pdf"), width=10, height=5)
  pheatmap(hmExp, 
           annotation=ann,
           annotation_colors = ann_colors,
           color = colorRampPalette(c(rep("#44B1C9",2), "white", rep("#DA4A35",2)))(50),
           cluster_cols =F,
           show_colnames = F,cluster_rows = T,
           gaps_col=as.vector(cumsum(table(Type))),
           scale="row",
           fontsize = 8,
           fontsize_row=6,
           fontsize_col=8)
  dev.off()
}




####8.1合并生存时间####


setwd("F:/1/5LUSC/8Model")   #设置工作目录
library(limma)       #引用包
expFile="diffgene.tcga.LUSC.txt"      #表达数据文件
cliFile="time_LUSC.txt"   #临床数据

#读取表达文件，并对输入文件整理
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
#data=t(data)

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(data))
data=t(data)
data=avereps(data)

#读取生存数据
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)     #读取临床文件

#数据合并并输出结果
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
write.table(out,file="expTime.tcga.LUSC.txt",sep="\t",row.names=F,quote=F)



library(limma)       #引用包
expFile="diffgene.tcga.LUAD.txt"      #表达数据文件
cliFile="time_LUAD.txt"   #临床数据

#读取表达文件，并对输入文件整理
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
#data=t(data)

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(data))
data=t(data)
data=avereps(data)

#读取生存数据
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)     #读取临床文件

#数据合并并输出结果
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
write.table(out,file="expTime.tcga.LUAD.txt",sep="\t",row.names=F,quote=F)
#gene = colnames(data)



library(limma)       #引用包
expFile="diffgene.geoGSE30219.txt"       #表达数据文件
cliFile="timeGSE30219.txt"   #临床数据

#读取表达文件，并对输入文件整理
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#data=data[rowMeans(data)>0,]
data=t(data)

#读取生存数据
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)     #读取临床文件

#数据合并并输出结果
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
#data=data[,TCGAgene]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)

write.table(out,file="expTime.GSE30219.txt",sep="\t",row.names=F,quote=F)



library(limma)       #引用包
expFile="diffgene.geoGSE37745.txt"       #表达数据文件
cliFile="timeGSE37745.txt"   #临床数据

#读取表达文件，并对输入文件整理
rt=read.table(expFile,sep="\t",header=T,check.names=F)

rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#data=data[rowMeans(data)>0,]
data=t(data)

#读取生存数据
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)     #读取临床文件

#数据合并并输出结果
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
#data=data[,TCGAgene]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
write.table(out,file="expTime.GSE37745.txt",sep="\t",row.names=F,quote=F)



library(limma)       #引用包
expFile="diffgene.geoGSE41271.txt"       #表达数据文件
cliFile="timeGSE41271.txt"   #临床数据

#读取表达文件，并对输入文件整理
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#data=data[rowMeans(data)>0,]
data=t(data)

#读取生存数据
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)     #读取临床文件

#数据合并并输出结果
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
#data=data[,TCGAgene]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
write.table(out,file="expTime.GSE41271.txt",sep="\t",row.names=F,quote=F)



library(limma)       #引用包
expFile="diffgene.geoGSE68465.txt"       #表达数据文件
cliFile="timeGSE68465.txt"   #临床数据

#读取表达文件，并对输入文件整理
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#data=data[rowMeans(data)>0,]
data=t(data)

#读取生存数据
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)     #读取临床文件

#数据合并并输出结果
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
#data=data[,TCGAgene]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
write.table(out,file="expTime.GSE68465.txt",sep="\t",row.names=F,quote=F)



library(limma)       #引用包
expFile="diffgene.geoGSE73403.txt"       #表达数据文件
cliFile="timeGSE73403.txt"   #临床数据

#读取表达文件，并对输入文件整理
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#data=data[rowMeans(data)>0,]
data=t(data)

#读取生存数据
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)     #读取临床文件

#数据合并并输出结果
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
#data=data[,TCGAgene]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
write.table(out,file="expTime.GSE73403.txt",sep="\t",row.names=F,quote=F)



library(limma)       #引用包
expFile="diffgene.geoGSE74777.txt"       #表达数据文件
cliFile="timeGSE74777.txt"   #临床数据

#读取表达文件，并对输入文件整理
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#data=data[rowMeans(data)>0,]
data=t(data)

#读取生存数据
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)     #读取临床文件

#数据合并并输出结果
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
#data=data[,TCGAgene]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
write.table(out,file="expTime.GSE74777.txt",sep="\t",row.names=F,quote=F)



library(limma)       #引用包
expFile="diffgene.geoGSE157010.txt"       #表达数据文件
cliFile="timeGSE157010.txt"   #临床数据

#读取表达文件，并对输入文件整理
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#data=data[rowMeans(data)>0,]
data=t(data)

#读取生存数据
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)     #读取临床文件

#数据合并并输出结果
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
#data=data[,TCGAgene]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
write.table(out,file="expTime.GSE157010.txt",sep="\t",row.names=F,quote=F)



####8.2model####


#引用包
library(survival)
library(caret)
library(glmnet)
library(survminer)
library(timeROC)

setwd("F:/1/5LUSC/8Model")   #设置工作目录
rt=read.table("uniSigExpTimeCAF.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件
rt$futime[rt$futime<=0]=0.003

#对数据进行分组, 构建模型
n=100000     #分组的数目
for(i in 1:n){
  #############对数据进行分组
  inTrain<-createDataPartition(y=rt[,2], p=0.7, list=F)
  train<-rt[inTrain,]
  test<-rt[-inTrain,]
  trainOut=cbind(id=row.names(train),train)
  testOut=cbind(id=row.names(test),test)
  
  #lasso回归分析
  x=as.matrix(train[,c(3:ncol(train))])
  y=data.matrix(Surv(train$futime,train$fustat))
  fit <- glmnet(x, y, family = "cox", maxit = 1000)
  cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
  coef <- coef(fit, s = cvfit$lambda.min)
  index <- which(coef != 0)
  actCoef <- coef[index]
  lassoGene=row.names(coef)[index]
  lassoSigExp=train[,c("futime", "fustat", lassoGene)]
  lassoSigExpOut=cbind(id=row.names(lassoSigExp), lassoSigExp)
  geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
  if(nrow(geneCoef)<2){next}
  
  #############构建COX模型
  multiCox <- coxph(Surv(futime, fustat) ~ ., data = lassoSigExp)
  #multiCox=step(multiCox,direction = "both")
  multiCoxSum=summary(multiCox)
  
  #输出模型的公式
  outMultiTab=data.frame()
  outMultiTab=cbind(
    coef=multiCoxSum$coefficients[,"coef"],
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
  outMultiTab=outMultiTab[,1:2]
  
  #输出train组风险文件
  riskScore=predict(multiCox,type="risk",newdata=train)      #利用train得到模型预测train样品风险
  coxGene=rownames(multiCoxSum$coefficients)
  coxGene=gsub("`","",coxGene)
  outCol=c("futime","fustat",coxGene)
  medianTrainRisk=median(riskScore)
  risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
  trainRiskOut=cbind(id=rownames(cbind(train[,outCol],riskScore,risk)),cbind(train[,outCol],riskScore,risk))
  
  #输出test组风险文件
  riskScoreTest=predict(multiCox,type="risk",newdata=test)     #利用train得到模型预测test样品风险
  riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
  testRiskOut=cbind(id=rownames(cbind(test[,outCol],riskScoreTest,riskTest)),cbind(test[,outCol],riskScore=riskScoreTest,risk=riskTest))
  
  #计算train组和test组中高低风险组生存差异的pvalue	
  diff=survdiff(Surv(futime, fustat) ~risk,data = train)
  pValue=1-pchisq(diff$chisq, df=1)
  diffTest=survdiff(Surv(futime, fustat) ~riskTest,data = test)
  pValueTest=1-pchisq(diffTest$chisq, df=1)
  
  #ROC曲线下面积
  predictTime=5   #预测的时间
  roc=timeROC(T=train$futime, delta=train$fustat,
              marker=riskScore, cause=1,
              times=c(predictTime), ROC=TRUE)
  rocTest=timeROC(T=test$futime, delta=test$fustat,
                  marker=riskScoreTest, cause=1,
                  times=c(predictTime), ROC=TRUE)	
  
  if((pValue<0.05) & (roc$AUC[2]>0.73) & (pValueTest<0.05) & (rocTest$AUC[2]>0.7)){
    #输出分组结果
    write.table(trainOut,file="data.train.txt",sep="\t",quote=F,row.names=F)
    write.table(testOut,file="data.test.txt",sep="\t",quote=F,row.names=F)
    #lasso结果
    write.table(lassoSigExpOut,file="lasso.SigExp.txt",sep="\t",row.names=F,quote=F)
    pdf("lasso.lambda.pdf")
    plot(fit, xvar = "lambda", label = TRUE)
    dev.off()
    pdf("lasso.cvfit.pdf")
    plot(cvfit)
    abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
    dev.off()
    #输出多因素结果
    write.table(outMultiTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)
    write.table(trainRiskOut,file="risk.train.txt",sep="\t",quote=F,row.names=F)
    write.table(testRiskOut,file="risk.test.txt",sep="\t",quote=F,row.names=F)
    #所有样品的风险文件
    allRiskOut=rbind(trainRiskOut, testRiskOut)
    write.table(allRiskOut,file="risk.all.txt",sep="\t",quote=F,row.names=F)
    break
  }
}



####8.3外部训练集####

#引用包
library(glmnet)
library(survival)
setwd("F:/1/5LUSC/8Model")
trainFile="lasso.SigExp.txt"      #train组输入文件
testFile1tcga.LUAD="expTime.tcga.LUAD.txt"  #test组输入文件
testFile1GSE30219="expTime.GSE30219.txt"  #test组输入文件
testFile1GSE37745="expTime.GSE37745.txt"  #test组输入文件
testFile1GSE41271="expTime.GSE41271.txt"  #test组输入文件
testFile1GSE68465="expTime.GSE68465.txt"  #test组输入文件
testFile1GSE73403="expTime.GSE73403.txt"  #test组输入文件
testFile1GSE74777="expTime.GSE74777.txt"  #test组输入文件
testFile1GSE157010="expTime.GSE157010.txt"  #test组输入文件

rt=read.table(trainFile, header=T, sep="\t", row.names=1,check.names=F)    #读取train组输入文件
#COX模型构建
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
#multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)
#输出模型相关信息
outMultiTab=data.frame()
outMultiTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
write.table(outMultiTab,file="multiCox2.txt",sep="\t",row.names=F,quote=F)
#输出train组风险值
trainScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
quantile(trainScore,seq(0.1,1,0.05))
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
#risk=as.vector(ifelse(trainScore>6.841980e-02,"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="trainRisk2.txt",sep="\t",quote=F,row.names=F)

#输出test组风险值
rt=read.table(testFile1tcga.LUAD, header=T, sep="\t", row.names=1, check.names=F)
rt$futime=rt$futime/365
testFinalGeneExp=rt[,coxGene]
testScore=predict(multiCox,type="risk",newdata=rt)
outCol=c("futime","fustat",coxGene)
quantile(testScore,seq(0.1,1,0.05))
risk=as.vector(ifelse(testScore>median(testScore),"high","low"))
#risk=as.vector(ifelse(testScore>6.868090,"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="testRisktcga.LUAD.txt",sep="\t",quote=F,row.names=F)

#输出test组风险值
rt=read.table(testFile1GSE41271, header=T, sep="\t", row.names=1, check.names=F)
rt$futime=rt$futime/365
testFinalGeneExp=rt[,coxGene]
testScore=predict(multiCox,type="risk",newdata=rt)
outCol=c("futime","fustat",coxGene)
quantile(testScore,seq(0.1,1,0.05))
#risk=as.vector(ifelse(testScore>median(testScore),"high","low"))
risk=as.vector(ifelse(testScore>26.07118,"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="testRiskGSE41271.txt",sep="\t",quote=F,row.names=F)

#输出test组风险值
rt=read.table(testFile1GSE68465, header=T, sep="\t", row.names=1, check.names=F)
rt$futime=rt$futime/365
testFinalGeneExp=rt[,coxGene]
testScore=predict(multiCox,type="risk",newdata=rt)
outCol=c("futime","fustat",coxGene)
quantile(testScore,seq(0.1,1,0.05))
#risk=as.vector(ifelse(testScore>median(testScore),"high","low"))
risk=as.vector(ifelse(testScore>7.607036,"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="testRiskGSE68465.txt",sep="\t",quote=F,row.names=F)

#输出test组风险值
rt=read.table(testFile1GSE30219, header=T, sep="\t", row.names=1, check.names=F)
rt$futime=rt$futime/365
testFinalGeneExp=rt[,coxGene]
testScore=predict(multiCox,type="risk",newdata=rt)
outCol=c("futime","fustat",coxGene)
quantile(testScore,seq(0.1,1,0.05))
risk=as.vector(ifelse(testScore>median(testScore),"high","low"))
#risk=as.vector(ifelse(testScore>31.39211 ,"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="testRiskGSE30219.txt",sep="\t",quote=F,row.names=F)

#输出test组风险值
rt=read.table(testFile1GSE73403, header=T, sep="\t", row.names=1, check.names=F)
rt$futime=rt$futime/365
testFinalGeneExp=rt[,coxGene]
testScore=predict(multiCox,type="risk",newdata=rt)
outCol=c("futime","fustat",coxGene)
quantile(testScore,seq(0.1,1,0.05))
risk=as.vector(ifelse(testScore>median(testScore),"high","low"))
#risk=as.vector(ifelse(testScore>31.39211 ,"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="testRiskGSE73403.txt",sep="\t",quote=F,row.names=F)

#输出test组风险值
rt=read.table(testFile1GSE37745, header=T, sep="\t", row.names=1, check.names=F)
rt$futime=rt$futime/365
testFinalGeneExp=rt[,coxGene]
testScore=predict(multiCox,type="risk",newdata=rt)
outCol=c("futime","fustat",coxGene)
quantile(testScore,seq(0.1,1,0.05))
risk=as.vector(ifelse(testScore>median(testScore),"high","low"))
#risk=as.vector(ifelse(testScore>31.39211 ,"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="testRiskGSE37745.txt",sep="\t",quote=F,row.names=F)

#输出test组风险值
rt=read.table(testFile1GSE74777, header=T, sep="\t", row.names=1, check.names=F)
rt$futime=rt$futime/365
testFinalGeneExp=rt[,coxGene]
testScore=predict(multiCox,type="risk",newdata=rt)
outCol=c("futime","fustat",coxGene)
quantile(testScore,seq(0.1,1,0.05))
#risk=as.vector(ifelse(testScore>median(testScore),"high","low"))
risk=as.vector(ifelse(testScore>3.867728 ,"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="testRiskGSE74777.txt",sep="\t",quote=F,row.names=F)

#输出test组风险值
rt=read.table(testFile1GSE157010, header=T, sep="\t", row.names=1, check.names=F)
rt$futime=rt$futime/365
testFinalGeneExp=rt[,coxGene]
testScore=predict(multiCox,type="risk",newdata=rt)
outCol=c("futime","fustat",coxGene)
quantile(testScore,seq(0.1,1,0.05))
#risk=as.vector(ifelse(testScore>median(testScore),"high","low"))
risk=as.vector(ifelse(testScore>3.867728 ,"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="testRiskGSE157010.txt",sep="\t",quote=F,row.names=F)

#绘制森林图函数
bioForest=function(coxFile=null,forestFile=null,forestCol=null){
  #读取输入文件
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  #输出图形
  pdf(file=forestFile, width = 6.3,height = 5.5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  #绘制森林图左边的临床信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  #绘制森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol[1], forestCol[2])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)
  axis(1)
  dev.off()
}
#绘制森林图函数
bioForest(coxFile="multiCox2.txt",forestFile="multiCox2.pdf",forestCol=c("Firebrick2","DodgerBlue1"))

####8.4cox画图####

setwd("F:/1/5LUSC/8Model")     #设置工作目录
inputFile="multiCox2.txt"        
outFile="multi.foreast22.pdf"         

rt=read.table(inputFile,header=T,sep="\t",row.names=1,check.names=F)
gene=rownames(rt)
hr=sprintf("%.3f",rt$"HR")
hrLow=sprintf("%.3f",rt$"HR.95L")
hrHigh=sprintf("%.3f",rt$"HR.95H")
Hazard.ratio=paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal=ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

#出图格式
pdf(file=outFile, width = 7, height =7)
n=nrow(rt)
nRow=n+1
ylim=c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))

#森林图左边的基因信息
xlim = c(0,3)
par(mar=c(4,2,1.5,1.5))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

#绘制森林图
par(mar=c(4,1,1.5,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.03,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, "DarkOrchid","Orange2")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()




setwd("F:/1/5LUSC/8Model")    #设置工作目录
inputFile="uni.trainCox.txt"        
outFile="uni.foreast22.pdf"         

rt=read.table(inputFile,header=T,sep="\t",row.names=1,check.names=F)
gene=rownames(rt)
hr=sprintf("%.3f",rt$"HR")
hrLow=sprintf("%.3f",rt$"HR.95L")
hrHigh=sprintf("%.3f",rt$"HR.95H")
Hazard.ratio=paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal=ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

#出图格式
pdf(file=outFile, width = 7, height =6)
n=nrow(rt)
nRow=n+1
ylim=c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))

#森林图左边的基因信息
xlim = c(0,3)
par(mar=c(4,2,1.5,1.5))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

#绘制森林图
par(mar=c(4,1,1.5,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.03,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, "DarkOrchid","Orange2")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()





####8.5生存和ROC####

#引用包
library(survival)
library(survminer)
setwd("F:/1/5LUSC/8Model")     #设置工作目录


svdata1=read.table("risk.all.txt",header=T,sep="\t",row.names = 1)
svdata = svdata1[,c(1,2,length(colnames(svdata1))-1)]
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.15) #默认组内sample不能低于30%
res.cat <- surv_categorize(res.cut)
##统计作图
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl<-list()
group <- res.cat[,"riskScore"] 
survival_dat <- data.frame(group = group)
fit <- survfit(my.surv ~ group)
group <- factor(group, levels = c("low", "high"))
data.survdiff <- survdiff(my.surv ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
if (p.val>0.05) next
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
svsort <- svdata[order(svdata[,"riskScore"]),]
pl[["riskScore"]]<-ggsurvplot(fit, data = survival_dat ,
                              ggtheme = theme_bw(), #想要网格就运行这行
                              conf.int = T, #不画置信区间，想画置信区间就把F改成T
                              #conf.int.style = "step",#置信区间的类型，还可改为ribbon
                              censor = F, #不显示观察值所在的位置
                              palette = c("DarkOrchid","Orange2"), #线的颜色对应高、低
                              legend.title = "All set-riskScore",#基因名写在图例题目的位置
                              font.legend = 11,#图例的字体大小
                              #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                              #在图例上标出高低分界点的表达量，和组内sample数量
                              #在左下角标出pvalue、HR、95% CI
                              #太小的p value标为p < 0.001
                              pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                         paste("p = ",round(p.val,3), sep = "")),
                                           HR, CI, sep = "\n"))
ggsave("Survall.pdf",width = 5,height = 5)






svdata1=read.table("risk.train.txt",header=T,sep="\t",row.names = 1)
svdata = svdata1[,c(1,2,length(colnames(svdata1))-1)]
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.15) #默认组内sample不能低于30%
res.cat <- surv_categorize(res.cut)
##统计作图
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl<-list()
group <- res.cat[,"riskScore"] 
survival_dat <- data.frame(group = group)
fit <- survfit(my.surv ~ group)
group <- factor(group, levels = c("low", "high"))
data.survdiff <- survdiff(my.surv ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
if (p.val>0.05) next
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
svsort <- svdata[order(svdata[,"riskScore"]),]
pl[["riskScore"]]<-ggsurvplot(fit, data = survival_dat ,
                              ggtheme = theme_bw(), #想要网格就运行这行
                              conf.int = T, #不画置信区间，想画置信区间就把F改成T
                              #conf.int.style = "step",#置信区间的类型，还可改为ribbon
                              censor = F, #不显示观察值所在的位置
                              palette = c("DarkOrchid","Orange2"), #线的颜色对应高、低
                              legend.title = "Training set-riskScore",#基因名写在图例题目的位置
                              font.legend = 11,#图例的字体大小
                              #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                              #在图例上标出高低分界点的表达量，和组内sample数量
                              #在左下角标出pvalue、HR、95% CI
                              #太小的p value标为p < 0.001
                              pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                         paste("p = ",round(p.val,3), sep = "")),
                                           HR, CI, sep = "\n"))
ggsave("Survtrain.pdf",width = 5,height = 5)




svdata1=read.table("risk.test.txt",header=T,sep="\t",row.names = 1)
svdata = svdata1[,c(1,2,length(colnames(svdata1))-1)]
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.15) #默认组内sample不能低于30%
res.cat <- surv_categorize(res.cut)
##统计作图
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl<-list()
group <- res.cat[,"riskScore"] 
survival_dat <- data.frame(group = group)
fit <- survfit(my.surv ~ group)
group <- factor(group, levels = c("low", "high"))
data.survdiff <- survdiff(my.surv ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
if (p.val>0.05) next
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
svsort <- svdata[order(svdata[,"riskScore"]),]
pl[["riskScore"]]<-ggsurvplot(fit, data = survival_dat ,
                              ggtheme = theme_bw(), #想要网格就运行这行
                              conf.int = T, #不画置信区间，想画置信区间就把F改成T
                              #conf.int.style = "step",#置信区间的类型，还可改为ribbon
                              censor = F, #不显示观察值所在的位置
                              palette = c("DarkOrchid","Orange2"), #线的颜色对应高、低
                              legend.title = "Testing set-riskScore",#基因名写在图例题目的位置
                              font.legend = 11,#图例的字体大小
                              #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                              #在图例上标出高低分界点的表达量，和组内sample数量
                              #在左下角标出pvalue、HR、95% CI
                              #太小的p value标为p < 0.001
                              pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                         paste("p = ",round(p.val,3), sep = "")),
                                           HR, CI, sep = "\n"))
ggsave("Survtest.pdf",width = 5,height = 5)




svdata1=read.table("testRisktcga.LUAD.txt",header=T,sep="\t",row.names = 1)
svdata = svdata1[,c(1,2,length(colnames(svdata1))-1)]
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.15) #默认组内sample不能低于30%
res.cat <- surv_categorize(res.cut)
##统计作图
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl<-list()
group <- res.cat[,"riskScore"] 
survival_dat <- data.frame(group = group)
fit <- survfit(my.surv ~ group)
group <- factor(group, levels = c("low", "high"))
data.survdiff <- survdiff(my.surv ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
if (p.val>0.05) next
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
svsort <- svdata[order(svdata[,"riskScore"]),]
pl[["riskScore"]]<-ggsurvplot(fit, data = survival_dat ,
                              ggtheme = theme_bw(), #想要网格就运行这行
                              conf.int = T, #不画置信区间，想画置信区间就把F改成T
                              #conf.int.style = "step",#置信区间的类型，还可改为ribbon
                              censor = F, #不显示观察值所在的位置
                              palette = c("DarkOrchid","Orange2"), #线的颜色对应高、低
                              legend.title = "TCGA-LUAD-riskScore",#基因名写在图例题目的位置
                              font.legend = 11,#图例的字体大小
                              #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                              #在图例上标出高低分界点的表达量，和组内sample数量
                              #在左下角标出pvalue、HR、95% CI
                              #太小的p value标为p < 0.001
                              pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                         paste("p = ",round(p.val,3), sep = "")),
                                           HR, CI, sep = "\n"))
ggsave("SurvTCGA-LUAD.pdf",width = 5,height = 5)




svdata1=read.table("testRiskGSE41271.txt",header=T,sep="\t",row.names = 1)
svdata = svdata1[,c(1,2,length(colnames(svdata1))-1)]
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.15) #默认组内sample不能低于30%
res.cat <- surv_categorize(res.cut)
##统计作图
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl<-list()
group <- res.cat[,"riskScore"] 
survival_dat <- data.frame(group = group)
fit <- survfit(my.surv ~ group)
group <- factor(group, levels = c("low", "high"))
data.survdiff <- survdiff(my.surv ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
if (p.val>0.05) next
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
svsort <- svdata[order(svdata[,"riskScore"]),]
pl[["riskScore"]]<-ggsurvplot(fit, data = survival_dat ,
                              ggtheme = theme_bw(), #想要网格就运行这行
                              conf.int = T, #不画置信区间，想画置信区间就把F改成T
                              #conf.int.style = "step",#置信区间的类型，还可改为ribbon
                              censor = F, #不显示观察值所在的位置
                              palette = c("DarkOrchid","Orange2"), #线的颜色对应高、低
                              legend.title = "GSE41271-riskScore",#基因名写在图例题目的位置
                              font.legend = 11,#图例的字体大小
                              #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                              #在图例上标出高低分界点的表达量，和组内sample数量
                              #在左下角标出pvalue、HR、95% CI
                              #太小的p value标为p < 0.001
                              pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                         paste("p = ",round(p.val,3), sep = "")),
                                           HR, CI, sep = "\n"))
ggsave("SurvGSE41271.pdf",width = 5,height = 5)



svdata1=read.table("testRiskGSE68465.txt",header=T,sep="\t",row.names = 1)
svdata = svdata1[,c(1,2,length(colnames(svdata1))-1)]
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.15) #默认组内sample不能低于30%
res.cat <- surv_categorize(res.cut)
##统计作图
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl<-list()
group <- res.cat[,"riskScore"] 
survival_dat <- data.frame(group = group)
fit <- survfit(my.surv ~ group)
group <- factor(group, levels = c("low", "high"))
data.survdiff <- survdiff(my.surv ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
if (p.val>0.05) next
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
svsort <- svdata[order(svdata[,"riskScore"]),]
pl[["riskScore"]]<-ggsurvplot(fit, data = survival_dat ,
                              ggtheme = theme_bw(), #想要网格就运行这行
                              conf.int = T, #不画置信区间，想画置信区间就把F改成T
                              #conf.int.style = "step",#置信区间的类型，还可改为ribbon
                              censor = F, #不显示观察值所在的位置
                              palette = c("DarkOrchid","Orange2"), #线的颜色对应高、低
                              legend.title = "GSE68465-riskScore",#基因名写在图例题目的位置
                              font.legend = 11,#图例的字体大小
                              #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                              #在图例上标出高低分界点的表达量，和组内sample数量
                              #在左下角标出pvalue、HR、95% CI
                              #太小的p value标为p < 0.001
                              pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                         paste("p = ",round(p.val,3), sep = "")),
                                           HR, CI, sep = "\n"))
ggsave("SurvGSE68465.pdf",width = 5,height = 5)





svdata1=read.table("testRiskGSE30219.txt",header=T,sep="\t",row.names = 1)
svdata = svdata1[,c(1,2,length(colnames(svdata1))-1)]
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.15) #默认组内sample不能低于30%
res.cat <- surv_categorize(res.cut)
##统计作图
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl<-list()
group <- res.cat[,"riskScore"] 
survival_dat <- data.frame(group = group)
fit <- survfit(my.surv ~ group)
group <- factor(group, levels = c("low", "high"))
data.survdiff <- survdiff(my.surv ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
if (p.val>0.05) next
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
svsort <- svdata[order(svdata[,"riskScore"]),]
pl[["riskScore"]]<-ggsurvplot(fit, data = survival_dat ,
                              ggtheme = theme_bw(), #想要网格就运行这行
                              conf.int = T, #不画置信区间，想画置信区间就把F改成T
                              #conf.int.style = "step",#置信区间的类型，还可改为ribbon
                              censor = F, #不显示观察值所在的位置
                              palette = c("DarkOrchid","Orange2"), #线的颜色对应高、低
                              legend.title = "GSE30219-riskScore",#基因名写在图例题目的位置
                              font.legend = 11,#图例的字体大小
                              #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                              #在图例上标出高低分界点的表达量，和组内sample数量
                              #在左下角标出pvalue、HR、95% CI
                              #太小的p value标为p < 0.001
                              pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                         paste("p = ",round(p.val,3), sep = "")),
                                           HR, CI, sep = "\n"))
ggsave("SurvGSE30219.pdf",width = 5,height = 5)




svdata1=read.table("testRiskGSE73403.txt",header=T,sep="\t",row.names = 1)
svdata = svdata1[,c(1,2,length(colnames(svdata1))-1)]
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.15) #默认组内sample不能低于30%
res.cat <- surv_categorize(res.cut)
##统计作图
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl<-list()
group <- res.cat[,"riskScore"] 
survival_dat <- data.frame(group = group)
fit <- survfit(my.surv ~ group)
group <- factor(group, levels = c("low", "high"))
data.survdiff <- survdiff(my.surv ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
if (p.val>0.05) next
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
svsort <- svdata[order(svdata[,"riskScore"]),]
pl[["riskScore"]]<-ggsurvplot(fit, data = survival_dat ,
                              ggtheme = theme_bw(), #想要网格就运行这行
                              conf.int = T, #不画置信区间，想画置信区间就把F改成T
                              #conf.int.style = "step",#置信区间的类型，还可改为ribbon
                              censor = F, #不显示观察值所在的位置
                              palette = c("DarkOrchid","Orange2"), #线的颜色对应高、低
                              legend.title = "GSE73403-riskScore",#基因名写在图例题目的位置
                              font.legend = 11,#图例的字体大小
                              #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                              #在图例上标出高低分界点的表达量，和组内sample数量
                              #在左下角标出pvalue、HR、95% CI
                              #太小的p value标为p < 0.001
                              pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                         paste("p = ",round(p.val,3), sep = "")),
                                           HR, CI, sep = "\n"))
ggsave("SurvGSE73403.pdf",width = 5,height = 5)




svdata1=read.table("testRiskGSE37745.txt",header=T,sep="\t",row.names = 1)
svdata = svdata1[,c(1,2,length(colnames(svdata1))-1)]
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.15) #默认组内sample不能低于30%
res.cat <- surv_categorize(res.cut)
##统计作图
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl<-list()
group <- res.cat[,"riskScore"] 
survival_dat <- data.frame(group = group)
fit <- survfit(my.surv ~ group)
group <- factor(group, levels = c("low", "high"))
data.survdiff <- survdiff(my.surv ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
if (p.val>0.05) next
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
svsort <- svdata[order(svdata[,"riskScore"]),]
pl[["riskScore"]]<-ggsurvplot(fit, data = survival_dat ,
                              ggtheme = theme_bw(), #想要网格就运行这行
                              conf.int = T, #不画置信区间，想画置信区间就把F改成T
                              #conf.int.style = "step",#置信区间的类型，还可改为ribbon
                              censor = F, #不显示观察值所在的位置
                              palette = c("DarkOrchid","Orange2"), #线的颜色对应高、低
                              legend.title = "GSE37745-riskScore",#基因名写在图例题目的位置
                              font.legend = 11,#图例的字体大小
                              #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                              #在图例上标出高低分界点的表达量，和组内sample数量
                              #在左下角标出pvalue、HR、95% CI
                              #太小的p value标为p < 0.001
                              pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                         paste("p = ",round(p.val,3), sep = "")),
                                           HR, CI, sep = "\n"))
ggsave("SurvGSE37745.pdf",width = 5,height = 5)





svdata1=read.table("testRiskGSE74777.txt",header=T,sep="\t",row.names = 1)
svdata = svdata1[,c(1,2,length(colnames(svdata1))-1)]
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.15) #默认组内sample不能低于30%
res.cat <- surv_categorize(res.cut)
##统计作图
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl<-list()
group <- res.cat[,"riskScore"] 
survival_dat <- data.frame(group = group)
fit <- survfit(my.surv ~ group)
group <- factor(group, levels = c("low", "high"))
data.survdiff <- survdiff(my.surv ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
if (p.val>0.05) next
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
svsort <- svdata[order(svdata[,"riskScore"]),]
pl[["riskScore"]]<-ggsurvplot(fit, data = survival_dat ,
                              ggtheme = theme_bw(), #想要网格就运行这行
                              conf.int = T, #不画置信区间，想画置信区间就把F改成T
                              #conf.int.style = "step",#置信区间的类型，还可改为ribbon
                              censor = F, #不显示观察值所在的位置
                              palette = c("DarkOrchid","Orange2"), #线的颜色对应高、低
                              legend.title = "GSE74777-riskScore",#基因名写在图例题目的位置
                              font.legend = 11,#图例的字体大小
                              #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                              #在图例上标出高低分界点的表达量，和组内sample数量
                              #在左下角标出pvalue、HR、95% CI
                              #太小的p value标为p < 0.001
                              pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                         paste("p = ",round(p.val,3), sep = "")),
                                           HR, CI, sep = "\n"))
ggsave("SurvGSE74777.pdf",width = 5,height = 5)



svdata1=read.table("testRiskGSE157010.txt",header=T,sep="\t",row.names = 1)
svdata = svdata1[,c(1,2,length(colnames(svdata1))-1)]
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.15) #默认组内sample不能低于30%
res.cat <- surv_categorize(res.cut)
##统计作图
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl<-list()
group <- res.cat[,"riskScore"] 
survival_dat <- data.frame(group = group)
fit <- survfit(my.surv ~ group)
group <- factor(group, levels = c("low", "high"))
data.survdiff <- survdiff(my.surv ~ group)
p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
if (p.val>0.05) next
HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
svsort <- svdata[order(svdata[,"riskScore"]),]
pl[["riskScore"]]<-ggsurvplot(fit, data = survival_dat ,
                              ggtheme = theme_bw(), #想要网格就运行这行
                              conf.int = T, #不画置信区间，想画置信区间就把F改成T
                              #conf.int.style = "step",#置信区间的类型，还可改为ribbon
                              censor = F, #不显示观察值所在的位置
                              palette = c("DarkOrchid","Orange2"), #线的颜色对应高、低
                              legend.title = "GSE157010-riskScore",#基因名写在图例题目的位置
                              font.legend = 11,#图例的字体大小
                              #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                              #在图例上标出高低分界点的表达量，和组内sample数量
                              #在左下角标出pvalue、HR、95% CI
                              #太小的p value标为p < 0.001
                              pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                         paste("p = ",round(p.val,3), sep = "")),
                                           HR, CI, sep = "\n"))
ggsave("SurvGSE157010.pdf",width = 5,height = 5)






#引用包
library(survival)
library(survminer)
library(timeROC)

riskFile="risk.all.txt"     #风险文件
cliFile="clinical1.txt"      #临床数据文件

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

#定义颜色
bioCol=c("DarkOrchid","Orange2","MediumSeaGreen","NavyBlue","#8B668B","#FF4500","#135612","#561214")

######绘制1 3 5年的ROC曲线
ROC_rt=timeROC(T=risk$futime,delta=risk$fustat,
               marker=risk$riskScore,cause=1,
               weighting='aalen',
               times=c(1,3,5),ROC=TRUE)
pdf(file="ROC.all.pdf", width=5, height=5)
plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=4)
plot(ROC_rt,time=3,col=bioCol[2],add=TRUE,title=FALSE,lwd=4)
plot(ROC_rt,time=5,col=bioCol[3],add=TRUE,title=FALSE,lwd=4)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=bioCol[1:3], lwd=4, bty = 'n',title = "All set")
dev.off()


######绘制临床的ROC曲线
predictTime=5     #定义预测年限
aucText=c()
pdf(file="cliROC.all.pdf", width=5.5, height=5.5)
#绘制风险得分的ROC曲线
i=3
ROC_rt=timeROC(T=risk$futime,
               delta=risk$fustat,
               marker=risk$riskScore, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=4)
aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
abline(0,1)
#对临床数据进行循环，绘制临床数据的ROC曲线
for(i in 4:ncol(rt)){
  ROC_rt=timeROC(T=rt$futime,
                 delta=rt$fustat,
                 marker=rt[,i], cause=1,
                 weighting='aalen',
                 times=c(predictTime),ROC=TRUE)
  plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=4, add=TRUE)
  aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
}
#绘制图例，得到ROC曲线下的面积
legend("bottomright", aucText,lwd=4,bty="n",col=bioCol[1:(ncol(rt)-1)],title = "All set")
dev.off()





#引用包
library(survival)
library(survminer)
library(timeROC)

riskFile="risk.train.txt"     #风险文件
cliFile="clinical1.txt"      #临床数据文件

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

#定义颜色
bioCol=c("DarkOrchid","Orange2","MediumSeaGreen","NavyBlue","#8B668B","#FF4500","#135612","#561214")

######绘制1 3 5年的ROC曲线
ROC_rt=timeROC(T=risk$futime,delta=risk$fustat,
               marker=risk$riskScore,cause=1,
               weighting='aalen',
               times=c(1,3,5),ROC=TRUE)
pdf(file="ROC.train.pdf", width=5, height=5)
plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=4)
plot(ROC_rt,time=3,col=bioCol[2],add=TRUE,title=FALSE,lwd=4)
plot(ROC_rt,time=5,col=bioCol[3],add=TRUE,title=FALSE,lwd=4)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=bioCol[1:3], lwd=4, bty = 'n',title = "Training set")
dev.off()


######绘制临床的ROC曲线
predictTime=5     #定义预测年限
aucText=c()
pdf(file="cliROC.train.pdf", width=5.5, height=5.5)
#绘制风险得分的ROC曲线
i=3
ROC_rt=timeROC(T=risk$futime,
               delta=risk$fustat,
               marker=risk$riskScore, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=4)
aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
abline(0,1)
#对临床数据进行循环，绘制临床数据的ROC曲线
for(i in 4:ncol(rt)){
  ROC_rt=timeROC(T=rt$futime,
                 delta=rt$fustat,
                 marker=rt[,i], cause=1,
                 weighting='aalen',
                 times=c(predictTime),ROC=TRUE)
  plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=4, add=TRUE)
  aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
}
#绘制图例，得到ROC曲线下的面积
legend("bottomright", aucText,lwd=4,bty="n",col=bioCol[1:(ncol(rt)-1)],title = "Training set")
dev.off()





#引用包
library(survival)
library(survminer)
library(timeROC)

riskFile="risk.test.txt"     #风险文件
cliFile="clinical1.txt"      #临床数据文件

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

#定义颜色
bioCol=c("DarkOrchid","Orange2","MediumSeaGreen","NavyBlue","#8B668B","#FF4500","#135612","#561214")

######绘制1 3 5年的ROC曲线
ROC_rt=timeROC(T=risk$futime,delta=risk$fustat,
               marker=risk$riskScore,cause=1,
               weighting='aalen',
               times=c(1,3,5),ROC=TRUE)
pdf(file="ROC.test.pdf", width=5, height=5)
plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=4)
plot(ROC_rt,time=3,col=bioCol[2],add=TRUE,title=FALSE,lwd=4)
plot(ROC_rt,time=5,col=bioCol[3],add=TRUE,title=FALSE,lwd=4)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=bioCol[1:3], lwd=4, bty = 'n',title = "Testing set")
dev.off()


######绘制临床的ROC曲线
predictTime=5     #定义预测年限
aucText=c()
pdf(file="cliROC.test.pdf", width=5.5, height=5.5)
#绘制风险得分的ROC曲线
i=3
ROC_rt=timeROC(T=risk$futime,
               delta=risk$fustat,
               marker=risk$riskScore, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=4)
aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
abline(0,1)
#对临床数据进行循环，绘制临床数据的ROC曲线
for(i in 4:ncol(rt)){
  ROC_rt=timeROC(T=rt$futime,
                 delta=rt$fustat,
                 marker=rt[,i], cause=1,
                 weighting='aalen',
                 times=c(predictTime),ROC=TRUE)
  plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=4, add=TRUE)
  aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
}
#绘制图例，得到ROC曲线下的面积
legend("bottomright", aucText,lwd=4,bty="n",col=bioCol[1:(ncol(rt)-1)],title = "Testing set")
dev.off()






####9.1indep####

library(survival)      #引用包
setwd("F:/1/5LUSC/9Application")      #设置工作目录

############绘制森林图函数
bioForest=function(coxFile=null, forestFile=null, forestCol=null){
  #读取输入文件
  rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #输出图形
  pdf(file=forestFile, width=6.6, height=4)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #绘制森林图左边的临床信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  #绘制右边的森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
  dev.off()
}
############绘制森林图函数

############独立预后分析函数
indep=function(riskFile=null, cliFile=null, project=null){
  risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #读取风险文件
  cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)      #读取临床文件
  
  #数据合并
  sameSample=intersect(row.names(cli),row.names(risk))
  risk=risk[sameSample,]
  cli=cli[sameSample,]
  rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])
  
  #单因素独立预后分析
  uniCoxFile=paste0(project,".uniCox.txt")
  uniCoxPdf=paste0(project,".uniCox.pdf")
  uniTab=data.frame()
  for(i in colnames(rt[,3:ncol(rt)])){
    cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
    coxSummary = summary(cox)
    uniTab=rbind(uniTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
  write.table(uniTab,file=uniCoxFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=uniCoxFile, forestFile=uniCoxPdf, forestCol="Orange2")
  
  
  #多因素独立预后分析
  multiCoxFile=paste0(project,".multiCox.txt")
  multiCoxPdf=paste0(project,".multiCox.pdf")
  uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<1,]
  rt1=rt[,c("futime","fustat",as.vector(uniTab[,"id"]))]
  multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
  multiCoxSum=summary(multiCox)
  multiTab=data.frame()
  multiTab=cbind(
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  multiTab=cbind(id=row.names(multiTab),multiTab)
  write.table(multiTab, file=multiCoxFile, sep="\t", row.names=F, quote=F)
  bioForest(coxFile=multiCoxFile, forestFile=multiCoxPdf, forestCol="DarkOrchid")
}
############独立预后分析函数

#独立预后分析
indep(riskFile="risk.all.txt", cliFile="clinical1.txt", project="all")

####9.2C-index####

#引用包
library(dplyr)
library(survival)
library(rms)
library(pec)

riskFile="risk.all.txt"     #风险文件
cliFile="clinical2.txt"      #临床数据文件
setwd("F:/1/5LUSC/9Application")    #修改工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

#定义颜色
bioCol=c("DarkOrchid","Orange2","MediumSeaGreen","NavyBlue","Firebrick3","#8B668B","#FF4500")

#C-index值计算
riskScore=cph(Surv(futime,fustat)~riskScore, data=rt, surv=TRUE)
Age=cph(Surv(futime,fustat)~Age, data=rt, surv=TRUE)
Gender=cph(Surv(futime,fustat)~Gender, data=rt, surv=TRUE)
Stage=cph(Surv(futime,fustat)~Stage, data=rt, surv=TRUE)
Subdivision=cph(Surv(futime,fustat)~Subdivision, data=rt, surv=TRUE)
c_index  <- cindex(list("Risk score"=riskScore, 
                        "Age"=Age,
                        "Gender"=Gender,
                        "Stage"=Stage,
                        "Subdivision"=Subdivision),
                   formula=Surv(futime,fustat)~ .,
                   data=rt,
                   eval.times=seq(0,15,1),
                   splitMethod="bootcv",
                   B=1000
)
#输出图形
pdf(file="C-index.pdf", width=6, height=6)
plot(c_index, xlim=c(0,15), ylim=c(0.4,0.8), col=bioCol, legend.x=6, legend.y=0.82, legend.cex=1)
dev.off()

####9.3cliGroupSur####

#引用包
library(survival)
library(survminer)

riskFile="risk.all.txt"      #风险文件
cliFile="clinicalage.txt"       #临床数据文件
setwd("F:/1/5LUSC/9Application")      #设置工作目录

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取临床文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cliName=colnames(cli)[1]

#数据合并
sameSample=intersect(row.names(cli), row.names(risk))
risk=risk[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(futime=risk[,1], fustat=risk[,2], cli, risk[,"risk",drop=F])
colnames(rt)=c("futime", "fustat", "clinical", "Risk")
tab=table(rt[,"clinical"])
tab=tab[tab!=0]

#对每个临床信息里面的每个分组进行循环
for(j in names(tab)){
  rt1=rt[(rt[,"clinical"]==j),]
  tab1=table(rt1[,"Risk"])
  tab1=tab1[tab1!=0]
  labels=names(tab1)
  if(length(labels)!=2){next}
  if((cliName=="age") | (cliName=="Age") | (cliName=="AGE")){
    titleName=paste0("age",j)
  }
  
  #计算高低风险组差异pvalue
  diff=survdiff(Surv(futime, fustat) ~Risk,data = rt1)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  
  #绘制生存曲线
  fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt1)
  surPlot=ggsurvplot(fit, 
                     data=rt1,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     title=paste0("Patients with ",j),
                     legend.title="Risk",
                     legend.labs=labels,
                     font.legend=12,
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette=c("DarkOrchid","Orange2"),
                     risk.table=F,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  
  #输出图片
  j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
  pdf(file=paste0("survival.",cliName,"_",j,".pdf"), onefile = FALSE,
      width = 6,#图片的宽度
      height =5)#图片的高度
  print(surPlot)
  dev.off()
}

#引用包
library(survival)
library(survminer)

riskFile="risk.all.txt"      #风险文件
cliFile="clinicalgender.txt"       #临床数据文件

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取临床文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cliName=colnames(cli)[1]

#数据合并
sameSample=intersect(row.names(cli), row.names(risk))
risk=risk[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(futime=risk[,1], fustat=risk[,2], cli, risk[,"risk",drop=F])
colnames(rt)=c("futime", "fustat", "clinical", "Risk")
tab=table(rt[,"clinical"])
tab=tab[tab!=0]

#对每个临床信息里面的每个分组进行循环
for(j in names(tab)){
  rt1=rt[(rt[,"clinical"]==j),]
  tab1=table(rt1[,"Risk"])
  tab1=tab1[tab1!=0]
  labels=names(tab1)
  if(length(labels)!=2){next}
  if((cliName=="age") | (cliName=="Age") | (cliName=="AGE")){
    titleName=paste0("age",j)
  }
  
  #计算高低风险组差异pvalue
  diff=survdiff(Surv(futime, fustat) ~Risk,data = rt1)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  
  #绘制生存曲线
  fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt1)
  surPlot=ggsurvplot(fit, 
                     data=rt1,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     title=paste0("Patients with ",j),
                     legend.title="Risk",
                     legend.labs=labels,
                     font.legend=12,
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette=c("DarkOrchid","Orange2"),
                     risk.table=F,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  
  #输出图片
  j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
  pdf(file=paste0("survival.",cliName,"_",j,".pdf"), onefile = FALSE,
      width = 6,#图片的宽度
      height =5)#图片的高度
  print(surPlot)
  dev.off()
}


#引用包
library(survival)
library(survminer)

riskFile="risk.all.txt"      #风险文件
cliFile="clinicalT.txt"       #临床数据文件

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取临床文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cliName=colnames(cli)[1]

#数据合并
sameSample=intersect(row.names(cli), row.names(risk))
risk=risk[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(futime=risk[,1], fustat=risk[,2], cli, risk[,"risk",drop=F])
colnames(rt)=c("futime", "fustat", "clinical", "Risk")
tab=table(rt[,"clinical"])
tab=tab[tab!=0]

#对每个临床信息里面的每个分组进行循环
for(j in names(tab)){
  rt1=rt[(rt[,"clinical"]==j),]
  tab1=table(rt1[,"Risk"])
  tab1=tab1[tab1!=0]
  labels=names(tab1)
  if(length(labels)!=2){next}
  if((cliName=="age") | (cliName=="Age") | (cliName=="AGE")){
    titleName=paste0("age",j)
  }
  
  #计算高低风险组差异pvalue
  diff=survdiff(Surv(futime, fustat) ~Risk,data = rt1)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  
  #绘制生存曲线
  fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt1)
  surPlot=ggsurvplot(fit, 
                     data=rt1,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     title=paste0("Patients with ",j),
                     legend.title="Risk",
                     legend.labs=labels,
                     font.legend=12,
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette=c("DarkOrchid","Orange2"),
                     risk.table=F,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  
  #输出图片
  j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
  pdf(file=paste0("survival.",cliName,"_",j,".pdf"), onefile = FALSE,
      width = 6,#图片的宽度
      height =5)#图片的高度
  print(surPlot)
  dev.off()
}


#引用包
library(survival)
library(survminer)

riskFile="risk.all.txt"      #风险文件
cliFile="clinicalN.txt"       #临床数据文件

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取临床文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cliName=colnames(cli)[1]

#数据合并
sameSample=intersect(row.names(cli), row.names(risk))
risk=risk[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(futime=risk[,1], fustat=risk[,2], cli, risk[,"risk",drop=F])
colnames(rt)=c("futime", "fustat", "clinical", "Risk")
tab=table(rt[,"clinical"])
tab=tab[tab!=0]

#对每个临床信息里面的每个分组进行循环
for(j in names(tab)){
  rt1=rt[(rt[,"clinical"]==j),]
  tab1=table(rt1[,"Risk"])
  tab1=tab1[tab1!=0]
  labels=names(tab1)
  if(length(labels)!=2){next}
  if((cliName=="age") | (cliName=="Age") | (cliName=="AGE")){
    titleName=paste0("age",j)
  }
  
  #计算高低风险组差异pvalue
  diff=survdiff(Surv(futime, fustat) ~Risk,data = rt1)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  
  #绘制生存曲线
  fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt1)
  surPlot=ggsurvplot(fit, 
                     data=rt1,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     title=paste0("Patients with ",j),
                     legend.title="Risk",
                     legend.labs=labels,
                     font.legend=12,
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette=c("DarkOrchid","Orange2"),
                     risk.table=F,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  
  #输出图片
  j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
  pdf(file=paste0("survival.",cliName,"_",j,".pdf"), onefile = FALSE,
      width = 6,#图片的宽度
      height =5)#图片的高度
  print(surPlot)
  dev.off()
}


#引用包
library(survival)
library(survminer)

riskFile="risk.all.txt"      #风险文件
cliFile="clinicalstage.txt"       #临床数据文件

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取临床文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cliName=colnames(cli)[1]

#数据合并
sameSample=intersect(row.names(cli), row.names(risk))
risk=risk[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(futime=risk[,1], fustat=risk[,2], cli, risk[,"risk",drop=F])
colnames(rt)=c("futime", "fustat", "clinical", "Risk")
tab=table(rt[,"clinical"])
tab=tab[tab!=0]

#对每个临床信息里面的每个分组进行循环
for(j in names(tab)){
  rt1=rt[(rt[,"clinical"]==j),]
  tab1=table(rt1[,"Risk"])
  tab1=tab1[tab1!=0]
  labels=names(tab1)
  if(length(labels)!=2){next}
  if((cliName=="age") | (cliName=="Age") | (cliName=="AGE")){
    titleName=paste0("age",j)
  }
  
  #计算高低风险组差异pvalue
  diff=survdiff(Surv(futime, fustat) ~Risk,data = rt1)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  
  #绘制生存曲线
  fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt1)
  surPlot=ggsurvplot(fit, 
                     data=rt1,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     title=paste0("Patients with ",j),
                     legend.title="Risk",
                     legend.labs=labels,
                     font.legend=12,
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette=c("DarkOrchid","Orange2"),
                     risk.table=F,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  
  #输出图片
  j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
  pdf(file=paste0("survival.",cliName,"_",j,".pdf"), onefile = FALSE,
      width = 6,#图片的宽度
      height =5)#图片的高度
  print(surPlot)
  dev.off()
}




#引用包
library(survival)
library(survminer)

riskFile="risk.all.txt"      #风险文件
cliFile="clinicalSubdivision.txt"       #临床数据文件

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取临床文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cliName=colnames(cli)[1]

#数据合并
sameSample=intersect(row.names(cli), row.names(risk))
risk=risk[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(futime=risk[,1], fustat=risk[,2], cli, risk[,"risk",drop=F])
colnames(rt)=c("futime", "fustat", "clinical", "Risk")
tab=table(rt[,"clinical"])
tab=tab[tab!=0]

#对每个临床信息里面的每个分组进行循环
for(j in names(tab)){
  rt1=rt[(rt[,"clinical"]==j),]
  tab1=table(rt1[,"Risk"])
  tab1=tab1[tab1!=0]
  labels=names(tab1)
  if(length(labels)!=2){next}
  if((cliName=="age") | (cliName=="Age") | (cliName=="AGE")){
    titleName=paste0("age",j)
  }
  
  #计算高低风险组差异pvalue
  diff=survdiff(Surv(futime, fustat) ~Risk,data = rt1)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  
  #绘制生存曲线
  fit <- survfit(Surv(futime, fustat) ~ Risk, data = rt1)
  surPlot=ggsurvplot(fit, 
                     data=rt1,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     title=paste0("Patients with ",j),
                     legend.title="Risk",
                     legend.labs=labels,
                     font.legend=12,
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette=c("DarkOrchid","Orange2"),
                     risk.table=F,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  
  #输出图片
  j=gsub(">=","ge",j);j=gsub("<=","le",j);j=gsub(">","gt",j);j=gsub("<","lt",j)
  pdf(file=paste0("survival.",cliName,"_",j,".pdf"), onefile = FALSE,
      width = 6,#图片的宽度
      height =5)#图片的高度
  print(surPlot)
  dev.off()
}


####9.4Nomo####


#引用包
library(survival)
library(regplot)
library(rms)

riskFile="risk.all.txt"      #风险输入文件
cliFile="clinical4.txt"       #临床数据文件
setwd("E:/6LUAD-CAF/1data/9application")   #修改工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]
cli$Age=as.numeric(cli$Age)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1[,c("futime", "fustat", "risk")], cli)

#绘制列线图
res.cox=coxph(Surv(futime, fustat) ~ . , data = rt)
nom1=regplot(res.cox,
             plots = c("density", "boxes"),
             clickable=F,
             title="",
             points=TRUE,
             droplines=TRUE,
             observation=rt[50,],
             rank="sd",
             failtime = c(1,3,5),
             prfail = F)

#列线图风险打分
nomoRisk=predict(res.cox, data=rt, type="risk")
rt=cbind(risk1, Nomogram=nomoRisk)
outTab=rbind(ID=colnames(rt), rt)
write.table(outTab, file="nomoRisk.txt", sep="\t", col.names=F, quote=F)

#校准曲线
pdf(file="calibration.pdf", width=6, height=6)
#1年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1)
cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=3, col="Firebrick2", sub=F)
#3年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=3)
cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=3, col="MediumSeaGreen", sub=F, add=T)
#5年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=5)
cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=3, col="NavyBlue", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("Firebrick3","MediumSeaGreen","NavyBlue"), lwd=3, bty = 'n')
dev.off()


####9.5riskDiff####


#引用包
library(limma)
expFile="TCGA_LUSC.txt"  #表达数据文件
riskFile="risk.all.txt"       #风险文件
logFCfilter=1 #logFC过滤条件
fdrFilter=0.05#fdr过滤条件
setwd("F:/1/5LUSC/9Application")  #设置工作目录

#读取表达数据文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

##去除正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)
data=t(data)

#读取risk文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(colnames(data), row.names(risk))
data=data[,sameSample]
risk=risk[sameSample,]

#提取low risk和high risk样品
riskLow=risk[risk$risk=="low",]
riskHigh=risk[risk$risk=="high",]
dataLow=data[,row.names(riskLow)]
dataHigh=data[,row.names(riskHigh)]
data=cbind(dataLow,dataHigh)
data=data[rowMeans(data)>1,]
conNum=ncol(dataLow)
treatNum=ncol(dataHigh)
Type=c(rep(1,conNum), rep(2,treatNum))

#差异分析
outTab=data.frame()
for(i in row.names(data)){
  rt=data.frame(expression=data[i,], Type=Type)
  wilcoxTest=wilcox.test(expression ~ Type, data=rt)
  pvalue=wilcoxTest$p.value
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(gene=i,lowMean=conGeneMeans,highMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
outTab=cbind(outTab, fdr=fdr)

#输出差异表格
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff, file="riskDiff.txt", sep="\t", row.names=F, quote=F)



####9.6GO####


#引用包
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)

pvalueFilter=0.05       #p值过滤条件
qvalueFilter=1       #矫正后的p值过滤条件

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}
ontology.col=c("#00AFBB", "#E7B800", "#90EE90")

setwd("F:/1/5LUSC/9Application")      #设置工作目录
rt=read.table("riskDiff.txt", header=T, sep="\t", check.names=F)     #读取输入文件

#基因名字转换为基因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]#去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO富集分析
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
#保存富集结果
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

#定义显示GO的数目
showNum=10
if(nrow(GO)<30){
  showNum=nrow(GO)
}

#柱状图
pdf(file="GObarplot.pdf", width=10, height=7)
bar=barplot(kk, drop=TRUE, showCategory=showNum, label_format=30, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()

#气泡图
pdf(file="GObubble.pdf", width=10, height=7)
bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=30, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

###########绘制GO圈图
data=GO[order(GO$pvalue),]
datasig=data[data$pvalue<0.05,,drop=F]
BP = datasig[datasig$ONTOLOGY=="BP",,drop=F]
CC = datasig[datasig$ONTOLOGY=="CC",,drop=F]
MF = datasig[datasig$ONTOLOGY=="MF",,drop=F]
BP = head(BP,6)
CC = head(CC,6)
MF = head(MF,6)
data = rbind(BP,CC,MF)
main.col = ontology.col[as.numeric(as.factor(data$ONTOLOGY))]

#整理圈图数据
BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(GO=data$ID,start=1,end=max(BgGene))
rownames(df) = df$GO
bed2 = data.frame(GO=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(GO=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(GO=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5

#绘制圈图主体部分
pdf("GO.circlize.pdf",width=10,height=10)
par(omi=c(0.1,0.1,0.1,1.5))
circos.par(track.min=c(0.01,0.01))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                      #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                    })
circos.clear()
#绘制圈图中间的图例
middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',ontology.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)
#绘制GO分类的图例
main.legend = Legend(
  labels = c("Biological Process", "Molecular Function","Cellular Component"),  type="points",pch=15,
  legend_gp = gpar(col=ontology.col), title_position = "topcenter",
  title = "ONTOLOGY", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)
#绘制pvalue的图例
logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()


####9.7KEGG####


#引用包
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)

pvalueFilter=0.05      #p值过滤条件
qvalueFilter=1      #矫正后的p值过滤条件

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}

setwd("F:/1/5LUSC/9Application")  #设置工作目录
rt=read.table("riskDiff.txt", header=T, sep="\t", check.names=F)     #读取输入文件

#提取差异基因的名称,将基因名字转换为基因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]#去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#KEGG富集分析
R.utils::setOption("clusterProfiler.download.method","auto")
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#保存显著富集的结果
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#定义显示通路的数目
showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

#柱状图
pdf(file="KEGGbarplot.pdf", width=9, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel)
dev.off()

#气泡图
pdf(file="KEGGbubble.pdf", width = 9, height = 7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel)
dev.off()

###########绘制KEGG圈图
Pathway.col=c("#90EE90", "#E7B800", "#00AFBB")
showNum=18
data=KEGG[order(KEGG$p.adjust),]
if(nrow(KEGG)>showNum){
  data=data[1:showNum,]
}
data$Pathway="KEGG"
main.col = Pathway.col[as.numeric(as.factor(data$Pathway))]

#整理圈图数据
BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(KEGG=data$ID,start=1,end=max(BgGene))
rownames(df) = df$KEGG
bed2 = data.frame(KEGG=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(KEGG=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(KEGG=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5

#绘制圈图主体部分
pdf(file="KEGG.circlize.pdf",width=10,height=10)
par(omi=c(0.1,0.1,0.1,1.5))
circos.par(track.min=c(0.01,0.01))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                      #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                    })
circos.clear()
#绘制圈图中间的图例
middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',Pathway.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)
#绘制KEGG分类的图例
main.legend = Legend(
  labels = c("KEGG"),  type="points",pch=15,
  legend_gp = gpar(col=Pathway.col), title_position = "topcenter",
  title = "Pathway", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)
#绘制pvalue的图例
logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()



####10.1mutation####


library(maftools)       #引用包
setwd("F:/1/5LUSC/10Immune")   #设置工作目录

#读取风险文件
risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F)
outTab=risk[,c(1, ncol(risk))]
colnames(outTab)=c("Tumor_Sample_Barcode", "Risk")
write.table(outTab, file="ann.txt", sep="\t", quote=F, row.names=F)

#读取基因突变的文件
geneNum=20
geneMut=read.table("geneMut.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneMut)[1:geneNum]

#注释的颜色
ann_colors=list()
col=c("Orange2","DarkOrchid")
names(col)=c("low", "high")
ann_colors[["Risk"]]=col

#绘制低风险组瀑布图
pdf(file="low.pdf", width=6, height=6)
maf=read.maf(maf="low.maf", clinicalData="ann.txt")
oncoplot(maf=maf, clinicalFeatures="Risk", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()

#绘制高风险组瀑布图
pdf(file="high.pdf", width=6, height=6)
maf=read.maf(maf="high.maf", clinicalData="ann.txt")
oncoplot(maf=maf, clinicalFeatures="Risk", genes=gene, annotationColor=ann_colors, keepGeneOrder=T)
dev.off()


####10.2TIDE####


#准备TIDE输入文件
setwd("F:/1/5LUSC/10Immune")
pairFile="TCGA_LUSC.txt"   

#读取表达文件，并对输入文件整理
rt=read.table(pairFile, header=T, sep="\t", check.names=F)
#rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
#exp=log2(exp+1)
#使用apply函数
#exp = as.numeric(unlist(exp))

Expr <- t(apply(exp, 1, function(x)x-(mean(x)))) 
write.table(data.frame(ID=rownames(Expr),Expr,check.names = F),'tcga_normalize.txt', sep="\t", quote=F, row.names = TRUE)


#准备TIDE输入文件
setwd("F:/1/5LUSC/10Immune")
pairFile="TCGA_LUSC.txt"   

#读取表达文件，并对输入文件整理
rt=read.table(pairFile, header=T, sep="\t", check.names=F)
#rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
exp=log2(exp+1)
#使用apply函数
#exp = as.numeric(unlist(exp))

Expr <- t(apply(exp, 1, function(x)x-(mean(x)))) 
write.table(data.frame(ID=rownames(Expr),Expr,check.names = F),'tcga_normalize_log.txt', sep="\t", quote=F, row.names = TRUE)


#引用包
library(limma)
library(ggpubr)
tideFile="TIDE_normalize_log.txt"          #TIDE的打分文件
riskFile="risk.all.txt"      #风险文件
setwd("F:/1/5LUSC/10Immune")   #设置工作目录


#读取TIDE数据
tide=read.table(tideFile, header=T, sep="\t", check.names=F, row.names=1)
group=sapply(strsplit(row.names(tide),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
tide=tide[group==0,,drop=F]

tide[,2] = gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(tide))

#dim(tide)
tide1=tide[!duplicated(tide$V2),]

rownames(tide1) = tide1[,2]

tide = as.data.frame(tide1[,-2])
rownames(tide) = rownames(tide1)
colnames(tide) = "TIDE"

#tide = tide1[,-2]



#row.names(tide)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(tide))
tide=avereps(tide)

#读取风险数据文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
sameSample=intersect(row.names(tide), row.names(risk))
tide=tide[sameSample, , drop=F]
risk=risk[sameSample, "risk", drop=F]
data=cbind(tide, risk)

#设置比较组
data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$risk))
data$risk=factor(data$risk, levels=c("Low-risk", "High-risk"))
group=levels(factor(data$risk))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#绘制小提琴图
gg1=ggviolin(data, x="risk", y="TIDE", fill = "risk", 
             xlab="", ylab="TIDE",
             palette=c("DarkOrchid","Orange2"),
             legend.title="Risk",
             add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons)

#输出图形	
pdf(file="TIDE.pdf", width=5, height=4.5)
print(gg1)
dev.off()


jco <- c("DarkOrchid","Orange2")
gg1= ggplot(data = data,aes(x = risk, y = TIDE, fill = risk))+
  scale_fill_manual(values = jco[2:1]) + 
  geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = TRUE, outlier.size = -1, color="black", lwd=0.8, alpha = 0.7)+ # 背景色透明化
  geom_point(shape = 21, size=2, position = position_jitterdodge(), color="black", alpha=1)+ # 边框线黑色
  theme_classic() +
  ylab(expression("TIDE")) +
  xlab("")  +
  stat_compare_means(comparisons = my_comparisons) + 
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=0.2), # 原图有框
    axis.ticks = element_line(size=0.2,color="black"),
    axis.ticks.length = unit(0.2,"cm"),
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10))+
  theme(axis.text= element_text(face = "bold.italic",colour = "#441718",size = 16),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        #axis.ticks = element_blank(),
        axis.title = element_text(face = "bold.italic",colour = "#441718",size = 16),#Font face ("plain", "italic", "bold", "bold.italic")
        #axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_blank(),
        plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16),
        #title = element_text(face = "bold.italic",colour = "#441718",size = 13),
        legend.text = element_text(face ="bold.italic"),
        panel.border = element_rect(fill=NA,color="#35A79D",size=1.5,linetype="solid"),
        panel.background = element_rect(fill = "#F1F6FC"),
        panel.grid.major = element_line(color = "#CFD3D6", size =.5,linetype ="dotdash" ),
        legend.title = element_text(face ="bold.italic",size = 13)
        #legend.position = "none",
        #axis.title = element_blank()
  )

#输出图形	
pdf(file="TIDE1.pdf", width=5, height=4.5)
print(gg1)
dev.off()


#引用包
library(survival)
library(survminer)
TIDEFile="TIDE_normalize.txt"            #肿瘤突变负荷文件
riskFile="risk.all.txt"      #风险文件
setwd("F:/1/5LUSC/10Immune")   #修改工作目录

#读取输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #读取风险文件
TIDE=read.table(TIDEFile, header=T, sep="\t", check.names=F, row.names=1)      #读取TIDE数据文件

#合并数据
TIDE[,2] = gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(TIDE))
#dim(TIDE)
TIDE1=TIDE[!duplicated(TIDE$V2),]
rownames(TIDE1) = TIDE1[,2]
TIDE = as.data.frame(TIDE1[,-2])
rownames(TIDE) = rownames(TIDE1)
colnames(TIDE) = "TIDE"
#row.names(TIDE)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(TIDE))
sameSample=intersect(row.names(TIDE), row.names(risk))
TIDE=TIDE[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(risk, TIDE)

#获取肿瘤突变负荷最优cutoff
res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("TIDE"))
cutoff=as.numeric(res.cut$cutpoint[1])
TIDEType=ifelse(data[,"TIDE"]<=cutoff, "L-TIDE", "H-TIDE")
scoreType=ifelse(data$risk=="low", "low risk", "high risk")
mergeType=paste0(TIDEType, "+", scoreType)

#定义生存分析函数
bioSurvival=function(surData=null, outFile=null){
  diff=survdiff(Surv(futime, fustat) ~ group, data=surData)
  length=length(levels(factor(surData[,"group"])))
  pValue=1-pchisq(diff$chisq, df=length-1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ group, data = surData)
  #print(surv_median(fit))
  
  #绘制生存曲线
  bioCol=c("DarkOrchid","Orange2","NavyBlue","MediumSeaGreen","Firebrick3")
  bioCol=bioCol[1:length]
  surPlot=ggsurvplot(fit, 
                     data=surData,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title="",
                     legend.labs=levels(factor(surData[,"group"])),
                     font.legend=10,
                     legend = c(0.8, 0.8),
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette = bioCol,
                     surv.median.line = "hv",
                     risk.table=F,
                     cumevents=F,
                     risk.table.height=.25)
  #输出图形
  pdf(file=outFile, onefile = FALSE, width=5.5, height=4.8)
  print(surPlot)
  dev.off()
}

#绘制肿瘤突变负荷的生存曲线
data$group=TIDEType
bioSurvival(surData=data, outFile="TIDE.survival.pdf")

#绘制肿瘤突变负荷联合高低风险的生存曲线
data$group=mergeType
bioSurvival(surData=data, outFile="TIDE-risk.survival.pdf")


####10.3TMB####


#引用包
library(limma)
library(ggpubr)
setwd("F:/1/5LUSC/10Immune")     #设置工作目录

#读取肿瘤突变负荷文件
tmb=read.table("TMB.txt", header=T, sep="\t", check.names=F, row.names=1)

#读取风险数据文件
risk=read.table("risk.all.txt", header=T, sep="\t", check.names=F, row.names=1)

#合并数据
sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(tmb, risk)
data$TMB=log2(data$TMB+1)

#设置比较组
data$risk=ifelse(data$risk=="high", "High-risk", "Low-risk")
group=levels(factor(data$risk))
data$risk=factor(data$risk, levels=c("Low-risk", "High-risk"))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#绘制小提琴图
boxplot=ggviolin(data, x="risk", y="TMB", fill="risk",
                 xlab="",
                 ylab="Tumor mutation burden (log2)",
                 legend.title="",
                 palette = c("DodgerBlue1","Firebrick2"),
                 add = "boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons)

#输出图片
pdf(file="riskTMB.pdf", width=5, height=4.5)
print(boxplot)
dev.off()


jco <- c("Firebrick2","DodgerBlue1")
gg1= ggplot(data = data,aes(x = risk, y = TMB, fill = risk))+
  scale_fill_manual(values = jco[2:1]) + 
  geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = TRUE, outlier.size = -1, color="black", lwd=0.8, alpha = 0.7)+ # 背景色透明化
  geom_point(shape = 21, size=2, position = position_jitterdodge(), color="black", alpha=1)+ # 边框线黑色
  theme_classic() +
  ylab(expression("TMB")) +
  xlab("")  +
  stat_compare_means(comparisons = my_comparisons) + 
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=0.2), # 原图有框
    axis.ticks = element_line(size=0.2,color="black"),
    axis.ticks.length = unit(0.2,"cm"),
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10))

#输出图片
pdf(file="riskTMB1.pdf", width=5, height=4.8)
print(gg1)
dev.off()


####10.4TMBsurv####
#引用包
library(survival)
library(survminer)
tmbFile="TMB.txt"            #肿瘤突变负荷文件
riskFile="risk.all.txt"      #风险文件
setwd("F:/1/5LUSC/10Immune")    #修改工作目录

#读取输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #读取风险文件
tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)      #读取TMB数据文件

#合并数据
sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
data=cbind(risk, tmb)

#获取肿瘤突变负荷最优cutoff
res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("TMB"))
cutoff=as.numeric(res.cut$cutpoint[1])
tmbType=ifelse(data[,"TMB"]<=cutoff, "L-TMB", "H-TMB")
scoreType=ifelse(data$risk=="low", "low risk", "high risk")
mergeType=paste0(tmbType, "+", scoreType)

#定义生存分析函数
bioSurvival=function(surData=null, outFile=null){
  diff=survdiff(Surv(futime, fustat) ~ group, data=surData)
  length=length(levels(factor(surData[,"group"])))
  pValue=1-pchisq(diff$chisq, df=length-1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ group, data = surData)
  #print(surv_median(fit))
  
  #绘制生存曲线
  bioCol=c("Firebrick2","DodgerBlue1","NavyBlue","MediumSeaGreen","Firebrick3")
  bioCol=bioCol[1:length]
  surPlot=ggsurvplot(fit, 
                     data=surData,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title="",
                     legend.labs=levels(factor(surData[,"group"])),
                     font.legend=10,
                     legend = c(0.8, 0.8),
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette = bioCol,
                     surv.median.line = "hv",
                     risk.table=F,
                     cumevents=F,
                     risk.table.height=.25)
  #输出图形
  pdf(file=outFile, onefile = FALSE, width=5.5, height=4.8)
  print(surPlot)
  dev.off()
}

#绘制肿瘤突变负荷的生存曲线
data$group=tmbType
bioSurvival(surData=data, outFile="TMB.survival.pdf")

#绘制肿瘤突变负荷联合高低风险的生存曲线
data$group=mergeType
bioSurvival(surData=data, outFile="TMB-risk.survival.pdf")



####10.5immuneCor####


#引用包
library(limma)
library(scales)
library(ggplot2)
library(ggtext)
riskFile="risk.all.txt"      #风险输入文件
immFile="infiltration_estimation_for_tcga.csv"     #免疫细胞浸润文件
setwd("F:/1/5LUSC/10Immune")     #设置工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取免疫细胞浸润文件
immune=read.csv(immFile, header=T, sep=",", check.names=F, row.names=1)
immune=as.matrix(immune)
rownames(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(immune))
immune=avereps(immune)

#对风险文件和免疫细胞浸润文件取交集，得到交集样品
sameSample=intersect(row.names(risk), row.names(immune))
risk=risk[sameSample, "riskScore"]
immune=immune[sameSample,]

#对风险打分和免疫细胞进行相关性分析
x=as.numeric(risk)
outTab=data.frame()
for(i in colnames(immune)){
  y=as.numeric(immune[,i])
  corT=cor.test(x, y, method="spearman")
  cor=corT$estimate
  pvalue=corT$p.value
  if(pvalue<0.05){
    outTab=rbind(outTab,cbind(immune=i, cor, pvalue))
  }
}
#输出相关性结果
write.table(file="corResult.txt", outTab, sep="\t", quote=F, row.names=F)

#绘制气泡图
corResult=read.table("corResult.txt", head=T, sep="\t")
corResult$Software=sapply(strsplit(corResult[,1],"_"), '[', 2)
corResult$Software=factor(corResult$Software,level=as.character(unique(corResult$Software[rev(order(as.character(corResult$Software)))])))
b=corResult[order(corResult$Software),]
b$immune=factor(b$immune,levels=rev(as.character(b$immune)))
colslabels=rep(hue_pal()(length(levels(b$Software))),table(b$Software))     #定义颜色
pdf(file="cor.pdf", width=10, height=10)       #保存图片
ggplot(data=b, aes(x=cor, y=immune, color=Software))+
  labs(x="Correlation coefficient",y="Immune cell")+
  geom_point(size=4.1)+
  theme(panel.background=element_rect(fill="white",size=1,color="black"),
        panel.grid=element_line(color="grey75",size=0.5),
        axis.ticks = element_line(size=0.5),
        axis.text.y = ggtext::element_markdown(colour=rev(colslabels)))
dev.off()


####10.6immFunction####


#引用包
library(limma)
library(GSVA)
library(GSEABase)
library(ggpubr)
library(reshape2)
expFile="TCGA_LUSC.txt"             #表达输入文件
gmtFile="immune.gmt"             #免疫功能数据集文件
riskFile="risk.all.txt"              #风险文件
socreFile="immFunScore.txt"      #免疫功能打分的输出文件
setwd("F:/1/5LUSC/10Immune")      #设置工作目录

#读取表达输入文件，并对输入文件处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0,]

#读取数据集文件
geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#ssgsea分析
ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#定义ssGSEA score矫正函数
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
#对ssGSEA score进行矫正
data=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(data), data)
write.table(ssgseaOut, file=socreFile, sep="\t", quote=F, col.names=F)

#去除正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=t(data[,group==0])
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)

#读取风险文件
risk=read.table(riskFile,header=T,sep="\t",row.names=1,check.names=F)

#合并数据
sameSample=intersect(row.names(data),row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,"risk",drop=F]
rt1=cbind(data, risk)

#对免疫相关功能绘制箱线图
data=melt(rt1,id.vars=c("risk"))
colnames(data)=c("Risk","Type","Score")
data$Risk=factor(data$Risk, levels=c("low","high"))
p=ggboxplot(data, x="Type", y="Score", color = "Risk",
            ylab="Score",add = "none",xlab="",palette = c("Orange2","DarkOrchid") )
p=p+rotate_x_text(50)
p=p+stat_compare_means(aes(group=Risk),symnum.s=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")

#输出图片文件
pdf(file="immFunction.pdf", width=10, height=5)
print(p)
dev.off()


jco <- c("DarkOrchid","Orange2")

boxplot=ggplot(data = data,aes(x = Type, y = Score, fill = Risk))+
  scale_fill_manual(values = jco[2:1]) + 
  geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = TRUE, outlier.size = -1, color="black", lwd=0.8, alpha = 0.7)+
  geom_point(shape = 21, size=0.5, position = position_jitterdodge(), color="black", alpha=0.05,)+
  theme_classic() +
  ylab(expression("Score")) +
  xlab("")  +
  
  rotate_x_text(50)+
  stat_compare_means(aes(group=Risk),method = "wilcox.test",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")+
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=0.2), # 原图有框
    axis.ticks = element_line(size=0.2,color="black"),
    axis.ticks.length = unit(0.2,"cm"),
    #legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10))+
  theme(axis.text= element_text(face = "bold.italic",colour = "#441718",size = 16),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        #axis.ticks = element_blank(),
        axis.title = element_text(face = "bold.italic",colour = "#441718",size = 16),#Font face ("plain", "italic", "bold", "bold.italic")
        #axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_blank(),
        plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16),
        #title = element_text(face = "bold.italic",colour = "#441718",size = 13),
        legend.text = element_text(face ="bold.italic"),
        panel.border = element_rect(fill=NA,color="#35A79D",size=1.5,linetype="solid"),
        panel.background = element_rect(fill = "#F1F6FC"),
        panel.grid.major = element_line(color = "#CFD3D6", size =.5,linetype ="dotdash" ),
        legend.title = element_text(face ="bold.italic",size = 13)
        #legend.position = "none",
        #axis.title = element_blank()
  )


#输出图片
pdf(file="immFunction.boxplot1.pdf", width=11, height=5.5)
print(boxplot)
dev.off()


####10.7checkpoint####


#引用包
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
expFile="TCGA_LUSC.txt"      #表达输入文件
riskFile="risk.all.txt"       #风险输入文件
geneFile="gene.txt"       #免疫检查点的基因文件
setwd("F:/1/5LUSC/10Immune")    #设置工作目录

#读取基因表达文件,并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#读取基因文件
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data),as.vector(gene[,1]))
data=t(data[sameGene,])
data=log2(data+1)

#删除正常样品
group=sapply(strsplit(row.names(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",row.names(data))
data=avereps(data)

#合并数据
risk=read.table(riskFile, sep="\t", header=T, check.names=F, row.names=1)
sameSample=intersect(row.names(data),row.names(risk))
rt1=cbind(data[sameSample,],risk[sameSample,])
rt1=rt1[,c(sameGene,"risk")]

#提取显著差异的基因
sigGene=c()
for(i in colnames(rt1)[1:(ncol(rt1)-1)]){
  if(sd(rt1[,i])<0.001){next}
  wilcoxTest=wilcox.test(rt1[,i] ~ rt1[,"risk"])
  pvalue=wilcoxTest$p.value
  if(wilcoxTest$p.value<0.05){
    sigGene=c(sigGene, i)
  }
}
sigGene=c(sigGene, "risk")
rt1=rt1[,sigGene]

#把数据转换成ggplot2输入文件
rt1=melt(rt1,id.vars=c("risk"))
colnames(rt1)=c("risk","Gene","Expression")

#设置比较组
group=levels(factor(rt1$risk))
rt1$risk=factor(rt1$risk, levels=c("low","high"))
comp=combn(group,2)
my_comparisons=list()
for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}

#绘制箱线图
boxplot=ggboxplot(rt1, x="Gene", y="Expression", fill="risk",
                  xlab="",
                  ylab="Gene expression",
                  legend.title="Risk",
                  width=0.8,
                  palette = c("DodgerBlue1","Firebrick2") )+
  rotate_x_text(50)+
  stat_compare_means(aes(group=risk),
                     method="wilcox.test",
                     symnum.s=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")

#输出图片
pdf(file="checkpoint.diff.pdf", width=10, height=4.5)
print(boxplot)
dev.off()



jco <- c("DarkOrchid","Orange2")

boxplot=ggplot(data = rt1,aes(x = Gene, y = Expression, fill = risk))+
  scale_fill_manual(values = jco[2:1]) + 
  geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = TRUE, outlier.size = -1, color="black", lwd=0.8, alpha = 0.7)+
  geom_point(shape = 21, size=0.5, position = position_jitterdodge(), color="black", alpha=0.05,)+
  theme_classic() +
  ylab(expression("Gene expression")) +
  xlab("")  +
  
  rotate_x_text(50)+
  stat_compare_means(aes(group=risk),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")+
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=0.2), # 原图有框
    axis.ticks = element_line(size=0.2,color="black"),
    axis.ticks.length = unit(0.2,"cm"),
    #legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10))+
  theme(axis.text= element_text(face = "bold.italic",colour = "#441718",size = 16),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        #axis.ticks = element_blank(),
        axis.title = element_text(face = "bold.italic",colour = "#441718",size = 16),#Font face ("plain", "italic", "bold", "bold.italic")
        #axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_blank(),
        plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16),
        #title = element_text(face = "bold.italic",colour = "#441718",size = 13),
        legend.text = element_text(face ="bold.italic"),
        panel.border = element_rect(fill=NA,color="#35A79D",size=1.5,linetype="solid"),
        panel.background = element_rect(fill = "#F1F6FC"),
        panel.grid.major = element_line(color = "#CFD3D6", size =.5,linetype ="dotdash" ),
        legend.title = element_text(face ="bold.italic",size = 13)
        #legend.position = "none",
        #axis.title = element_blank()
  )

#输出图片
pdf(file="checkpoint.diff1.pdf", width=10, height=4.5)
print(boxplot)
dev.off()


####10.8pancancer####


#使用国内镜像安装包
options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
#加载包
library(ggplot2)
library(data.table)
library(cowplot)
library(ggpubr)
library(GSVA)
library(SimDesign)
library(tidyr)
setwd("F:/1/5LUSC/10Immune/pancancer")
#自定义函数将gmt文件读取为list
gmt2list <- function(annofile){
  if (!file.exists(annofile)) {
    stop("There is no such gmt file.")
  }
  
  if (tools::file_ext(annofile) == "xz") {
    annofile <- xzfile(annofile)
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
    close(annofile)
  } else if (tools::file_ext(annofile) == "gmt") {
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
  } else {
    stop ("Only gmt and gmt.xz are accepted for gmt2list")
  }
  
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  
  annoList <- lapply(y, `[`, c(-1,-2))
}

# 读取风险基因以及对应系数（来自原文补充材料表格Table S6）
risk.coeff <- read.table("multiCox.txt",sep = "\t", row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
# 读取肿瘤注释文件
rawAnno <- read.delim("merged_sample_quality_annotations.tsv",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
rawAnno$simple_barcode <- substr(rawAnno$aliquot_barcode,1,15)
samAnno <- rawAnno[!duplicated(rawAnno$simple_barcode),c("cancer type", "simple_barcode")]
samAnno <- samAnno[which(samAnno$`cancer type` != ""),]
write.table(samAnno,"output_simple_sample_annotation.txt",sep = "\t",row.names = F,col.names = T,quote = F)
# 快速读取表达谱数据并做数据预处理
expr <- fread("EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv",sep = "\t",stringsAsFactors = F,check.names = F,header = T)
expr <- as.data.frame(expr); rownames(expr) <- expr[,1]; expr <- expr[,-1]
gene <- sapply(strsplit(rownames(expr),"|",fixed = T), "[",1)
expr$gene <- gene
expr <- expr[!duplicated(expr$gene),]
rownames(expr) <- expr$gene; expr <- expr[,-ncol(expr)]
expr[expr < 0] <- 0 # 对于这份泛癌数据，将略小于0的数值拉到0，否则不能取log（其他途径下载的泛癌数据可能不需要此操作）
colnames(expr) <- substr(colnames(expr),1,15)
gc()
# 去掉对于风险基因存在NA值的样本
expr.sub <- expr[risk.coeff$Gene, ] # 提取仅有风险相关基因的表达谱子集
expr.sub <- as.data.frame(t(na.omit(t(expr.sub)))) # 对列做去空值，而非对行做
keepSam <- colnames(expr.sub) # 提取被保留的样本
expr <- expr[,keepSam] # 重构表达谱

# 读取生存数据(虽然在本代码中没有用到，但是原文使用的样本是具有生存数据的)
surv <- read.delim("Survival_SupplementalTable_S1_20171025_xena_sp", sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) 
# 确定肿瘤样本以及对应肿瘤类型
sam <- samAnno[which(samAnno$`cancer type` != "LAML"),"simple_barcode"] # 去掉白血病样本
comsam <- intersect(intersect(colnames(expr), sam), rownames(surv)) # 得到与表达谱以及生存的共有样本
tumsam <- comsam[substr(comsam,14,14) == "0"] # 仅提取肿瘤样本
tumAnno <- samAnno[which(samAnno$simple_barcode %in% tumsam),] # 获取这些肿瘤样本的注释信息
tumAnno <- tumAnno[order(tumAnno$`cancer type`),] # 根据肿瘤类型排序
tumors <- unique(tumAnno$`cancer type`) # 得到32个肿瘤

# 在所有样本中计算CAFRGs得分(在本代码中仅仅是为了确定根据cox-based CAFRGs score确定肿瘤的level)
CAFRGs.score <- list() # 初始化列表
CAFRGs.mean <- c() # 初始化得分均值向量
outTab <- NULL

colnames(risk.coeff)[1] = "Gene"

colnames(risk.coeff)[2] = "Coefficient"
for (i in tumors) {
  sam <- tumAnno[which(tumAnno$`cancer type` == i),"simple_barcode"] # 提取当前肿瘤类型的肿瘤样本
  expr.sub <- log2(expr[risk.coeff$Gene,sam] + 1) # 提取表达谱子集并对数化
  CAFRGs <- scale(apply(expr.sub,2,function(x) {x %*% risk.coeff$Coefficient})) # 计算经过z-score的CAFRGs得分
  CAFRGs.score[[i]] <- CAFRGs
  CAFRGs.mean <- c(CAFRGs.mean, mean(CAFRGs))
  outTab <- rbind.data.frame(outTab, # 保存得分的计算结果
                             data.frame(tumor = i, # 肿瘤类型
                                        CAFRGs = as.numeric(CAFRGs), # 当前得分
                                        row.names = sam,
                                        stringsAsFactors = F),
                             stringsAsFactors = F)
}
sapply(CAFRGs.score, range) # 不存在空值
write.table(outTab, file = "output_CAFRGs score of all tumor sample acCAFs 32 tumor types.txt",sep = "\t",row.names = T,col.names = F,quote = F)
names(CAFRGs.mean) <- tumors
CAFRGs.mean <- sort(CAFRGs.mean, decreasing = T) # 根据均值对肿瘤进行排序
tumor.level <- names(CAFRGs.mean) # 将排序结果作为肿瘤因子的等级

# 在所有样本中通过z-score计算致癌通路以及CAFRGs得分（注意此时CAFRGs得分不再是由cox系数计算，而是由zscore算法下的单样本富集得到）
oncosig <- gmt2list("oncogenic.gmt") # 将原文补充材料S4以及S6的基因制作成gmt文件，并将gmt文件读取为list
oncosig[[4]][1:13] = risk.coeff$Gene
oncosig[[4]][14:19] = ""
oncosig <- sapply(oncosig, function(x) setdiff(x,"")) # 去掉list中的空值
zscore.list <- list()
outSig <- NULL
for (i in tumors) {
  message(i)
  sam <- tumAnno[which(tumAnno$`cancer type` == i),"simple_barcode"] # 提取当前肿瘤类型的肿瘤样本
  expr.sub <- log2(expr[,sam] + 1) # 提取表达谱子集并对数化
  zscore.list[[i]] <- quiet(gsva(as.matrix(expr.sub), gset.idx.list = oncosig, method = "zscore")) # 方法选择zscore
  outSig <- rbind.data.frame(outSig, # 保存得分的计算结果
                             cbind.data.frame(tumor = i,
                                              as.data.frame(t(zscore.list[[i]]))),
                             stringsAsFactors = F)
}
write.table(outSig, file = "output_oncogenic and CAFRGs score of all tumor sample acCAFs 32 tumor types.txt",sep = "\t",row.names = T,col.names = F,quote = F)



# 设置颜色
mycol <- c("#A6CEE3",
           "#1F78B4",
           "#B2DF8A",
           "#33A02C",
           "#FB9A99",
           "#E31A1C",
           "#FDBF6F",
           "#FF7F00",
           "#CAB2D6",
           "#6A3D9A",
           "#B15928",
           "#8DD3C7",
           "#BEBADA",
           "#FB8072",
           "#80B1D3",
           "#FDB462",
           "#B3DE69",
           "#FCCDE5",
           "#D9D9D9",
           "#BC80BD",
           "#CCEBC5",
           "#FFED6F",
           "#8C510A",
           "#BF812D",
           "#DFC27D",
           "#F6E8C3",
           "#80CDC1",
           "#35978F",
           "#01665E",
           "#003C30",
           "#8E0152",
           "#C51B7D")

# 制作绘图数据并绘图
plotdata <- outSig
plotdata <- gather(plotdata, oncogenic, zscore, Angiogenesis:`Cell cycle`, factor_key=TRUE)
plotdata$tumor <- factor(plotdata$tumor, levels = tumor.level)

p1 <- ggplot(data = plotdata, aes(x = zscore, y = CAFRGs)) + 
  geom_point(aes(color=tumor),size=1.5,alpha = 0.5) +
  scale_color_manual(values = mycol) + 
  geom_smooth(method = "lm", se = FALSE) +
  geom_ribbon(stat = "smooth",method = "lm",se = TRUE,alpha = 0,linetype = "dashed") + 
  xlab("Oncogenic (z-score)") + ylab("CAFRGs (z-score)") + 
  stat_cor(method = "pearson", label.x = -40, label.y = 10) + 
  facet_wrap(.~oncogenic, nrow = 1) + 
  theme_bw() + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 12, color = "black"),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "right",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))
## Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
## ℹ Please use the `linewidth` argument instead.
ggsave(filename = "correlation scatter plot of zscored oncogenic and CAFRGs in pancancer.pdf", width = 15,height = 5)
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'


tmp1 <- plotdata[which(plotdata$oncogenic == "Angiogenesis"),]
p2 <- ggplot(data = tmp1, aes(x = zscore, y = CAFRGs)) + 
  geom_point(aes(color=tumor),size=1.2,alpha = 0.5) +
  scale_color_manual(values = mycol) + 
  geom_smooth(method = "lm", se = FALSE) +
  geom_ribbon(stat = "smooth",method = "lm",se = TRUE,alpha = 0,linetype = "dashed") + 
  stat_cor(method = "pearson", label.x = -40, label.y = 10) + 
  xlab("Angiogenesis (z-score)") + ylab("CAFRGs (z-score)") + 
  facet_wrap(.~tumor, ncol = 8) + 
  theme_bw() + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 12, color = "black"),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))
ggsave(filename = "correlation scatter plot of zscored Angiogenesis and CAFRGs in pancancer.pdf", width = 15,height = 8)
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'


tmp2 <- plotdata[which(plotdata$oncogenic == "EMT"),]
p3 <- ggplot(data = tmp2, aes(x = zscore, y = CAFRGs)) + 
  geom_point(aes(color=tumor),size=1.2,alpha = 0.5) +
  scale_color_manual(values = mycol) + 
  geom_smooth(method = "lm", se = FALSE) +
  geom_ribbon(stat = "smooth",method = "lm",se = TRUE,alpha = 0,linetype = "dashed") + 
  stat_cor(method = "pearson", label.x = -40, label.y = 10) + 
  xlab("EMT (z-score)") + ylab("CAFRGs (z-score)") + 
  facet_wrap(.~tumor, ncol = 8) + 
  theme_bw() + 
  theme(axis.text.x = element_text(vjust = 0.5, size = 12, color = "black"),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))
ggsave(filename = "correlation scatter plot of zscored EMT and CAFRGs in pancancer.pdf", width = 15,height = 8)
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
p <- plot_grid(p1,p2,p3, align = "v", ncol = 1)
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## `geom_smooth()` using formula = 'y ~ x'
## Warning: Graphs cannot be vertically aligned unless the axis parameter is set.
## Placing graphs unaligned.
ggsave(filename = "combined correlation scatter plot of zscored oncogenic and CAFRGs in pancancer.pdf", width = 18,height = 23)


####10.9oncoPredict####


#引用包
library(limma)
library(oncoPredict)
library(parallel)
set.seed(12345)

expFile="TCGA_LUSC.txt"     #表达数据文件
setwd("F:/1/5LUSC/10Immune")     #设置工作目录

#读取表达输入文件,并对数据进行处理
rt = read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
#data=data[rowMeans(data)>0.5,]
colnames(data)=gsub("(.*?)\\_(.*?)", "\\2", colnames(data))

#读取药物敏感性文件
GDSC2_Expr=readRDS(file='GDSC2_Expr.rds')
GDSC2_Res=readRDS(file = 'GDSC2_Res.rds')
GDSC2_Res=exp(GDSC2_Res) 

#药物敏感性
calcPhenotype(trainingExprData = GDSC2_Expr,    #train组的表达数据
              trainingPtype = GDSC2_Res,        #train组的药物数据
              testExprData = data,              #test组的表达数据
              batchCorrect = 'eb',              #批次矫正的方法
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,      #去除波动小的基因
              minNumSamples = 10,               #最小的样品数目
              printOutput = TRUE,               #是否输出结果
              removeLowVaringGenesFrom = 'rawData')

####10.10boxplot####


#引用包
library(limma)
library(ggplot2)
library(ggpubr)
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)

pFilter=0.05                      #pvalue过滤条件
riskFile="risk.all.txt"            #风险文件
drugFile="DrugPredictions.csv"     #药物敏感性文件
setwd("F:/1/5LUSC/10Immune")      #设置工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取药物敏感性文件
senstivity=read.csv(drugFile, header=T, sep=",", check.names=F, row.names=1)
colnames(senstivity)=gsub("(.*)\\_(\\d+)", "\\1", colnames(senstivity))
#dim(senstivity)
senstivity[,dim(senstivity)[2]+1] = gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(senstivity))
#dim(senstivity)

#dim(tide)
senstivity1=senstivity[!duplicated(senstivity$V199),]
rownames(senstivity1) = senstivity1[,dim(senstivity)[2]]
senstivity = as.data.frame(senstivity1[,-199])
rownames(senstivity) = rownames(senstivity1)
#colnames(senstivity) = "TIDE"
#tide = tide1[,-2]
#colnames(senstivity)
#row.names(senstivity)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(senstivity))

#数据合并
sameSample=intersect(row.names(risk), row.names(senstivity))
risk=risk[sameSample, "risk",drop=F]
senstivity=senstivity[sameSample,,drop=F]
senstivity[is.na(senstivity)] = 0
senstivity = log2(senstivity++1)
rt=cbind(risk, senstivity)

#设置比较组
rt$risk=factor(rt$risk, levels=c("low", "high"))
type=levels(factor(rt[,"risk"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}


#提取显著差异的基因
sigGene=c()
for(i in colnames(rt)[2:(ncol(rt))]){
  if(sd(rt[,i])<0.05){next}
  wilcoxTest=wilcox.test(rt[,i] ~ rt[,"risk"])
  pvalue=wilcoxTest$p.value
  if(wilcoxTest$p.value<0.05){
    sigGene=c(sigGene, i)
  }
}
sigGene=c(sigGene, "risk")
rt=rt[,sigGene]

colnames(rt)
rt = rt[,c("Camptothecin","Vinblastine","Cisplatin","Cytarabine","Docetaxel","Gefitinib","Vorinostat","Nilotinib","Olaparib","Afatinib","Doramapimod","Alisertib","Palbociclib","Dactolisib","Pictilisib","5-Fluorouracil","Dasatinib","Paclitaxel","Crizotinib","Rapamycin","Sorafenib","Irinotecan","Oxaliplatin","Tozasertib","Erlotinib","Gemcitabine","Bortezomib","Tamoxifen","Fulvestrant","Daporinad","Talazoparib","Dabrafenib","Temozolomide","Ruxolitinib","Linsitinib","Epirubicin","Cyclophosphamide","Pevonedistat","Sapitinib","Uprosertib","Lapatinib","Alpelisib","Taselisib","Leflunomide","Entinostat","Ribociclib","Picolinici-acid","Selumetinib","Ibrutinib","Zoledronate","Oxaliplatin.1","Carmustine","Mitoxantrone","Dactinomycin","Fulvestrant.1","Vincristine","Docetaxel.1","Podophyllotoxin bromide","Dihydrorotenone","Gallibiscoquinazole","Elephantin","Sinularin","Sabutoclax","Buparlisib","Ulixertinib","Venetoclax","Dactinomycin.1","Afuresertib","Osimertinib","Cediranib","Ipatasertib","Savolitinib","Sepantronium bromide","Foretinib","Pyridostatin","Ulixertinib.1","Vinorelbine","Uprosertib.1","risk")]
#rt = rt[,c("AZD7762","SB216763","KU-55933","PLX-4720","NU7441","Wee1 Inhibitor","Nutlin-3a (-)","Mirin","PD173074","ZM447439","RO-3306","MK-2206","BI-2536","BMS-536924","GSK1904529A","MK-1775","SB505124","EPZ004777","YK-4-279","BMS-345541","IAP_5620","AZD1208","LCL161","EPZ5676","SCH772984","IWP-2","OSI-027","VE-822","WZ4003","CZC24832","AZD5582","PFI3","PCI-34051","Wnt-C59","I-BET-762","RVX-208","OTX015","GSK343","ML323","AGI-6780","AZD5153","Eg5_9814","ERK_6604","IRAK4_4710","JAK1_8709","AZD5991","PAK_5339","TAF1_5496","ULK1_4989","VSP34_8731","IGF1R_3801","AZD4547","LY2109761","OF-1","KRAS (G12C) Inhibitor-12","MG-132","BDP-00009066","AGI-5198","AZD3759","AZD5363","AZD6738","GDC0810","GNE-317","GSK2578215A","I-BRD9","Telomerase Inhibitor IX","MIRA-1","NVP-ADW742","P22077","MIM1","BPD-00008900","BIBR-1532","MK-8776","VX-11e","LJI308","AZ6102","GSK591","VE821","AZD6482","AT13148","BMS-754807","risk")]

#把数据转换成ggplot2输入文件
rt=melt(rt,id.vars=c("risk"))
colnames(rt)=c("risk","Gene","Expression")

#设置比较组
group=levels(factor(rt$risk))
rt$risk=factor(rt$risk, levels=c("low","high"))
comp=combn(group,2)
my_comparisons=list()
for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
rt$Expression=log2(rt$Expression+1)

#绘制箱线图
boxplot=ggboxplot(rt, x="Gene", y="Expression", fill="risk",
                  xlab="",
                  ylab="Drug Senstivity",
                  legend.title="Risk",
                  width=0.8,
                  palette = c("DodgerBlue1","Firebrick2") )+
  rotate_x_text(50)+
  stat_compare_means(aes(group=risk),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")+
  theme(axis.text= element_text(face = "bold.italic",colour = "#441718",size = 16),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        #axis.ticks = element_blank(),
        axis.title = element_text(face = "bold.italic",colour = "#441718",size = 16),#Font face ("plain", "italic", "bold", "bold.italic")
        #axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_blank(),
        plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16),
        #title = element_text(face = "bold.italic",colour = "#441718",size = 13),
        legend.text = element_text(face ="bold.italic"),
        panel.border = element_rect(fill=NA,color="#35A79D",size=1.5,linetype="solid"),
        panel.background = element_rect(fill = "#F1F6FC"),
        panel.grid.major = element_line(color = "#CFD3D6", size =.5,linetype ="dotdash" ),
        legend.title = element_text(face ="bold.italic",size = 13)
        #legend.position = "none",
        #axis.title = element_blank()
  )

#输出图片
pdf(file="drugSenstivity1.pdf", width=14, height=6)
print(boxplot)
dev.off()



jco <- c("DarkOrchid","Orange2")

boxplot=ggplot(data = rt,aes(x = Gene, y = Expression, fill = risk))+
  scale_fill_manual(values = jco[2:1]) + 
  geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = TRUE, outlier.size = -1, color="black", lwd=0.8, alpha = 0.7)+
  geom_point(shape = 21, size=0.5, position = position_jitterdodge(), color="black", alpha=0.05,)+
  theme_classic() +
  ylab(expression("The half inhibitory concentration(IC50)")) +
  xlab("")  +
  
  rotate_x_text(50)+
  stat_compare_means(aes(group=risk),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")+
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=0.2), # 原图有框
    axis.ticks = element_line(size=0.2,color="black"),
    axis.ticks.length = unit(0.2,"cm"),
    #legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10))

#输出图片
pdf(file="drugSenstivity11.pdf", width=15.5, height=5.5)
print(boxplot)
dev.off()



#引用包
library(limma)
library(ggplot2)
library(ggpubr)
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)

pFilter=0.05                      #pvalue过滤条件
riskFile="risk.all.txt"            #风险文件
drugFile="DrugPredictions.csv"     #药物敏感性文件
setwd("F:/1/5LUSC/10Immune")      #设置工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取药物敏感性文件
senstivity=read.csv(drugFile, header=T, sep=",", check.names=F, row.names=1)
colnames(senstivity)=gsub("(.*)\\_(\\d+)", "\\1", colnames(senstivity))
#dim(senstivity)
senstivity[,dim(senstivity)[2]+1] = gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(senstivity))
#dim(senstivity)

#dim(tide)
senstivity1=senstivity[!duplicated(senstivity$V199),]
rownames(senstivity1) = senstivity1[,dim(senstivity)[2]]
senstivity = as.data.frame(senstivity1[,-199])
rownames(senstivity) = rownames(senstivity1)
#colnames(senstivity) = "TIDE"
#tide = tide1[,-2]
#colnames(senstivity)
#row.names(senstivity)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(senstivity))

#数据合并
sameSample=intersect(row.names(risk), row.names(senstivity))
risk=risk[sameSample, "risk",drop=F]
senstivity=senstivity[sameSample,,drop=F]
senstivity[is.na(senstivity)] = 0
senstivity = log2(senstivity++1)
rt=cbind(risk, senstivity)

#设置比较组
rt$risk=factor(rt$risk, levels=c("low", "high"))
type=levels(factor(rt[,"risk"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}


#提取显著差异的基因
sigGene=c()
for(i in colnames(rt)[2:(ncol(rt))]){
  if(sd(rt[,i])<0.05){next}
  wilcoxTest=wilcox.test(rt[,i] ~ rt[,"risk"])
  pvalue=wilcoxTest$p.value
  if(wilcoxTest$p.value<0.05){
    sigGene=c(sigGene, i)
  }
}
sigGene=c(sigGene, "risk")
rt=rt[,sigGene]

colnames(rt)
#rt = rt[,c("Camptothecin","Vinblastine","Cisplatin","Cytarabine","Docetaxel","Gefitinib","Vorinostat","Nilotinib","Olaparib","Afatinib","Doramapimod","Alisertib","Palbociclib","Dactolisib","Pictilisib","5-Fluorouracil","Dasatinib","Paclitaxel","Crizotinib","Rapamycin","Sorafenib","Irinotecan","Oxaliplatin","Tozasertib","Erlotinib","Gemcitabine","Bortezomib","Tamoxifen","Fulvestrant","Daporinad","Talazoparib","Dabrafenib","Temozolomide","Ruxolitinib","Linsitinib","Epirubicin","Cyclophosphamide","Pevonedistat","Sapitinib","Uprosertib","Lapatinib","Alpelisib","Taselisib","Leflunomide","Entinostat","Ribociclib","Picolinici-acid","Selumetinib","Ibrutinib","Zoledronate","Oxaliplatin.1","Carmustine","Mitoxantrone","Dactinomycin","Fulvestrant.1","Vincristine","Docetaxel.1","Podophyllotoxin bromide","Dihydrorotenone","Gallibiscoquinazole","Elephantin","Sinularin","Sabutoclax","Buparlisib","Ulixertinib","Venetoclax","Dactinomycin.1","Afuresertib","Osimertinib","Cediranib","Ipatasertib","Savolitinib","Sepantronium bromide","Foretinib","Pyridostatin","Ulixertinib.1","Vinorelbine","Uprosertib.1","risk")]
rt = rt[,c("AZD7762","SB216763","KU-55933","PLX-4720","NU7441","Wee1 Inhibitor","Nutlin-3a (-)","Mirin","PD173074","ZM447439","RO-3306","MK-2206","BI-2536","BMS-536924","GSK1904529A","MK-1775","SB505124","EPZ004777","YK-4-279","BMS-345541","IAP_5620","AZD1208","LCL161","EPZ5676","SCH772984","IWP-2","OSI-027","VE-822","WZ4003","CZC24832","AZD5582","PFI3","PCI-34051","Wnt-C59","I-BET-762","RVX-208","OTX015","GSK343","ML323","AGI-6780","AZD5153","Eg5_9814","ERK_6604","IRAK4_4710","JAK1_8709","AZD5991","PAK_5339","TAF1_5496","ULK1_4989","VSP34_8731","IGF1R_3801","AZD4547","LY2109761","OF-1","KRAS (G12C) Inhibitor-12","MG-132","BDP-00009066","AGI-5198","AZD3759","AZD5363","AZD6738","GDC0810","GNE-317","GSK2578215A","I-BRD9","Telomerase Inhibitor IX","MIRA-1","NVP-ADW742","P22077","MIM1","BPD-00008900","BIBR-1532","MK-8776","VX-11e","LJI308","AZ6102","GSK591","VE821","AZD6482","AT13148","BMS-754807","risk")]

#把数据转换成ggplot2输入文件
rt=melt(rt,id.vars=c("risk"))
colnames(rt)=c("risk","Gene","Expression")

#设置比较组
group=levels(factor(rt$risk))
rt$risk=factor(rt$risk, levels=c("low","high"))
comp=combn(group,2)
my_comparisons=list()
for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
rt$Expression=log2(rt$Expression+1)

#绘制箱线图
boxplot=ggboxplot(rt, x="Gene", y="Expression", fill="risk",
                  xlab="",
                  ylab="Drug Senstivity",
                  legend.title="Risk",
                  width=0.8,
                  palette = c("DodgerBlue1","Firebrick2") )+
  rotate_x_text(50)+
  stat_compare_means(aes(group=risk),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")+
  theme(axis.text= element_text(face = "bold.italic",colour = "#441718",size = 16),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        #axis.ticks = element_blank(),
        axis.title = element_text(face = "bold.italic",colour = "#441718",size = 16),#Font face ("plain", "italic", "bold", "bold.italic")
        #axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_blank(),
        plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16),
        #title = element_text(face = "bold.italic",colour = "#441718",size = 13),
        legend.text = element_text(face ="bold.italic"),
        panel.border = element_rect(fill=NA,color="#35A79D",size=1.5,linetype="solid"),
        panel.background = element_rect(fill = "#F1F6FC"),
        panel.grid.major = element_line(color = "#CFD3D6", size =.5,linetype ="dotdash" ),
        legend.title = element_text(face ="bold.italic",size = 13)
        #legend.position = "none",
        #axis.title = element_blank()
  )

#输出图片
pdf(file="drugSenstivity2.pdf", width=14, height=6)
print(boxplot)
dev.off()


jco <- c("DarkOrchid","Orange2")

boxplot=ggplot(data = rt,aes(x = Gene, y = Expression, fill = risk))+
  scale_fill_manual(values = jco[2:1]) + 
  geom_violin(alpha=0.4, position = position_dodge(width = .75),size=0.8, color="black") + # 边框线黑色
  geom_boxplot(notch = TRUE, outlier.size = -1, color="black", lwd=0.8, alpha = 0.7)+
  geom_point(shape = 21, size=0.5, position = position_jitterdodge(), color="black", alpha=0.05,)+
  theme_classic() +
  ylab(expression("The half inhibitory concentration(IC50)")) +
  xlab("")  +
  
  rotate_x_text(50)+
  stat_compare_means(aes(group=risk),
                     method="wilcox.test",
                     symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")+
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=0.2), # 原图有框
    axis.ticks = element_line(size=0.2,color="black"),
    axis.ticks.length = unit(0.2,"cm"),
    #legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10))
#输出图片
pdf(file="drugSenstivity22.pdf", width=15.5, height=5.5)
print(boxplot)
dev.off()


####10.11桑葚图和箱线图####


#引用包
library(ggalluvial)
library(ggplot2)
library(dplyr)
library(ggpubr)

cluFile="CAFcluster.txt"      #分型的结果文件
riskFile="risk.all.txt"       #风险文件
setwd("F:/1/5LUSC/10Immune")     #设置工作目录

#读取输入文件
cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
rownames(cluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(cluster))
sameSample=intersect(row.names(cluster), row.names(risk))
data=cbind(risk[sameSample,,drop=F], cluster[sameSample,,drop=F])

#######箱线图
#设置比较组
data$CAFcluster=factor(data$CAFcluster, levels=levels(factor(data$CAFcluster)))
group=levels(factor(data$CAFcluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#定义图形颜色
bioCol=c("#DA4A35","#44B1C9","#7CC767","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$CAFcluster)))]

#绘制箱线图
boxplot=ggboxplot(data, x="CAFcluster", y="riskScore", color="CAFcluster",
                  xlab="CAFcluster",
                  ylab="Risk score",
                  legend.title="CAFcluster",
                  palette=bioCol,
                  add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons)

#输出图形
pdf(file="clusterRisk.pdf", width=5.5, height=4.5)
print(boxplot)
dev.off()

#######箱线图
#准备桑基图输入文件
rt=data[,c("CAFcluster", "risk", "fustat")]
colnames(rt)=c("CAFcluster", "Risk", "Fustat")
rt[,"Fustat"]=ifelse(rt[,"Fustat"]==0, "Alive", "Dead")
corLodes=to_lodes_form(rt, axes = 1:ncol(rt), id = "Cohort")

#绘制桑基图
pdf(file="ggalluvial.pdf", width=6, height=5.5)
mycol=rep(c("#DA4A35","#44B1C9","#7CC767","DarkOrchid","Orange2","DodgerBlue1","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) +  
  #设置线条颜色，forward说明线条颜色与前面的柱状图一致，backward说明线条颜色与后面的柱状图一致。
  geom_flow(width = 2/10,aes.flow = "forward") + 
  geom_stratum(alpha = .9,width = 2/10) +
  scale_fill_manual(values = mycol) +
  #size=3代表字体大小
  geom_text(stat = "stratum", size = 3,color="black") +
  xlab("") + ylab("") + theme_bw() + 
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + #去掉坐标轴
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_blank()) + 
  ggtitle("") + guides(fill = FALSE)                            
dev.off()

