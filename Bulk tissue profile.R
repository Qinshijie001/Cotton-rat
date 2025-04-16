library(dplyr)
library(readxl)
library(xlsx)
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(limma)
library(RColorBrewer)
library(factoextra)
library(cluster)

my36colors <-c( '#53A85F', '#E5D2DD', '#AA9A59', '#E95C59', '#D6E7A3', '#57C3F3',
                '#F4A460', '#F3B1A0', '#9ACD32', '#AB3282', '#DDA0DD', '#91D0BE', 
                '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA',
                '#58A4C3', '#E4C755', '#F7F398', '#F1BB72', '#E63863', '#E39A35',
                '#C1E6F3', '#6778AE', '#BD956A', '#B53E2B', '#712820', '#DCC1DD',
                '#EE6A50', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')
data1 <- read_excel('Tissue_special.xlsx')
data1 <- as.data.frame(data1)

exp1 <- data1[,4:24]
rownames(exp1) <- data1$gene_id
exp1 <- as.data.frame(normalizeBetweenArrays(exp1))

exp2 <- exp1[rowSums(exp1)>10,]
data2 <- subset(data1,gene_id %in% rownames(exp2))

data3 <- subset(data2,TAU_Index>0.9)
data3 <- data3[order(data3$Preferred_Sample),]
rownames(data3) <- data3$gene_id

#
data4 <- data.frame()
for (i in (names(sort(table(data3$Preferred_Sample),decreasing = T)))){
  a=subset(data3,Preferred_Sample==i)
  data4 <- rbind(data4,a)
}

data4 <- data4[,names(sort(table(data3$Preferred_Sample),decreasing = T))]
p <- pheatmap(data4,scale = 'row',cluster_rows = F,cluster_cols = F)
ggsave('special_pheatmap.pdf',p,width = 5,height = 5)

#
data6 <- data3
data6$Preferred_Sample <- factor(data6$Preferred_Sample,levels = c(names(sort(table(data3$Preferred_Sample)))))
p=ggbarplot(data6,x='Preferred_Sample',y='TAU_Index')+coord_flip()
ggsave('special_counts.pdf',p,width = 3,height = 4.5)

#
gene_anno <- read.csv('../gene_annotation/mianshu_Blast_gene_name.csv')
data1[,4:24] <- exp1
exp3 <- merge(gene_anno,data1,by.x = 'gene_id1',by.y = 'gene_id',all = T)
exp4 <- exp3
rownames(exp4) <- exp4$blast_name1
recepter <- grep("CX3CR1|HSPG|TLR4|NCL$|EGFR|IGF1R|INSR|ICAM1",rownames(exp4),ignore.case = T,value = T)
recepter
recepter <- recepter[c(2,5:11)]
recepter
exp5 <- as.data.frame(t(exp4[recepter,8:28]))
exp5 <- log2(exp5+1)
exp5 <- exp5[names(sort(table(data3$Preferred_Sample),decreasing = T)),]
p <- pheatmap(t(exp5),scale = 'row',cluster_rows = T,cluster_cols = F)
ggsave('recepter.pdf',p,width = 8,height = 3)

##
gtex <- read.csv('/Gtex_exp.csv',row.names = 1)
mianshu <- as.data.frame(read_excel('mianshu_exp.xlsx'))
rownames(mianshu) <- mianshu$gene_id
mianshu <- mianshu[,-1]

homo <- read.table('mianshu_human_data_subset.txt',header = F,sep = '\t')
homo$mianshu_id <- gsub('mianshu|','',homo$V2,fixed = T)
homo$human_id <- gsub('gene-','',homo$V4,fixed = T)
homo$mianshu_id1 <- gsub('.1','',homo$mianshu_id,fixed = T)
homo1 <- distinct(homo,by=human_id,.keep_all = T)
rownames(homo1) <- homo1$human_id

colnames(mianshu) <- c("Spleen","Stomach","Pancreas","Trachea","Blood_Vessel","Bladder","Brain","Breast","Small_Intestine","Adipose_Tissue",
                       "Liver","Bone_Marrow","Heart","Lung","Muscle","Nose","Ovary","Kidney","Skin","Thymus","Tongue")
tissue_jiaoji <- intersect(colnames(mianshu),colnames(gtex))
mianshu_tissue <- mianshu[,tissue_jiaoji]
human_tissue <- gtex[,tissue_jiaoji] 
name_jiao <- intersect(rownames(human_tissue),homo1$human_id)

homo2 <- homo1[name_jiao,]
human_tissue1 <- human_tissue[name_jiao,]
mianshu_tissue1 <- mianshu_tissue[homo2$mianshu_id1,]
human_tissue2 <- as.data.frame(t(human_tissue1))  
mianshu_tissue2 <- as.data.frame(t(mianshu_tissue1))  

human_mianshu <- data.frame()
for (i in rownames(human_tissue2)){
  temp1 <- human_tissue2[i,]
  temp2 <- mianshu_tissue2[i,]
  colnames(temp2) <- colnames(temp1)
  temp3 <- rbind(temp1,temp2)
  human_mianshu <- rbind(human_mianshu,temp3)
}

human_mianshu1 <- as.data.frame(apply(human_mianshu,2,as.numeric))
rownames(human_mianshu1) <- rownames(human_mianshu)
human_mianshu1 <- log2(human_mianshu1+1)
human_mianshu2 <- as.data.frame(normalizeBetweenArrays(t(human_mianshu1)))
human_mianshu3 <- as.data.frame(t(human_mianshu1))

#
library(sva)
#human_mianshu4 <- human_mianshu2[,grep('Brain|Heart|Spleen|Muscle|Liver|Small|Ovary|Lung|Stomach|Bladder|Kidney|Pancreas',colnames(human_mianshu2),value = T)]
data_homo <- human_mianshu3[,grep('Skin|Blood|Adipose|Breast',colnames(human_mianshu3),value = T,invert = T)]
data_homo <- data_homo[rowSums(data_homo)>1,]
pheno <- data.frame(sample=colnames(data_homo))
edata <- data_homo
rownames(pheno) <- pheno$sample
pheno$batch <- ifelse(grepl('1',pheno$sample),'SH','HS')
#
batch = pheno$batch
modcombat = model.matrix(~1, data=pheno)
#
combat_edata1 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE,  prior.plots=FALSE)
dat <- as.data.frame(combat_edata1)
dat1 <- as.data.frame(scale(dat))

color1 <- colorRampPalette( rev(brewer.pal(9, "Blues")))(255)  
color2 = colorRampPalette(c('#58A4C3', "white", "firebrick3"))(100)

hclust and Heatmap
sampleDists <- dist(t(dat))   #dist
sampleDistMatrix <- as.matrix(sampleDists)  
pheatmap(sampleDistMatrix,
         fontsize=8,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         angle_col=90,
         col=color1)


sampleDists <- dist(t(dat1))   
sampleDistMatrix <- as.matrix(sampleDists)  
pheatmap(sampleDistMatrix,
         fontsize=8,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         angle_col=90,
         col=color1)

##kmeans k
fviz_nbclust(dat1,kmeans,method = "wss",k.max=50)
set.seed(123)  #
cluster1 <-pheatmap(dat1,scale = 'row',kmeans_k = 13,cluster_cols = F,color = color2)
cluster1_mean <- as.data.frame(cluster1$kmeans$centers)
cluster1_mean$size <- as.numeric(cluster1[["kmeans"]][["size"]])
cluster1_mean$size1 <- paste0('cluster',' ',rownames(cluster1_mean),':',cluster1_mean$size)
rownames(cluster1_mean) <- cluster1_mean$size1
cluster1_mean <- cluster1_mean[,1:26]

pheatmap(cluster1_mean,cluster_rows = F,cluster_cols = F,color = color2)
pheatmap(dat1,cluster_rows = T,cluster_cols = F,color = color2,scale = 'row')

#
kmeans_cluster <- data.frame(names(cluster1[["kmeans"]][["cluster"]]),cluster1[["kmeans"]][["cluster"]])
colnames(kmeans_cluster) <- c('gene_name','cluster')
lung_list <- list(cluster4=subset(kmeans_cluster,cluster==4)$gene_name,
                  cluster12=subset(kmeans_cluster,cluster==12)$gene_name)

#metascape format
lung_data_metascape <- as.data.frame(sapply(lung_list, '[', seq(max(sapply(lung_list, length)))))
write.csv(lung_data_metascape,'lung_data_metascape.csv',quote = F,row.names = F,na = '')

