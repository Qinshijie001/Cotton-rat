setwd('./')
library(Seurat)
library(tidyverse)
library(patchwork)
library(cowplot)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggrepel)
library(ggpubr)
library(ggbreak)
library(readxl)
library(SoupX)

set.seed(123)
my36colors <-c( '#53A85F', '#E5D2DD', '#AA9A59', '#E95C59', '#D6E7A3', '#57C3F3',
                '#F4A460', '#F3B1A0', '#9ACD32', '#AB3282', '#DDA0DD', '#91D0BE', 
                '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA',
                '#58A4C3', '#E4C755', '#F7F398', '#F1BB72', '#E63863', '#E39A35',
                '#C1E6F3', '#6778AE', '#BD956A', '#B53E2B', '#712820', '#DCC1DD',
                '#EE6A50', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')
my36colors1 <- paste0(my36colors,'50')
my36colors2 <- paste0(my36colors,'80')
group_color <- c('#53A85F','#E5D2DD','#E95C59', '#F3B1A0')

blast_name <- as.data.frame(read_excel('../blast_name.xlsx'))
rownames(blast_name) <- blast_name$gene_id

dir1=c("matrix")
list.dirs(dir1,recursive = F)
nasal_dir <- grep('Nasal',list.dirs(dir1,recursive = F),value = T)
nasal_dir

dir2=c("raw_feature_matrix")
list.dirs(dir2,recursive = F)
nasal_dir1 <- grep('Nasal',list.dirs(dir2,recursive = F),value = T)
nasal_dir1

#Adjustment of viral gene expression
meta_all <- data.frame()
for (i in c(5:16)){        ##
  toc <- Read10X(data.dir = paste0(nasal_dir[i],'/filtered_feature_bc_matrix'))
  tod <- Read10X(data.dir = paste0(nasal_dir1[i],'/raw_feature_bc_matrix'))
  
  rownames(toc) <- gsub('-|_','.',rownames(toc))
  rownames(tod) <- gsub('-|_','.',rownames(tod))
  
  tod <- tod[rownames(toc),]
  
  all <- toc
  project_name <- str_split(gsub('.+/','',nasal_dir[i]),'_',n = 2)[[1]][2]
  all <- CreateSeuratObject(all,project=project_name)
  all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
  all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
  all <- ScaleData(all, features = rownames(all))
  all <- RunPCA(all, features = VariableFeatures(all), npcs = 50, verbose = F)
  all <- FindNeighbors(all, dims = 1:40)
  all <- FindClusters(all, resolution = 0.8)
  all <- RunUMAP(all, dims = 1:40)
  
  #
  matx <- all@meta.data
  
  sc = SoupChannel(tod, toc)
  sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
  #Specify the genes to be adjusted
  rsv_gene <- tail(rownames(all),10)
  ute <- estimateNonExpressingCells(sc = sc,nonExpressedGeneList = list(rsv=rsv_gene))
  sc <- calculateContaminationFraction(sc,nonExpressedGeneList = list(rsv=rsv_gene),useToEst = ute,forceAccept = TRUE)
  #
  out = adjustCounts(sc)
  
  rsv_exp <- as.data.frame(t(as.matrix(tail(all@assays$RNA@counts,10))))
  rsv_exp$RSV_load <- rowSums(rsv_exp)
  all@meta.data$RSV_load <- rsv_exp$RSV_load
  all@meta.data$RSV_load_log <- log2(all@meta.data$RSV_load+1)
  
  rsv_exp1 <- as.data.frame(t(as.matrix(tail(out,10))))
  rsv_exp1$RSV_load1 <- rowSums(rsv_exp1)
  all@meta.data$RSV_load_soupx <- rsv_exp1$RSV_load1
  all@meta.data$RSV_load_log_soupx <- log2(all@meta.data$RSV_load_soupx+1)
  meta1 <- cbind(all@meta.data,rsv_exp1)
  rownames(meta1) <- paste(project_name,rownames(meta1),sep = '_')
  meta_all <- rbind(meta_all,meta1)
  
  table(meta1$RSV_load>0)
  table(meta1$RSV_load_soupx>0)
}
write.csv(meta_all,'meta_all_soupx.csv',quote = F)

#
scRNA_merge2 <- readRDS('../scRNA_merge2.RDS')
meta2 <- scRNA_merge2@meta.data
meta2$barcode <- rownames(meta2)
meta2$rsv_soupx <- 0
meta_all$barcode <- rownames(meta_all)
meta_all1 <- subset(meta_all,select = c('barcode','RSV_load_soupx','RSV_load_log_soupx'))
meta3 <- merge(meta2,meta_all1,by='barcode',all.x = T)
meta3$rsv_load_log_soupx <- meta3$rsv_soupx+meta3$RSV_load_log_soupx
meta3[is.na(meta3)] <- 0
rownames(meta3) <- meta3$barcode
meta3 <- meta3[rownames(meta2),]

scRNA_merge1 <- readRDS('../scRNA_list.RDS')
scRNA_merge1 <- merge(scRNA_merge1[[1]], y=scRNA_merge1[2:length(scRNA_merge1)],add.cell.id = names(scRNA_merge1)) 
dim(scRNA_merge1)
rsv_exp_all <- as.data.frame(t(as.matrix(tail(scRNA_merge1@assays$RNA@counts,10))))
rsv_exp_all <- rsv_exp_all[rownames(meta3),]
scRNA_merge2@meta.data <- cbind(scRNA_merge2@meta.data,rsv_exp_all)
rsv_exp_two <- rsv_exp_all
rsv_exp_two[rsv_exp_two > 0 ] = 1
rsv_exp_two <- rsv_exp_two[rowSums(rsv_exp_two)>1,]

meta4 <- subset(meta3,barcode %in% rownames(rsv_exp_two) & rsv_load_log_soupx>0) #RSV positivity is defined as an adjusted viral load greater than 0 and at least two transcripts detected
#
meta3$rsv_load_log_soupx[!(meta3$barcode %in% rownames(meta4))] = 0
scRNA_merge2@meta.data$rsv_load_log_soupx <- meta3$rsv_load_log_soupx

##
DimPlot(scRNA_merge2,cells.highlight = rownames(meta4),sizes.highlight = 0.1,raster = F,pt.size = 0.05)
ggsave('RSV_umap_all.pdf',width = 5,height = 3.3)

DimPlot(scRNA_merge2,cells.highlight = rownames(meta4),split.by = 'group_state',sizes.highlight = 0.1,raster = F,pt.size = 0.05)
ggsave('RSV_umap_group.pdf',width = 11,height = 3.2)

FeaturePlot(scRNA_merge2,features = 'rsv_load_log_soupx',cols = c('lightgray','red'))
ggsave('RSV_umap_viral_load.pdf',width = 4,height = 3.3)

scRNA_merge2_rsv <- subset(scRNA_merge2,subset = rsv_load_log_soupx>0)
saveRDS(scRNA_merge2_rsv,'scRNA_merge2_rsv.RDS')

RidgePlot(scRNA_merge2_rsv,features = 'rsv_load_log_soupx',cols = my36colors)
ggsave('RSV_umap_viral_load_pos_cell1.pdf',width = 6,height = 9)

FeaturePlot(scRNA_merge2_rsv,features = 'rsv_load_log_soupx',cols = c('lightgray','red'),
            label = T,label.size = 2.5)
ggsave('RSV_umap_viral_load1.pdf',width = 4,height = 3.3)

meta_rsv <- scRNA_merge2_rsv@meta.data
data1 <- data.frame(orig.ident=names(table(meta_rsv$orig.ident)),number=table(meta_rsv$orig.ident))
ggbarplot(data1,x = 'orig.ident',y = 'number.Freq',fill = '#CCCC33',label = data1$number.Freq,lab.vjust=0.5,lab.hjust = 0.1)+coord_flip()
ggsave('rsv pos cell of sample.pdf',width = 2.6,height = 3)

data2 <- data.frame(cell_name=names(table(meta_rsv$cell_name)),number=table(meta_rsv$cell_name),
        percent=round(table(meta_rsv$cell_name)/table(subset(meta3,group_state != 'Control')$cell_name)*100,2))
ggbarplot(data2,x = 'cell_name',y = 'number.Freq',fill = 'cell_name',label = data2$number.Freq,lab.vjust=0.5,palette = my36colors)+
  theme(axis.text.x = element_text(angle = 45,hjust=1))
ggsave('rsv pos cell of celltype.pdf',width = 10,height = 5)

ggbarplot(data2,x = 'cell_name',y = 'percent.Freq',fill = 'cell_name',label = data2$percent.Freq,lab.vjust=0.5,palette = my36colors)+
  theme(axis.text.x = element_text(angle = 45,hjust=1))
ggsave('rsv pos percent of celltype.pdf',width = 10,height = 5)

#
scRNA_merge2_rsv@meta.data$RSV.NS1.log <- log2(scRNA_merge2_rsv@meta.data$RSV.NS1+1)
scRNA_merge2_rsv@meta.data$RSV.NS2.log <- log2(scRNA_merge2_rsv@meta.data$RSV.NS2+1)
scRNA_merge2_rsv@meta.data$RSV.N.log <- log2(scRNA_merge2_rsv@meta.data$RSV.N+1)
scRNA_merge2_rsv@meta.data$RSV.P.log <- log2(scRNA_merge2_rsv@meta.data$RSV.P+1)
scRNA_merge2_rsv@meta.data$RSV.M.log <- log2(scRNA_merge2_rsv@meta.data$RSV.M+1)
scRNA_merge2_rsv@meta.data$RSV.SH.log <- log2(scRNA_merge2_rsv@meta.data$RSV.SH+1)
scRNA_merge2_rsv@meta.data$RSV.G.log <- log2(scRNA_merge2_rsv@meta.data$RSV.G+1)
scRNA_merge2_rsv@meta.data$RSV.F.log <- log2(scRNA_merge2_rsv@meta.data$RSV.F+1)
scRNA_merge2_rsv@meta.data$RSV.M2.log <- log2(scRNA_merge2_rsv@meta.data$RSV.M2+1)
scRNA_merge2_rsv@meta.data$RSV.L.log <- log2(scRNA_merge2_rsv@meta.data$RSV.L+1)

FeaturePlot(scRNA_merge2_rsv,features = tail(colnames(scRNA_merge2_rsv@meta.data),10),cols = c('lightgray','red'),ncol = 2)
ggsave('umap_rsv_gene.pdf',width = 6.5,height = 12)

DotPlot(scRNA_merge2_rsv,features =tail(colnames(scRNA_merge2_rsv@meta.data),10), group.by = 'cell_name',cols = c('#336699', '#E95C59'))+coord_flip()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave('dot_rsv_gene.pdf',width = 11,height = 4)

#candidate receptor
FeaturePlot(scRNA_merge2_rsv,features = c('Chr07G000072.CX3CR1','Chr07G001098.ICAM1','Chr01G000461.EGFR','Chr09G000267.HSPG2',
                                          'Chr09G000914.TLR4','Chr21G000176.NCL','Chr15G000361.INSR','Chr10G000614.INSRR'),cols = c('lightgray','red'),ncol = 2)
ggsave('umap_recepter_rsv.pdf',width = 7.2,height = 11)
