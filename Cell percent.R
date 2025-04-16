setwd('./')
library(Seurat)
library(tidyverse)
library(patchwork)
library(cowplot)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggbreak)
library(readxl)

my36colors <-c( '#53A85F', '#E5D2DD', '#AA9A59', '#E95C59', '#D6E7A3', '#57C3F3',
                '#F4A460', '#F3B1A0', '#9ACD32', '#AB3282', '#DDA0DD', '#91D0BE', 
                '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA',
                '#58A4C3', '#E4C755', '#F7F398', '#F1BB72', '#E63863', '#E39A35',
                '#C1E6F3', '#6778AE', '#BD956A', '#B53E2B', '#712820', '#DCC1DD',
                '#EE6A50', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')

#
my36colors1 <- paste0(my36colors,'50')
group_color <- c('#53A85F','#E5D2DD','#E95C59', '#F3B1A0')

scRNA_merge2 <- readRDS('../scRNA_merge2.RDS')
meta2 <- scRNA_merge2@meta

#Cell percentage
library(RColorBrewer)
ggplot(meta2, aes(group_state)) + labs(y='Cell percent',x='Group state')+
  geom_bar(aes(fill=cell_name), position="fill")+
  scale_fill_manual(values = my36colors)+theme_bw()

ggplot(meta2, aes(orig.ident)) + labs(y='Cell percent',x='Sample')+
  geom_bar(aes(fill=cell_name), position="fill")+
  scale_fill_manual(values = my36colors)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust=0.5))

df_percent=data.frame()
for (i in as.character(unique(meta2$cell_name))){
  df_test <- meta2 %>% group_by(orig.ident) %>% summarise(count_rows = n(),prop = sum(cell_name==i)/count_rows)
  df_test$cell_name <- rep(i,nrow(df_test))
  df_percent <- rbind(df_percent,df_test)
}

df_percent$group_state <- gsub('_.+','',df_percent$orig.ident)
df_percent$cell_name <- factor(df_percent$cell_name,levels =level_cell)

my_compare <- list(c('C','T2'),c('C','T5'),c('C','T8')) 

ggboxplot(df_percent, x = "group_state", y = 'prop',fill='group_state', palette = group_color)+
  stat_compare_means(aes(group=group_state),method = 'kruskal.test',label = "p.format",size=3)+
  theme(axis.text.x = element_text(angle = 0,hjust = 1))+
  facet_wrap(~cell_name,scales="free",nrow = 2)+theme_bw()

#Milo
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)

scRNA_merge2@meta.data$group_state1 <- as.character(scRNA_merge2@meta.data$group_state)
da_results_data <- data.frame()

target_group <- c('Control','2dpi')
test1 <- subset(scRNA_merge2,group_state1==target_group[1] | group_state1==target_group[2])
test1@meta.data$group_state1 <- factor(test1@meta.data$group_state1,levels = target_group)
scRNA_merge1_sce <- as.SingleCellExperiment(test1)

scRNA_merge1_sce <- runPCA(scRNA_merge1_sce, ncomponents=40)
scRNA_merge1_sce <- runUMAP(scRNA_merge1_sce)
plotUMAP(scRNA_merge1_sce,colour_by="cell_name")

scRNA_merge1_milo <- Milo(scRNA_merge1_sce)
plotUMAP(scRNA_merge1_milo)

scRNA_merge1_milo <- buildGraph(scRNA_merge1_milo, k = 10, d = 40)
scRNA_merge1_milo <- makeNhoods(scRNA_merge1_milo, prop = 0.1, k = 10, d=40, refined = TRUE,refinement_scheme="graph")
plotNhoodSizeHist(scRNA_merge1_milo)

scRNA_merge1_milo@colData$Sample <- scRNA_merge1_milo@colData$orig.ident
scRNA_merge1_milo@colData$Condition <- scRNA_merge1_milo@colData$group_state1

scRNA_merge1_milo <- countCells(scRNA_merge1_milo, meta.data = data.frame(colData(scRNA_merge1_milo)), samples="Sample")
head(nhoodCounts(scRNA_merge1_milo))

scRNA_merge1_design <- data.frame(colData(scRNA_merge1_milo))[,c("Sample", "Condition")]
scRNA_merge1_design <- distinct(scRNA_merge1_design)
rownames(scRNA_merge1_design) <- scRNA_merge1_design$Sample
## Reorder rownames to match columns of nhoodCounts(milo)
scRNA_merge1_design <- scRNA_merge1_design[colnames(nhoodCounts(scRNA_merge1_milo)), ,drop=FALSE]

scRNA_merge1_milo <- calcNhoodDistance(scRNA_merge1_milo, d=40)
rownames(scRNA_merge1_design) <- scRNA_merge1_design$Sample
da_results <- testNhoods(scRNA_merge1_milo, design = ~ Condition, design.df = scRNA_merge1_design)
da_results %>% arrange(- SpatialFDR) %>% head() 

scRNA_merge1_milo <- buildNhoodGraph(scRNA_merge1_milo)
plotUMAP(scRNA_merge1_milo,colour_by="cell_name") + plotNhoodGraphDA(scRNA_merge1_milo, da_results, alpha=0.05) +
  plot_layout(guides="collect")
ggsave(paste0('umap_nhood',target_group[1],'vs',target_group[2],'.pdf'),width=8.3,height=4)

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

da_results <- annotateNhoods(scRNA_merge1_milo, da_results, coldata_col = "cell_name")
ggplot(da_results, aes(cell_name_fraction)) + geom_histogram(bins=50)
da_results$cell_name <- ifelse(da_results$cell_name_fraction < 0.7, "Mixed", da_results$cell_name)
da_results$cell_name <- factor(da_results$cell_name,levels = c(intersect(level_cell,unique(da_results$cell_name)),"Mixed"))

plotDAbeeswarm(da_results, group.by = "cell_name")+geom_jitter(width = 0.1,size=0.01)+
  scale_color_gradientn(colours = c('blue','white','red'))+theme_bw()
ggsave(paste0('umap_nhood',target_group[1],'vs',target_group[2],'logFC.pdf'),width=2.5,height=5)

da_results$vs <- rep(paste0(target_group[1],'vs',target_group[2]),nrow(da_results))
write.csv(da_results,paste0('umap_nhood',target_group[1],'vs',target_group[2],'.csv'),quote = F)
da_results_data <- rbind(da_results_data,da_results)

ggplot(da_results,aes(x = cell_name,y=logFC,color=SpatialFDR))+geom_jitter(width = 0.1)+
  coord_flip()+scale_color_continuous(low='#E95C59', high='#53A85F')+theme_bw()

#
target_group <- c('2dpi','5dpi')
test1 <- subset(scRNA_merge2,group_state1==target_group[1] | group_state1==target_group[2])
test1@meta.data$group_state1 <- factor(test1@meta.data$group_state1,levels = target_group)
scRNA_merge1_sce <- as.SingleCellExperiment(test1)

scRNA_merge1_sce <- runPCA(scRNA_merge1_sce, ncomponents=40)
scRNA_merge1_sce <- runUMAP(scRNA_merge1_sce)
plotUMAP(scRNA_merge1_sce,colour_by="cell_name")

scRNA_merge1_milo <- Milo(scRNA_merge1_sce)
plotUMAP(scRNA_merge1_milo)

scRNA_merge1_milo <- buildGraph(scRNA_merge1_milo, k = 10, d = 40)
scRNA_merge1_milo <- makeNhoods(scRNA_merge1_milo, prop = 0.1, k = 10, d=40, refined = TRUE,refinement_scheme="graph")
plotNhoodSizeHist(scRNA_merge1_milo)

scRNA_merge1_milo@colData$Sample <- scRNA_merge1_milo@colData$orig.ident
scRNA_merge1_milo@colData$Condition <- scRNA_merge1_milo@colData$group_state1

scRNA_merge1_milo <- countCells(scRNA_merge1_milo, meta.data = data.frame(colData(scRNA_merge1_milo)), samples="Sample")
head(nhoodCounts(scRNA_merge1_milo))

scRNA_merge1_design <- data.frame(colData(scRNA_merge1_milo))[,c("Sample", "Condition")]
scRNA_merge1_design <- distinct(scRNA_merge1_design)
rownames(scRNA_merge1_design) <- scRNA_merge1_design$Sample
## Reorder rownames to match columns of nhoodCounts(milo)
scRNA_merge1_design <- scRNA_merge1_design[colnames(nhoodCounts(scRNA_merge1_milo)), ,drop=FALSE]

scRNA_merge1_milo <- calcNhoodDistance(scRNA_merge1_milo, d=40)
rownames(scRNA_merge1_design) <- scRNA_merge1_design$Sample
da_results <- testNhoods(scRNA_merge1_milo, design = ~ Condition, design.df = scRNA_merge1_design)
da_results %>% arrange(- SpatialFDR) %>% head() 

scRNA_merge1_milo <- buildNhoodGraph(scRNA_merge1_milo)
plotUMAP(scRNA_merge1_milo,colour_by="cell_name") + plotNhoodGraphDA(scRNA_merge1_milo, da_results, alpha=0.05) +
  plot_layout(guides="collect")
ggsave(paste0('umap_nhood',target_group[1],'vs',target_group[2],'.pdf'),width=8.3,height=4)

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

da_results <- annotateNhoods(scRNA_merge1_milo, da_results, coldata_col = "cell_name")
ggplot(da_results, aes(cell_name_fraction)) + geom_histogram(bins=50)
da_results$cell_name <- ifelse(da_results$cell_name_fraction < 0.7, "Mixed", da_results$cell_name)
da_results$cell_name <- factor(da_results$cell_name,levels = c(intersect(level_cell,unique(da_results$cell_name)),"Mixed"))

plotDAbeeswarm(da_results, group.by = "cell_name")+geom_jitter(width = 0.1,size=0.01)+
  scale_color_gradientn(colours = c('blue','white','red'))+theme_bw()
ggsave(paste0('umap_nhood',target_group[1],'vs',target_group[2],'logFC.pdf'),width=2.5,height=5)

da_results$vs <- rep(paste0(target_group[1],'vs',target_group[2]),nrow(da_results))
write.csv(da_results,paste0('umap_nhood',target_group[1],'vs',target_group[2],'.csv'),quote = F)
da_results_data <- rbind(da_results_data,da_results)

ggplot(da_results,aes(x = cell_name,y=logFC,color=SpatialFDR))+geom_jitter(width = 0.1)+
  coord_flip()+scale_color_continuous(low='#E95C59', high='#53A85F')+theme_bw()
  
#
target_group <- c('5dpi','8dpi')
test1 <- subset(scRNA_merge2,group_state1==target_group[1] | group_state1==target_group[2])
test1@meta.data$group_state1 <- factor(test1@meta.data$group_state1,levels = target_group)
scRNA_merge1_sce <- as.SingleCellExperiment(test1)

scRNA_merge1_sce <- runPCA(scRNA_merge1_sce, ncomponents=40)
scRNA_merge1_sce <- runUMAP(scRNA_merge1_sce)
plotUMAP(scRNA_merge1_sce,colour_by="cell_name")

scRNA_merge1_milo <- Milo(scRNA_merge1_sce)
plotUMAP(scRNA_merge1_milo)

scRNA_merge1_milo <- buildGraph(scRNA_merge1_milo, k = 10, d = 40)
scRNA_merge1_milo <- makeNhoods(scRNA_merge1_milo, prop = 0.1, k = 10, d=40, refined = TRUE,refinement_scheme="graph")
plotNhoodSizeHist(scRNA_merge1_milo)

scRNA_merge1_milo@colData$Sample <- scRNA_merge1_milo@colData$orig.ident
scRNA_merge1_milo@colData$Condition <- scRNA_merge1_milo@colData$group_state1

scRNA_merge1_milo <- countCells(scRNA_merge1_milo, meta.data = data.frame(colData(scRNA_merge1_milo)), samples="Sample")
head(nhoodCounts(scRNA_merge1_milo))

scRNA_merge1_design <- data.frame(colData(scRNA_merge1_milo))[,c("Sample", "Condition")]
scRNA_merge1_design <- distinct(scRNA_merge1_design)
rownames(scRNA_merge1_design) <- scRNA_merge1_design$Sample
## Reorder rownames to match columns of nhoodCounts(milo)
scRNA_merge1_design <- scRNA_merge1_design[colnames(nhoodCounts(scRNA_merge1_milo)), ,drop=FALSE]

scRNA_merge1_milo <- calcNhoodDistance(scRNA_merge1_milo, d=40)
rownames(scRNA_merge1_design) <- scRNA_merge1_design$Sample
da_results <- testNhoods(scRNA_merge1_milo, design = ~ Condition, design.df = scRNA_merge1_design)
da_results %>% arrange(- SpatialFDR) %>% head() 

scRNA_merge1_milo <- buildNhoodGraph(scRNA_merge1_milo)
plotUMAP(scRNA_merge1_milo,colour_by="cell_name") + plotNhoodGraphDA(scRNA_merge1_milo, da_results, alpha=0.05) +
  plot_layout(guides="collect")
ggsave(paste0('umap_nhood',target_group[1],'vs',target_group[2],'.pdf'),width=8.3,height=4)

ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

da_results <- annotateNhoods(scRNA_merge1_milo, da_results, coldata_col = "cell_name")
ggplot(da_results, aes(cell_name_fraction)) + geom_histogram(bins=50)
da_results$cell_name <- ifelse(da_results$cell_name_fraction < 0.7, "Mixed", da_results$cell_name)
da_results$cell_name <- factor(da_results$cell_name,levels = c(intersect(level_cell,unique(da_results$cell_name)),"Mixed"))

plotDAbeeswarm(da_results, group.by = "cell_name")+geom_jitter(width = 0.1,size=0.01)+
  scale_color_gradientn(colours = c('blue','white','red'))+theme_bw()
ggsave(paste0('umap_nhood',target_group[1],'vs',target_group[2],'logFC.pdf'),width=2.5,height=5)

da_results$vs <- rep(paste0(target_group[1],'vs',target_group[2]),nrow(da_results))
write.csv(da_results,paste0('umap_nhood',target_group[1],'vs',target_group[2],'.csv'),quote = F)
da_results_data <- rbind(da_results_data,da_results)

ggplot(da_results,aes(x = cell_name,y=logFC,color=SpatialFDR))+geom_jitter(width = 0.1)+
  coord_flip()+scale_color_continuous(low='#E95C59', high='#53A85F')+theme_bw()
