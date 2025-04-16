setwd('./')
library(Seurat)
library(tidyverse)
library(patchwork)
library(cowplot)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(readxl)

my36colors <-c( '#53A85F', '#E5D2DD', '#AA9A59', '#E95C59', '#D6E7A3', '#57C3F3',
                '#F4A460', '#F3B1A0', '#9ACD32', '#AB3282', '#DDA0DD', '#91D0BE', 
                '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA',
                '#58A4C3', '#E4C755', '#F7F398', '#F1BB72', '#E63863', '#E39A35',
                '#C1E6F3', '#6778AE', '#BD956A', '#B53E2B', '#712820', '#DCC1DD',
                '#EE6A50', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')
my36colors1 <- paste0(my36colors,'50')
my36colors2 <- paste0(my36colors,'80')
group_color <- c('#53A85F','#E5D2DD','#E95C59', '#F3B1A0')

scRNA_merge2 <- readRDS('../scRNA_merge2.RDS')
level_cell <- c('Basal','Club','Ciliated','Goblet',"Secretory1","Secretory2",'Sus','OEC','OSN','iOSN','iNeuron',
                'Astrocytes','Glial',"EC",'LEC','SMC','Myocyte','Fib','Myofib','CD4T','Treg','CD8T','NK','B',
                'Plasma','DC','Mono_Mac', 'pDC','Neu','Pro.T','Pro.B','Pro.M','RBC')

#Pseudo differential expression analysis
library(Libra)
DefaultAssay(scRNA_merge2) <- 'RNA'
DEG_data_frame_pseudo <- data.frame()

##
group_vs <- c('Control','2dpi')
temp_seurat <- subset(scRNA_merge2,subset = group_state %in% group_vs)
expr <- temp_seurat@assays$RNA@counts
meta <- temp_seurat@meta.data[,c('orig.ident','cell_name','group_state')]
colnames(meta) <- c('replicate','cell_type','label')
meta$label <- as.character(meta$label)
matrices = to_pseudobulk(expr, meta = meta,min_cells = 20,min_features = 0) 
DE = run_de(expr, meta = meta,replicate_col = "replicate",cell_type_col = "cell_type",label_col = "label",min_cells=20,min_features=0) 

#Filter low expression genes
for (i in names(matrices)){
  temp_matrix <- matrices[[i]]
  temp_matrix$gene_sum <- rowSums(temp_matrix)
  per_num <- table(temp_seurat$cell_name)[[i]]*0.1
  temp_matrix <- subset(temp_matrix,gene_sum >=per_num)
  temp_FC <- subset(DE,cell_type == i & gene %in% rownames(temp_matrix))
  temp_FC$vs <- paste(group_vs[1],group_vs[2],sep = 'vs')
  DEG_data_frame_pseudo <- rbind(DEG_data_frame_pseudo,temp_FC)
}

##
group_vs <- c('Control','5dpi')
temp_seurat <- subset(scRNA_merge2,subset = group_state %in% group_vs)
expr <- temp_seurat@assays$RNA@counts
meta <- temp_seurat@meta.data[,c('orig.ident','cell_name','group_state')]
colnames(meta) <- c('replicate','cell_type','label')
meta$label <- as.character(meta$label)
matrices = to_pseudobulk(expr, meta = meta,min_cells = 20,min_features = 0) 
DE = run_de(expr, meta = meta,replicate_col = "replicate",cell_type_col = "cell_type",label_col = "label",min_cells=20,min_features=0) 

#Filter low expression genes
for (i in names(matrices)){
  temp_matrix <- matrices[[i]]
  temp_matrix$gene_sum <- rowSums(temp_matrix)
  per_num <- table(temp_seurat$cell_name)[[i]]*0.1
  temp_matrix <- subset(temp_matrix,gene_sum >=per_num)
  temp_FC <- subset(DE,cell_type == i & gene %in% rownames(temp_matrix))
  temp_FC$vs <- paste(group_vs[1],group_vs[2],sep = 'vs')
  DEG_data_frame_pseudo <- rbind(DEG_data_frame_pseudo,temp_FC)
}

##
group_vs <- c('Control','8dpi')
temp_seurat <- subset(scRNA_merge2,subset = group_state %in% group_vs)
expr <- temp_seurat@assays$RNA@counts
meta <- temp_seurat@meta.data[,c('orig.ident','cell_name','group_state')]
colnames(meta) <- c('replicate','cell_type','label')
meta$label <- as.character(meta$label)
matrices = to_pseudobulk(expr, meta = meta,min_cells = 20,min_features = 0) 
DE = run_de(expr, meta = meta,replicate_col = "replicate",cell_type_col = "cell_type",label_col = "label",min_cells=20,min_features=0) 

#Filter low expression genes
for (i in names(matrices)){
  temp_matrix <- matrices[[i]]
  temp_matrix$gene_sum <- rowSums(temp_matrix)
  per_num <- table(temp_seurat$cell_name)[[i]]*0.1
  temp_matrix <- subset(temp_matrix,gene_sum >=per_num)
  temp_FC <- subset(DE,cell_type == i & gene %in% rownames(temp_matrix))
  temp_FC$vs <- paste(group_vs[1],group_vs[2],sep = 'vs')
  DEG_data_frame_pseudo <- rbind(DEG_data_frame_pseudo,temp_FC)
}

#
DEG_data_frame_pseudo$gene1 <- sub('.','_',DEG_data_frame_pseudo$gene,fixed = T)
DEG_data_frame_pseudo$gene1 <- gsub('.+_','',DEG_data_frame_pseudo$gene1)

DEG_data_frame_pseudo1 <- subset(DEG_data_frame_pseudo,p_val_adj < 0.05 & abs(avg_logFC) >= 0.5849625)
DEG_data_frame_pseudo1$cell_type <- factor(DEG_data_frame_pseudo1$cell_type,levels = level_cell)
write.csv(DEG_data_frame_pseudo1,file = 'DEG_data_frame_pseudo1.csv')

#
ggplot(data= DEG_data_frame_pseudo1, aes(cell_type,fill=vs))+geom_bar(position="dodge")+theme_classic()+
  scale_fill_manual(values = my36colors)+theme(axis.text.x = element_text(angle = 45,hjust=1))
ggsave('DEG_data_frame_pseudo_count1.pdf',width = 9,height = 3)
DEG_data_frame_pseudo1$up_down <- ifelse(DEG_data_frame_pseudo1$avg_logFC>0,'up','down')


##Functional enrichment analysis
library(clusterProfiler)
library(enrichplot)

GO <- read.table('GO_enrich_blast_name.txt',header = T,sep = '\t',quote = '')
GO <- subset(GO,V2 %in% names(table(GO$V2)[table(GO$V2)>5]))
BP <- subset(GO,V4=='biological_process',select = c('V3','blast_name1'))
colnames(BP) <- c("TERM","GENE")

KEGG <- read.table('KEGG_enrich_blast_name.txt',header = T,sep = '\t',quote = '')
KEGG <- subset(KEGG,V2 %in% names(table(KEGG$V2)[table(KEGG$V2)>5])) 
KEGG$pathway <-  paste(KEGG$V2,KEGG$V3,sep = ' ')
KEGG <- KEGG[,c('pathway','blast_name1')]
colnames(KEGG) <- c("TERM","GENE")

#
DEG <- as.matrix(subset(DEG_data_frame_pseudo1,vs=='Controlvs5dpi' & avg_logFC>0 & cell_type != 'RBC'))
DEG.list <- split(DEG[,2], DEG[,1])
DEG.go <- data.frame()
for (i in names(DEG.list)){
  go_temp <- enricher(DEG.list[[i]],TERM2GENE = BP,pAdjustMethod = "BH",pvalueCutoff=0.05,qvalueCutoff = 0.1)
  go_temp1 <- as.data.frame(go_temp)
  go_temp1$cell_name <- rep(i,nrow(go_temp1))
  DEG.go <- rbind(DEG.go,go_temp1)
}

DEG.go1 <- DEG.go %>% group_by(cell_name) %>% top_n(n=-5,wt=pvalue)
DEG.go1$Description <- factor(DEG.go1$Description,levels = unique(DEG.go1$Description))
DEG.go1$cell_name <- factor(DEG.go1$cell_name,levels = level_cell)
ggplot(DEG.go1,aes(cell_name,Description))+geom_point(aes(size=Count,color=pvalue))+theme_bw()+
  scale_colour_gradient(low="#CD3700",high="gray")+theme(axis.title= element_text(size=7, color="black"),
  axis.text.x=element_text(size=7,color = "black",angle = 90,hjust = 1,vjust = 0.4), 
  axis.text.y=element_text(size=7,color = "black"))+scale_size(range=c(1, 3))

#
DEG <- as.matrix(subset(DEG_data_frame_pseudo1,vs=='Controlvs5dpi' & avg_logFC<0 & cell_type != 'RBC'))
DEG.list <- split(DEG[,2], DEG[,1])
DEG.go <- data.frame()
for (i in names(DEG.list)){
  go_temp <- enricher(DEG.list[[i]],TERM2GENE = BP,pAdjustMethod = "BH",pvalueCutoff=0.1,qvalueCutoff = 1)
  go_temp1 <- as.data.frame(go_temp)
  go_temp1$cell_name <- rep(i,nrow(go_temp1))
  DEG.go <- rbind(DEG.go,go_temp1)
}

DEG.go1 <- DEG.go %>% group_by(cell_name) %>% top_n(n=-5,wt=pvalue)
DEG.go1$Description <- factor(DEG.go1$Description,levels = unique(DEG.go1$Description))
DEG.go1$cell_name <- factor(DEG.go1$cell_name,levels = level_cell)

ggplot(DEG.go1,aes(cell_name,Description))+geom_point(aes(size=Count,color=pvalue))+theme_bw()+
  scale_colour_gradient(low="blue",high="gray")+theme(axis.title= element_text(size=7, color="black"),
  axis.text.x=element_text(size=7,color = "black",angle = 90,hjust = 1,vjust = 0.4), 
  axis.text.y=element_text(size=7,color = "black"))+scale_size(range=c(1, 3))


#DEG frequency
library(ggalluvial)

vs_group <- 'Controlvs5dpi'
deg_up <- subset(DEG_data_frame_pseudo1,vs==vs_group & avg_logFC>0 & cell_type != 'RBC')
deg_up_freq <- data.frame(gene=names(table(deg_up$gene)),counts=as.numeric(table(deg_up$gene)))
deg_up_freq1 <- merge(deg_up,deg_up_freq,by = 'gene')
deg_up_freq1 <- deg_up_freq1[order(deg_up_freq1$counts,decreasing = T),]
write.table(deg_up_freq1,paste0(vs_group,'.up.freq.txt'),row.names = F,quote = F)

#Sankey diagram
deg_up_freq2 <- subset(deg_up_freq1,counts %in% unique(deg_up_freq1$counts)[1:5]) 
df <- to_lodes_form(deg_up_freq2[,c('cell_type','gene1')],axes = 1:2,id = "value")
df$stratum <- factor(df$stratum,levels = c(level_cell,unique(deg_up_freq2$gene1)))
col <- c(my36colors[1:length(unique(deg_up_freq2$cell_type))],rep('#F08080',1000))
ggplot(df, aes(x = x, fill=stratum, label=stratum,stratum = stratum, alluvium  = value))+
geom_flow(width = 0.3,curve_type = "sine",alpha = 0.5,color = 'white', size = 0.1)+
  geom_stratum(width = 0.6)+geom_text(stat = 'stratum', size = 2.5, color = 'black')+
  scale_fill_manual(values = col)+theme_void()+theme(legend.position = 'none')
ggsave(paste0(vs_group,'.up.freq.pdf'),width = 2.8,height = 6)

vs_group <- 'Controlvs5dpi'
deg_down <- subset(DEG_data_frame_pseudo1,vs==vs_group & avg_logFC<0 & cell_type != 'RBC')
deg_down_freq <- data.frame(gene=names(table(deg_down$gene)),counts=as.numeric(table(deg_down$gene)))
deg_down_freq1 <- merge(deg_down,deg_down_freq,by = 'gene')
deg_down_freq1 <- deg_down_freq1[order(deg_down_freq1$counts,decreasing = T),]
write.table(deg_down_freq1,paste0(vs_group,'.down.freq.txt'),row.names = F,quote = F)

#
deg_down_freq2 <- subset(deg_down_freq1,counts %in% unique(deg_down_freq1$counts)[1:5]) 
df <- to_lodes_form(deg_down_freq2[,c('cell_type','gene1')],axes = 1:2,id = "value")
df$stratum <- factor(df$stratum,levels = c(level_cell,unique(deg_down_freq2$gene1)))
col <- c(my36colors[1:length(unique(deg_down_freq2$cell_type))],rep('#57C3F3',1000))
ggplot(df, aes(x = x, fill=stratum, label=stratum,stratum = stratum, alluvium  = value))+
geom_flow(width = 0.3,curve_type = "sine",alpha = 0.5,color = 'white', size = 0.1)+
  geom_stratum(width = 0.6)+geom_text(stat = 'stratum', size = 2.5, color = 'black')+
  scale_fill_manual(values = col)+theme_void()+theme(legend.position = 'none')
ggsave(paste0(vs_group,'.down.freq.pdf'),width = 2.8,height = 6)

#Differentially expressed cytokines
cytokine <- read.table('../cytokine.txt',header = T,sep = '\t')
cytokine_deg <- subset(DEG_data_frame_pseudo1,gene %in% cytokine$gene_name & cell_type != 'RBC')
cytokine_deg$cell_type <- factor(cytokine_deg$cell_type,levels = level_cell)
ggplot(cytokine_deg,aes(cell_type,gene1))+geom_point(aes(color=up_down,shape=vs,size=abs(cytokine_deg$avg_logFC)),alpha=0.5)+
  scale_color_manual(values = c('blue', 'red'))+scale_size(range = c(2,4))+
  theme_bw(base_line_size = 0.2)+theme(axis.text.y = element_text(color='black',size = 8))+
  theme(axis.text.x = element_text(color='black',size = 8,angle = 45,hjust = 0.7,vjust = 1))
ggsave('cytokine_deg.pdf',width = 8,height = 11)

##Differentially expressed TFs
TF_list <- read.table('TF.txt',sep = '\t',header = T)
TF_list$id <- gsub('.1','',TF_list$TF,fixed = T)
TF_list <- subset(TF_list,Full.sequence.e.value < 0.00001 & Best.domain.e.value < 0.00001)
DEG_data_frame_pseudo1$id <- gsub('.','_',DEG_data_frame_pseudo1$gene,fixed = T)
DEG_data_frame_pseudo1$id <- gsub('_.+','',DEG_data_frame_pseudo1$id)
TF_list_deg <- merge(TF_list,DEG_data_frame_pseudo1,by = 'id')
TF_list_deg <- subset(TF_list_deg,cell_type != 'RBC')
write.table(TF_list_deg,'TF_list_deg.txt',quote = F,row.names = F,sep = '\t')

