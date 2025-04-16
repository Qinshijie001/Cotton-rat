setwd('./')
library(Seurat)
library(tidyverse)
library(patchwork)
library(cowplot)
library(dplyr)
library(ggplot2)
library(DoubletFinder)
library(viridis)
library(ggrepel)
library(harmony)
library(ggpubr)
library(ggbreak)
library(readxl)

my36colors <-c( '#53A85F', '#E5D2DD', '#AA9A59', '#E95C59', '#D6E7A3', '#57C3F3',
                '#F4A460', '#F3B1A0', '#9ACD32', '#AB3282', '#DDA0DD', '#91D0BE', 
                '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA',
                '#58A4C3', '#E4C755', '#F7F398', '#F1BB72', '#E63863', '#E39A35',
                '#C1E6F3', '#6778AE', '#BD956A', '#B53E2B', '#712820', '#DCC1DD',
                '#EE6A50', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')

##
my36colors1 <- paste0(my36colors,'50')
group_color <- c('#53A85F','#E5D2DD','#E95C59', '#F3B1A0')

blast_name <- as.data.frame(read_excel('blast_name.xlsx'))
rownames(blast_name) <- blast_name$gene_id

dir1=c("../matrix")
list.dirs(dir1,recursive = F)
Nasal_dir <- grep('Nasal',list.dirs(dir1,recursive = F),value = T)
Nasal_dir

#Load data
scRNAlist <- list()
for(i in Nasal_dir){
  print(i)
  counts <- Read10X(data.dir = paste0(i,'/filtered_feature_bc_matrix'))
  new_name <- blast_name[rownames(counts),]
  rownames(counts) <- new_name$blast_name
  counts@Dimnames[[1]] <- new_name$blast_name
  project_name <- str_split(i,'_',n = 2)[[1]][2]
  scRNAlist[[project_name]] <- CreateSeuratObject(counts,project = project_name)
}

#Remove doublet
for (i in c(1:length(scRNAlist))){
  scRNAlist[[i]]$percent.mt <- PercentageFeatureSet(scRNAlist[[i]], pattern = "mito")
  scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = percent.mt <= 30)
  scRNAlist[[i]] <- NormalizeData(scRNAlist[[i]])
  scRNAlist[[i]] <- FindVariableFeatures(scRNAlist[[i]], selection.method = "vst", nfeatures = 2000)
  scRNAlist[[i]] <- ScaleData(scRNAlist[[i]]) 
  scRNAlist[[i]] <- RunPCA(scRNAlist[[i]]) 
  scRNAlist[[i]] <- RunUMAP(scRNAlist[[i]], dims = 1:30)
  scRNAlist[[i]] <- FindNeighbors(scRNAlist[[i]], reduction = "pca", dims = 1:30)
  scRNAlist[[i]] <- FindClusters(scRNAlist[[i]], resolution = 0.5)
  
  sweep.res.list <- paramSweep_v3(scRNAlist[[i]], PCs = 1:30, sct = F)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  homotypic.prop <- modelHomotypic(scRNAlist[[i]]$seurat_clusters)
  nExp_poi <- round(0.075*ncol(scRNAlist[[i]]))
  nExp_poi.adj<- round(nExp_poi*( 1-homotypic.prop))
  
  scRNAlist[[i]] <- doubletFinder_v3(scRNAlist[[i]], PCs = 1:30, pN = 0.25, pK = mpK, 
                                     nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
  
  scRNAlist[[i]]@meta.data$Doublts_score <- scRNAlist[[i]]@meta.data[,7]  
  scRNAlist[[i]]@meta.data$Doublts_anno <- scRNAlist[[i]]@meta.data[,8]  
  
}

#Merge data and filter
scRNA_merge <- merge(scRNAlist[[1]], y=scRNAlist[2:length(scRNAlist)],add.cell.id = names(scRNAlist)) 
scRNA_merge <- subset(scRNA_merge,subset = Doublts_anno == 'Singlet')
scRNA_merge@meta.data <- scRNA_merge@meta.data[,c('orig.ident','nCount_RNA','nFeature_RNA','percent.mt','Doublts_score')]
#
HB.genes <- HB.genes <- c("Chr11G001546.HBA","Chr11G001548.HBA","Chr11G001549.HBZ","Chr24G000283.HBB")
scRNA_merge[["percent.HB"]]<-PercentageFeatureSet(scRNA_merge, features=HB.genes) 
##
minGene=300
maxGene=4500
minRNA=500
maxRNA=15000
pctHB=15
scRNA_merge1 <- subset(scRNA_merge, subset = nFeature_RNA >= minGene & nFeature_RNA <= maxGene & percent.HB <= pctHB & nCount_RNA <= maxRNA & nCount_RNA >= minRNA)

VlnPlot(scRNA_merge1, group.by = "orig.ident",  
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB"), 
        cols =rainbow(16), 
        pt.size = 0, 
        ncol = 4) 
ggsave('quality.pdf',width = 20,height = 4)

scRNA_merge1@meta.data$group_state <- scRNA_merge1@meta.data$orig.ident
scRNA_merge1@meta.data$group_state[grepl('C',scRNA_merge1@meta.data$group_state)] <- 'Control'
scRNA_merge1@meta.data$group_state[grepl('T2',scRNA_merge1@meta.data$group_state)] <- '2dpi'
scRNA_merge1@meta.data$group_state[grepl('T5',scRNA_merge1@meta.data$group_state)] <- '5dpi'
scRNA_merge1@meta.data$group_state[grepl('T8',scRNA_merge1@meta.data$group_state)] <- '8dpi'
scRNA_merge1@meta.data$group_state <- factor(scRNA_merge1@meta.data$group_state,levels = c('Control','2dpi', '5dpi','8dpi'))

#Load virus load
rsv_exp <- as.data.frame(t(as.matrix(tail(scRNA_merge1@assays$RNA@counts,10))))
rsv_exp$RSV_load <- rowSums(rsv_exp)
scRNA_merge1@meta.data$RSV_load <- rsv_exp$RSV_load
scRNA_merge1@meta.data$RSV_load_log <- log2(scRNA_merge1@meta.data$RSV_load+1)
scRNA_merge1 <- scRNA_merge1[1:21657,] #

#seurat pipline
DefaultAssay(scRNA_merge1) <- "RNA"
scRNA_merge1 <- NormalizeData(object = scRNA_merge1, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA_merge1 <- FindVariableFeatures(object = scRNA_merge1, selection.method = "vst", nfeatures = 2000)
scRNA_merge1 <- ScaleData(scRNA_merge1, verbose = FALSE,vars.to.regress = c("percent.mt"))
scRNA_merge1 <- RunPCA(scRNA_merge1, pc.genes = scRNA_merge1@var.genes, npcs = 80, verbose = FALSE)
scRNA_merge1 <- ProjectDim(object = scRNA_merge1)
ElbowPlot(object = scRNA_merge1,ndims = 80)
DimPlot(scRNA_merge1, reduction = "pca", group.by="orig.ident")
DimHeatmap(scRNA_merge1, dims = 40:50, cells = 500, balanced = TRUE)

pc.num=1:50
scRNA_merge1 <- scRNA_merge1 %>%
  RunHarmony(c("orig.ident",'group_state'), plot_convergence = F)
scRNA_merge1 <- scRNA_merge1 %>% 
  RunUMAP(reduction = "harmony", dims = pc.num) %>% 
  FindNeighbors(reduction = "harmony", dims = pc.num) %>% 
  FindClusters(resolution = c(0.3,0.5,1,1.5,2,3,4,5)) %>% 
  identity()
  
#Calculate cluster specific genes
DefaultAssay(scRNA_merge1) <- 'RNA'
Idents(scRNA_merge1) <- 'RNA_snn_res.3'
sc.markers <- FindAllMarkers(scRNA_merge1, only.pos =TRUE,logfc.threshold =1,min.pct = 0.7)
top50 <- sc.markers %>% group_by(cluster) %>% top_n(n=50,wt=avg_log2FC)
write.csv(sc.markers,'cluster_3_DEG.csv',quote = F,row.names = F)

Idents(scRNA_merge1) <- 'RNA_snn_res.0.5'
sc.markers.1 <- FindAllMarkers(scRNA_merge1, only.pos =TRUE,logfc.threshold =1,min.pct = 0.7)
top50.1 <- sc.markers.1 %>% group_by(cluster) %>% top_n(n=50,wt=avg_log2FC)
write.csv(sc.markers.1,'cluster_0.5_DEG.csv',quote = F,row.names = F)

Idents(scRNA_merge1) <- 'RNA_snn_res.2'
scRNA_merge2 <- RenameIdents(scRNA_merge1, 
                             '0' = "CD4T",   
                             '1' = "Fib",     
                             '2' = "OSN",   
                             '3' = "NK",   
                             '4' = "Myofib",  #
                             '5' = "CD8T",
                             '6' = "OSN", 
                             '7' = "Fib", 
                             '8' = "iOSN",  
                             '9' = "OSN",    
                             '10' = "B",                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                             '11' = "Mono_Mac", 
                             '12' = "Mono_Mac",  
                             '13' = "Fib", 
                             '14' = "Fib", 
                             '15' = "EC", 
                             '16' = "DC", 
                             '17' = "Secretory1",
                             '18' = "Treg", 
                             '19' = "Secretory1", 
                             '20' = "Sus", 
                             '21' = "Pro.T",    
                             '22' = "CD4T", 
                             '23' = "Mono_Mac", 
                             '24' = "Ciliated", 
                             '25' = "OSN", 
                             '26' = "Secretory2", 
                             '27' = "Fib", 
                             '28' = "OEC", 
                             '29' = "Sus", 
                             '30' = "pDC", 
                             '31' = "Glial",       
                             '32' = "Fib",   
                             '33' = "Plasma", 
                             '34' = "iOSN", 
                             '35' = "Myocyte", 
                             '36' = "SMC", 
                             '37' = "SMC",  
                             '38' = "B",
                             '39' = "CD4T",
                             '40' = "Basal", 
                             '41' = "Club", 
                             '42' = "EC",
                             '43' = "iNeuron",  
                             '44' = "Neu", 
                             '45' = "LEC",  
                             '46' = "Sus", 
                             '47' = "Astrocytes", 
                             '48' = "Goblet", 
                             '49' = "OSN", 
                             '50' = "Secretory1",
                             '51' = "RBC",
                             '52' = "Goblet", 
                             '53' = "Fib",
                             '54' = "Secretory1", 
                             '55' = "Pro.B", 
                             '56' = "Pro.M", 
                             '57' = "Fib",
                             '58' = "OSN",
                             '59' = "OSN",
                             '60' = "NK",
                             '61' = "Secretory1")

scRNA_merge2@meta.data$cell_name <- scRNA_merge2@active.ident
level_cell <- c('Basal','Club','Ciliated','Goblet',"Secretory1","Secretory2",'Sus','OEC','OSN','iOSN','iNeuron',
                'Astrocytes','Glial',"EC",'LEC','SMC','Myocyte','Fib','Myofib','CD4T','Treg','CD8T','NK','B',
                'Plasma','DC','Mono_Mac', 'pDC','Neu','Pro.T','Pro.B','Pro.M','RBC')
scRNA_merge2@meta.data$cell_name <- factor(scRNA_merge2$cell_name,levels = level_cell)
Idents(scRNA_merge2) <- 'cell_name'
DimPlot(object = scRNA_merge2, reduction = 'umap',label = T,cols = my36colors)
saveRDS(scRNA_merge2,'scRNA_merge2.RDS')

DimPlot(object = scRNA_merge2, reduction = 'umap',label = T,cols = my36colors,label.size = 2)

DimPlot(object = scRNA_merge2, reduction = 'umap',label = F,cols = group_color,label.size = 2,group.by = 'group_state',split.by = 'group_state')

#Calculate the number of cells
meta2 <- scRNA_merge2@meta.data
a=data.frame(names(table(meta2$cell_name)),table(meta2$cell_name))
ggbarplot(a,x = 'Var1',y = 'Freq',fill = 'Var1',label = a$Freq,palette = my36colors,lab.vjust=0.5)+coord_flip()

#Cell type heatmap
DefaultAssay(scRNA_merge2) <- 'RNA'
Idents(scRNA_merge2) <- 'cell_name'
sc.markers.cell <- FindAllMarkers(scRNA_merge2, only.pos =TRUE,logfc.threshold =0.4,min.pct = 0.4)
top50.cell <- sc.markers.cell %>% group_by(cluster) %>% top_n(n=50,wt=avg_log2FC)
write.csv(sc.markers.cell,'sc.markers.cell.csv',quote = F,row.names = F)

library(pheatmap)
library(ggplotify)
color1 = colorRampPalette(c('#58A4C3', "#F8F8FF", "firebrick3"))(50)

scRNA_merge_exp <- as.data.frame(AverageExpression(scRNA_merge2,group.by = 'cell_name'))
deg_exp <- scRNA_merge_exp[top50.cell$gene,]
colnames(deg_exp) <- gsub('RNA.','',colnames(deg_exp),fixed = T)
plot1=pheatmap(deg_exp,scale = 'row',cluster_cols = F,cluster_rows = F,annotation_names_col = T)

plot1=pheatmap(deg_exp,scale = 'row',cluster_cols = F,cluster_rows = F,annotation_names_col = T,color = color1)

#Visual marker
markers <- c("Chr26G000173.KRT18",
             "Chr11G000430.KRT17" ,"Chr11G000434.KRT14",
             "Chr25G000004.LYPD2","Chr02G000821.SCGB1C1",
             "Chr10G000506.MLF1","Chr06G000934.TPPP3",
             "Chr07G000764.PLET1","Chr24G000565.KLK15",
             "Chr11G001001.VMO1","Chr02G000774.MUC5B",
             "Chr23G000424.REG3G","Chr04G000197.WFDC2",
             "Chr01G000547.SEC14L3","Chr18G000196.CYP2A5",
             "Chr14G000303.S100B",'Chr17G000287.DBI',
             "Chr11G001509.GNG13","Chr03G000081.CNGA2",
             'Chr02G001356.GAP43',
             'Chr18G000348.GNG8',
             "Chr10G000238.PIP",
             'Chr18G000279.APOE','Chr07G000750.CRYAB',
             "Chr11G000259.PECAM1",
             "Chr09G001129.CCL21A","Chr11G001639.FLT4",
             "Chr02G000270.ACTA2",
             "Chr04G000187.TNNC2",'Chr02G000979.MYL11',
             "Chr15G000082.DCN",
             "Chr21G000306.FN1","Chr12G000356.MGP",'Chr16G000223.COMP',
             "Chr07G000815.CD3D","Chr01G000287.CD4",
             "Chr06G000503.IL2RA","Chr07G000521.RORA",
             "Chr23G000379.CD8A", "Chr23G000378.CD8B",
             "Chr01G000375.KLRK1", "Chr24G000539.NKG7",
             "Chr18G000229.CD79A","Chr02G000445.MS4A1",
             'Chr13G000364.MZB1',"Chr01G000879.JCHAIN",
             "Chr14G000095.H2-EB1", "Chr14G000094.H2-AA",
             "Chr15G000360.LYZ","Chr13G000340.CD14",
             "Chr04G000555.IL1B","Chr09G000275.C1QA", "Chr15G000627.FCER1G",
             "Chr13G000168.CD74","Chr01G001593.PLD4",
             "Chr14G000706.TREM1", "Chr03G000779.CYBB",
             "Chr02G000866.MKI67","Chr11G000482.TOP2A",
             'Chr11G001548.HBA')
			 
FeaturePlot(scRNA_merge2, features =markers,cols = c('lightgray','orange','red'),ncol = 7,label = F)
DotPlot(scRNA_merge2,features =markers,cols = c('#FCFCFC','#ff3333'),group.by = 'cell_name')+scale_size(range = c(0,5))+theme(axis.text.x = element_text(angle = 45,hjust = 1))

#Cell type specific biological processes
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
top50.cell.list <- split(as.matrix(top50.cell)[,7], top50.cell[,6])
top50.cell.go <- data.frame()
for (i in names(top50.cell.list)){
  go_temp <- enricher(top50.cell.list[[i]],TERM2GENE = BP,pAdjustMethod = "BH",pvalueCutoff=1,qvalueCutoff = 1)
  go_temp1 <- as.data.frame(go_temp)
  go_temp1$cell_name <- rep(i,nrow(go_temp1))
  top50.cell.go <- rbind(top50.cell.go,go_temp1)
}

top50.cell.go1 <- top50.cell.go %>% group_by(cell_name) %>% top_n(n=-100,wt=pvalue)
write.table(top50.cell.go1,'target.item.txt',sep = '\t',quote = F,row.names = F)
target.go <- read_excel('target.item.select.xlsx')
target.go$Description <- factor(target.go$Description,levels = rev(unique(target.go$Description)))
target.go$cell_name <- factor(target.go$cell_name,levels = level_cell)

ggplot(target.go,aes(cell_name,Description))+geom_point(aes(size=Count,color=pvalue))+theme_bw()+
  scale_colour_gradient(low="#CD3700",high="#58A4C3")+theme(axis.title= element_text(size=7, color="black"),
  axis.text.x=element_text(size=7,color = "black",angle = 90,hjust = 1,vjust = 0.4), 
  axis.text.y=element_text(size=7,color = "black"))+scale_size(range=c(1, 3))

ggbarplot(target.go,x = 'Description',y = 'Count',fill = 'cell_name',palette = my36colors)+
  coord_flip()+theme_test()+theme(axis.title= element_text(size=7, color="black"),
        axis.text.x=element_text(size=7,color = "black"), 
        axis.text.y=element_text(size=7,color = "black"))
