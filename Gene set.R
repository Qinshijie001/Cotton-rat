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

scRNA_merge2 <- readRDS('scRNA_merge2.RDS')
my_comparisons <- list(c('Control','2dpi'),c('Control','5dpi'),c('Control','8dpi'))
my_comparisons1 <- list(c('Control','2dpi'),c('Control','5dpi'),c('Control','8dpi'),
                        c('2dpi','5dpi'),c('2dpi','8dpi'),c('5dpi','8dpi'))
#
target_signature <- "response to virus"
GES_signature <- BP$GENE[which(BP$TERM==target_signature)]
scRNA_merge2 <- AddModuleScore(scRNA_merge2,features = list(unique(GES_signature)),name = target_signature)
a <- scRNA_merge2@meta.data
min1=min(a[,ncol(a)])
max1=max(a[,ncol(a)])
ggboxplot(a, x = "group_state", y = colnames(a)[ncol(a)],fill = "group_state")+
  stat_compare_means(comparisons = my_comparisons1)+guides(fill=F)+scale_fill_manual(values = group_color)
ggsave(paste0(colnames(a)[ncol(a)],'.pdf'),width = 3,height = 4)

ggplot(a,aes(x=cell_name,y=get(colnames(a)[ncol(a)]),fill=group_state))+ylim(min1,max1)+
  labs(y=target_signature)+geom_boxplot(outlier.color = 'gray',outlier.size = 0.2)+ scale_fill_manual(values = group_color)+
  stat_compare_means(aes(group=group_state), method = 'kruskal.test',label = "p.signif",label.y =max1)+
  theme_test()+theme(axis.text.x = element_text(angle = 45,hjust=1,color = 'black'),axis.text.y = element_text(color = 'black'))
ggsave(paste0(colnames(a)[ncol(a)],'group.pdf'),width = 10,height = 2.5)

#
target_signature <- "response to type II interferon"
GES_signature <- BP$GENE[which(BP$TERM==target_signature)]
scRNA_merge2 <- AddModuleScore(scRNA_merge2,features = list(unique(GES_signature)),name = target_signature)
a <- scRNA_merge2@meta.data
min1=min(a[,ncol(a)])
max1=max(a[,ncol(a)])
ggboxplot(a, x = "group_state", y = colnames(a)[ncol(a)],fill = "group_state")+
  stat_compare_means(comparisons = my_comparisons1)+guides(fill=F)+scale_fill_manual(values = group_color)
ggsave(paste0(colnames(a)[ncol(a)],'.pdf'),width = 3,height = 4)

ggplot(a,aes(x=cell_name,y=get(colnames(a)[ncol(a)]),fill=group_state))+ylim(min1,max1)+
  labs(y=target_signature)+geom_boxplot(outlier.color = 'gray',outlier.size = 0.2)+ scale_fill_manual(values = group_color)+
  stat_compare_means(aes(group=group_state), method = 'kruskal.test',label = "p.signif",label.y =max1)+
  theme_test()+theme(axis.text.x = element_text(angle = 45,hjust=1,color = 'black'),axis.text.y = element_text(color = 'black'))
ggsave(paste0(colnames(a)[ncol(a)],'group.pdf'),width = 10,height = 2.5)

#
target_signature <- "response to type I interferon"
GES_signature <- BP$GENE[which(BP$TERM==target_signature)]
scRNA_merge2 <- AddModuleScore(scRNA_merge2,features = list(unique(GES_signature)),name = target_signature)
a <- scRNA_merge2@meta.data
min1=min(a[,ncol(a)])
max1=max(a[,ncol(a)])
ggboxplot(a, x = "group_state", y = colnames(a)[ncol(a)],fill = "group_state")+
  stat_compare_means(comparisons = my_comparisons1)+guides(fill=F)+scale_fill_manual(values = group_color)
ggsave(paste0(colnames(a)[ncol(a)],'.pdf'),width = 3,height = 4)

ggplot(a,aes(x=cell_name,y=get(colnames(a)[ncol(a)]),fill=group_state))+ylim(min1,max1)+
  labs(y=target_signature)+geom_boxplot(outlier.color = 'gray',outlier.size = 0.2)+ scale_fill_manual(values = group_color)+
  stat_compare_means(aes(group=group_state), method = 'kruskal.test',label = "p.signif",label.y =max1)+
  theme_test()+theme(axis.text.x = element_text(angle = 45,hjust=1,color = 'black'),axis.text.y = element_text(color = 'black'))
ggsave(paste0(colnames(a)[ncol(a)],'group.pdf'),width = 10,height = 2.5)

#
target_signature <- "antigen processing and presentation"
GES_signature <- BP$GENE[which(BP$TERM==target_signature)]
scRNA_merge2 <- AddModuleScore(scRNA_merge2,features = list(unique(GES_signature)),name = target_signature)
a <- scRNA_merge2@meta.data
min1=min(a[,ncol(a)])
max1=max(a[,ncol(a)])
ggboxplot(a, x = "group_state", y = colnames(a)[ncol(a)],fill = "group_state")+
  stat_compare_means(comparisons = my_comparisons1)+guides(fill=F)+scale_fill_manual(values = group_color)
ggsave(paste0(colnames(a)[ncol(a)],'.pdf'),width = 3,height = 4)

ggplot(a,aes(x=cell_name,y=get(colnames(a)[ncol(a)]),fill=group_state))+ylim(min1,max1)+
  labs(y=target_signature)+geom_boxplot(outlier.color = 'gray',outlier.size = 0.2)+ scale_fill_manual(values = group_color)+
  stat_compare_means(aes(group=group_state), method = 'kruskal.test',label = "p.signif",label.y =max1)+
  theme_test()+theme(axis.text.x = element_text(angle = 45,hjust=1,color = 'black'),axis.text.y = element_text(color = 'black'))
ggsave(paste0(colnames(a)[ncol(a)],'group.pdf'),width = 10,height = 2.5)

#
scRNA_merge2@meta.data$rsv_group <- 'RSV-'
scRNA_merge2@meta.data$rsv_group[scRNA_merge2@meta.data$RSV_load>0 & scRNA_merge2@meta.data$RSV_load<10 & scRNA_merge2@meta.data$rsv_load_log_soupx>0]='low RSV+'
scRNA_merge2@meta.data$rsv_group[scRNA_merge2@meta.data$RSV_load>=10 & scRNA_merge2@meta.data$rsv_load_log_soupx>0]='high RSV+'

target_signature <- "transcription by RNA polymerase II"
GES_signature <- BP$GENE[which(BP$TERM==target_signature)]
scRNA_merge2 <- AddModuleScore(scRNA_merge2,features = list(unique(GES_signature)),name = target_signature)
a <- scRNA_merge2@meta.data
min1=min(a[,ncol(a)])
max1=max(a[,ncol(a)])
ggboxplot(a, x = "rsv_group", y = colnames(a)[ncol(a)],fill = "rsv_group")+
  stat_compare_means(comparisons = my_comparisons1)+guides(fill=F)+scale_fill_manual(values = group_color)
ggsave(paste0(colnames(a)[ncol(a)],'.pdf'),width = 3,height = 4)

ggplot(a,aes(x=cell_name,y=get(colnames(a)[ncol(a)]),fill=rsv_group))+ylim(min1,max1)+
  labs(y=target_signature)+geom_boxplot(outlier.color = 'gray',outlier.size = 0.2)+ scale_fill_manual(values = group_color)+
  stat_compare_means(aes(group=rsv_group), method = 'kruskal.test',label = "p.signif",label.y =max1)+
  theme_test()+theme(axis.text.x = element_text(angle = 45,hjust=1,color = 'black'),axis.text.y = element_text(color = 'black'))
ggsave(paste0(colnames(a)[ncol(a)],'group.pdf'),width = 10,height = 2.5)

target_signature <- "translational initiation"
GES_signature <- BP$GENE[which(BP$TERM==target_signature)]
scRNA_merge2 <- AddModuleScore(scRNA_merge2,features = list(unique(GES_signature)),name = target_signature)
a <- scRNA_merge2@meta.data
min1=min(a[,ncol(a)])
max1=max(a[,ncol(a)])
ggboxplot(a, x = "rsv_group", y = colnames(a)[ncol(a)],fill = "rsv_group")+
  stat_compare_means(comparisons = my_comparisons1)+guides(fill=F)+scale_fill_manual(values = group_color)
ggsave(paste0(colnames(a)[ncol(a)],'.pdf'),width = 3,height = 4)

ggplot(a,aes(x=cell_name,y=get(colnames(a)[ncol(a)]),fill=rsv_group))+ylim(min1,max1)+
  labs(y=target_signature)+geom_boxplot(outlier.color = 'gray',outlier.size = 0.2)+ scale_fill_manual(values = group_color)+
  stat_compare_means(aes(group=rsv_group), method = 'kruskal.test',label = "p.signif",label.y =max1)+
  theme_test()+theme(axis.text.x = element_text(angle = 45,hjust=1,color = 'black'),axis.text.y = element_text(color = 'black'))
ggsave(paste0(colnames(a)[ncol(a)],'group.pdf'),width = 10,height = 2.5)

target_signature <- "translational termination"
GES_signature <- BP$GENE[which(BP$TERM==target_signature)]
scRNA_merge2 <- AddModuleScore(scRNA_merge2,features = list(unique(GES_signature)),name = target_signature)
a <- scRNA_merge2@meta.data
min1=min(a[,ncol(a)])
max1=max(a[,ncol(a)])
ggboxplot(a, x = "rsv_group", y = colnames(a)[ncol(a)],fill = "rsv_group")+
  stat_compare_means(comparisons = my_comparisons1)+guides(fill=F)+scale_fill_manual(values = group_color)
ggsave(paste0(colnames(a)[ncol(a)],'.pdf'),width = 3,height = 4)

ggplot(a,aes(x=cell_name,y=get(colnames(a)[ncol(a)]),fill=rsv_group))+ylim(min1,max1)+
  labs(y=target_signature)+geom_boxplot(outlier.color = 'gray',outlier.size = 0.2)+ scale_fill_manual(values = group_color)+
  stat_compare_means(aes(group=rsv_group), method = 'kruskal.test',label = "p.signif",label.y =max1)+
  theme_test()+theme(axis.text.x = element_text(angle = 45,hjust=1,color = 'black'),axis.text.y = element_text(color = 'black'))
ggsave(paste0(colnames(a)[ncol(a)],'group.pdf'),width = 10,height = 2.5)
