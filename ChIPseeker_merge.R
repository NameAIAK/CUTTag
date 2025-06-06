library(ChIPseeker)
# BiocManager::install('ChIPseeker')
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db)
# BiocManager::install("EnsDb.Hsapiens.v75")
library(EnsDb.Hsapiens.v75)
# install.packages('ggimage')
library(ggplot2)
library(dplyr)
edb <- EnsDb.Hsapiens.v75
seqlevelsStyle(edb) <- "UCSC"

# files <- getSampleFiles()
# peak <- readPeakFile(files[[1]])
# peak
library(dplyr)
file1 <- './N87_peak.anno.txt'
file2 <- './AGS_peak.anno.txt'
# 读取两个TXT文件（假设文件没有表头）
df1 <- read.csv(file1,sep = '\t')
df1 <- df1[,c(2,3,4,16,18,19)]
df2 <- read.csv(file2,sep = '\t')
df2 <- df2[,c(2,3,4,16,18,19)]
colnames(df2)
# 根据 chr, pos, pos2 合并数据（保留三列完全匹配的行）
merged_data <- merge(df1, df2, by = c("Gene.Name","Gene.Description","Gene.Type"))
merged_data <- merged_data[,c(4,5,6,1,2,3)]
colnames(merged_data) <- colnames(df1)
data <- distinct(merged_data)
# 查看合并后的前几行
head(data)
write.table(data, "merged_result.txt", row.names = FALSE, quote = FALSE,sep = '\t')
# install.packages ('writexl')
library (writexl)
write_xlsx(data,'merged_result.xlsx ')

data3 <- data[,1:3]
write.table(data3, "mergepeak.txt", row.names = FALSE, col.names = FALSE,quote = FALSE,sep = '\t')


file <- './mergepeak.txt'
peak_file <- readPeakFile(file)
peak_file

peakAnno <- annotatePeak(file, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
# peakAnno.edb <- annotatePeak(file, tssRegion=c(-3000, 3000),
                             # TxDb=edb, annoDb="org.Hs.eg.db")

pdf(file="merge_peakAnnoPie.pdf",width=5.5,height=3)
plotAnnoPie(peakAnno)
dev.off()
pdf(file="merge_peakAnnoUpsetplot.pdf",width=8,height=6)
upsetplot(peakAnno, vennpie=TRUE)
dev.off()

# file <- './N87.peak_summits.bed'
# 
# peakAnno <- annotatePeak(file, tssRegion=c(-3000, 3000),
#                          TxDb=txdb, annoDb="org.Hs.eg.db")
# # peakAnno.edb <- annotatePeak(file, tssRegion=c(-3000, 3000),
# # TxDb=edb, annoDb="org.Hs.eg.db")
# 
# plotAnnoPie(peakAnno)
# upsetplot(peakAnno, vennpie=TRUE)
# plotAnnoBar(peakAnno)

####enrichment####
# BiocManager::install("ReactomePA")
library(ReactomePA)


peakAnno_df <- as.data.frame(peakAnno)
table(peakAnno_df$annotation)
peak_Promoter <- peakAnno_df[grepl("Promoter", peakAnno_df$annotation), ]
table(peak_Promoter$annotation)
# pathway <- enrichPathway(peak_Promoter$geneId)
# pdf(file="N87_peakAnnodotplot.pdf",width=4,height=4)
# dotplot(pathway)
# dev.off()

ego.bp = enrichGO(gene=peak_Promoter$geneId, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE)                    
# write.csv(data.frame(ego.bp@result),file = 'TM_up_GO.csv',row.names = F)
#气泡图#
ego.bp_up <- data.frame(ego.bp@result)
ego.bp_up <- ego.bp_up[ego.bp_up$pvalue<0.05,]
# ego.bp_up <- ego.bp_up[c(1:10),]
# # ego.bp_up$GeneRatio_numeric <- sapply(ego.bp_up$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_up[order(ego.bp_up$GeneRatio_numeric),]$Description
ego.bp_up$Description<- factor(ego.bp_up$Description,levels =ego.bp_up[order(ego.bp_up$Count),]$Description)
write.csv(file="merge_Promter_go_p0.05.csv",data.frame(ego.bp_up),row.names=F)

pdf(file="merge_Promter_go_bar.pdf",width = 6,height = 6)
ggplot(ego.bp_up[1:10,],
       aes(y=Description,x=Count))+
  geom_point(aes(size=Count,color=p.adjust))+
  # facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Count",y="",title="Promter GO Enrichment")+
  theme_bw()
dev.off()


####MC-kegg_up####
kk <- enrichKEGG(gene= peak_Promoter$geneId, organism= 'hsa',pvalueCutoff = 0.05)
kk_up <- data.frame(kk@result)
kk_up <- kk_up[kk_up$pvalue<0.05,]
kk_up <- kk_up %>%
  dplyr::filter(!grepl("disease", Description, ignore.case = TRUE))
# kk_up <- ego.bp_up[c(1:10),]
# # ego.bp_up$GeneRatio_numeric <- sapply(ego.bp_up$GeneRatio, convert_fraction_to_numeric)  
# Kego.bp_up[order(ego.bp_up$GeneRatio_numeric),]$Description

kk_up$Description<- factor(kk_up$Description,levels =kk_up[order(kk_up$Count),]$Description)
write.csv(file="merge_Promter_kegg_p0.05.csv",data.frame(kk_up),row.names=F)

pdf(file="merge_Promter_kegg_bar.pdf",width = 6,height = 6)

ggplot(data = kk_up[1:10,],
       aes(x = Count, y = Description, fill = pvalue)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7) +
  labs(x = "GeneRatio",
       y = "Description",
       title = "Promter_kegg") +
  geom_text(aes(x = 0.002, #用数值向量控制文本标签起始位置
                label = Description),
            hjust = 0)+ #hjust = 0,左对齐
  theme_classic() + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        axis.text.y = element_blank())

dev.off()












