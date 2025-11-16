


rm(list=ls())
gc()
library(Seurat)
library(patchwork)
library(harmony)
library(tidyverse)
library(cowplot)
library(qs)
library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)
options(future.globals.maxSize = 1000 * 1024^3)
sce = qread("findcluster.qs")

sce@meta.data$RNA_snn_res.1.2 <- factor(sce@meta.data$RNA_snn_res.1.2, levels = as.character(0:37))
Idents(sce) = "RNA_snn_res.1.2"

genes_to_check = c('CD19', 'CD79A', 'CD79B', 'MS4A1',"CD22","IGHM",'BANK1', #B
                   "MZB1","PRDM1","JCHAIN","SDC1",'IGHG1','IGLL5', #B Plasma cells
                   'PTPRC', ##immune
                   'CD3D', 'CD3E', "CD2",'TRAC','IL32', #T cells
                   'IL7R','CD4','CD40LG',"LEF1","CCR7",'RTKN2', # CD4 T
                   'CD8A','CD8B',"GZMK","GZMA","GZMH", #CD8 T
                   'SELL','IL2RA','IL7R','FOXP3',  #Foxp3
                   'GNLY', 'NKG7',"KLRF1","KLRD1","PRF1","GZMB","IFIT1",'KLRB1','NCR1','FGFBP2','CX3CR1', #NK cells
                   'NCAM1','FCGR3A',   #NK2 cells
                   'CST3', 'LYZ', 'CD68', 'CD163','MMP19',"ITGAX","AIF1","FCER1G", #mono
                   "CD14", "S100A8","S100A9","FCN1", #CD14+ monocyte
                   "FCGR3A","CTSS","CDKN1C", #CD16+ monocyte
                   "APOE","C1QA",'C1QB', #Macrophages
                   "CD1C", 'HLA-DQA1', #cDC
                   "CLEC9A","IRF8","FLT3", #cCD1
                   "CLEC10A","CLEC4C","CD1E","CD1C","FCGR1A", #cCD2
                   "CCR7","LAMP3", "IDO1", 'IDO2', #cDC3
                   "LILRA4","SHD","LRRC26","PACSIN1","IL3RA","SLC32A1","IRF7","IFI44L","PLD4","SERPINF1", #pDCs
                   "PRSS57","SOX4", #progenitor cells (Progen) in PBMC PMID: 35389781
                   "MKI67","TOP2A",'KLRC1','STMN1') #Prolif

genes_to_check= unique(c(genes_to_check))
genes_to_check=str_to_upper(genes_to_check)

genes_to_check
p_mark <- DotPlot(sce, features = unique(genes_to_check), assay='RNA',group.by ="RNA_snn_res.1.2") + coord_flip()
ggsave('plot2.pdf',height = 15,width = 30,limitsize = FALSE)

p1_umap_dim=DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.1.2", label = T, raster=F)
p1_umap_dim|p_mark
ggsave('1.2.pdf',width = 30,height = 15,limitsize = FALSE)

celltype=data.frame(ClusterID=c(0:37),celltype='NA')

celltype[celltype$ClusterID %in% c(1,2,3,4,5,6,8,9,10,11,14,16,17,18,22,25,29,31,35,36,37),2]='T'
celltype[celltype$ClusterID %in% c(0,15,31),2]='NK'
celltype[celltype$ClusterID %in% c(7,12,20,21,26),2]='B'
celltype[celltype$ClusterID %in% c(24),2]='PB'
celltype[celltype$ClusterID %in% c(13,19,30,32,34),2]='Mono'
celltype[celltype$ClusterID %in% c(23,27),2]='DC'
celltype[celltype$ClusterID %in% c(28),2]='HSP'
celltype[celltype$ClusterID %in% c(33),2]='PLT'


colnames(celltype) = c("ClusterID","celltype_main")
table(celltype$celltype_main)
dim(celltype)

sce@meta.data$CT_res1.2 = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@active.ident == celltype$ClusterID[i]),'CT_res1.2'] <- celltype$celltype[i]}
table(sce@meta.data$CT_res1.2)

sce$CT_res1.2 = factor(sce$CT_res1.2, levels = c("B", "PB", 
                                                 "T", "NK",
                                                 "Mono", "DC",
                                                 "HSP", "PLT"))


qsave(sce, file = "all_res_1.2.qs")

table(sce$CT_res1.2)
seurat.Mye = subset(sce, CT_res1.2 %in% c("Mono","DC"))
dir.create("1.Mye.sub/Rawdata/",recursive = T)
qsave(seurat.Mye, file = "./1.Mye.sub/Rawdata/Step4.Mye.data.qs")

seurat.T = subset(sce, CT_res1.2 %in% c("T","NK"))
dir.create("2.T.sub/Rawdata/",recursive = T)
qsave(seurat.T, file = "./2.T.sub/Rawdata/Step4.T.data.qs")

seurat.B = subset(sce, CT_res1.2 %in% c("B","PB"))
dir.create("3.B.sub/Rawdata/",recursive = T)
qsave(seurat.B, file = "./3.B.sub/Rawdata/Step4.B.data.qs")

seurat.HSP = subset(sce, CT_res1.2 %in% c("HSP"))
dir.create("4.HSP.sub/Rawdata/",recursive = T)
qsave(seurat.HSP, file = "./4.HSP.sub/Rawdata/Step4.HSP.data.qs")

seurat.PLT = subset(sce, CT_res1.2 %in% c("PLT"))
dir.create("6.PLT.sub/Rawdata/",recursive = T)
qsave(seurat.PLT, file = "./6.PLT.sub/Rawdata/Step4.PLT.data.qs")

table(sce$CT_res1.2)
seurat.Mye
seurat.T
seurat.B
seurat.HSP
seurat.DL
seurat.PLT

write.csv(sce@meta.data, file = "metadata_annotation.csv")
