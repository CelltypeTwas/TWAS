

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
sce = readRDS("sce_filtered.rds")
seurat.qc <- subset(sce, subset = nFeature_RNA >= 500 & 
                      nCount_RNA>= 1000 & percent.mt <= 20)
table(is.na(seurat.qc@meta.data$POOL_sampleID))
samples_to_remove <- c("B24@319_320", "B24@325_326", "B77@193_194")
cells_to_keep <- rownames(seurat.qc@meta.data)[
  !(seurat.qc@meta.data$POOL_sampleID %in% samples_to_remove)
]
seurat.qc.filtered <- subset(seurat.qc, cells = cells_to_keep)
qsave(seurat.qc.filtered, "seurat.qc.filtered.qs")
