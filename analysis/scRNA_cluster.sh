
cat > step1_harmony_1step.R

library(optparse)
option_list = list(
  make_option("--dimss", action="store", default=0, type='integer', help="dimss for number")
)

opt = parse_args(OptionParser(option_list=option_list))

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
sce <- qread("seurat.qc.filtered.qs")
sce <- NormalizeData(object = sce, normalization.method = "LogNormalize", scale.factor = 1e4)
sce <- FindVariableFeatures(object = sce, selection.method = 'mean.var.plot', mean.cutoff = c(0.1, Inf), dispersion.cutoff = c(0.5, Inf), nfeatures = 2000)
sce <- ScaleData(object = sce, features = rownames(sce))
sce <- RunPCA(object = sce, features = VariableFeatures(object = sce), npcs = 50)
sce <- sce %>% 
  RunHarmony(group.by.vars = "POOL")

qsave(sce, file = "sce_after_harmony.qs")

sce <- sce %>% 
  FindNeighbors(reduction = "harmony", dims = 1:opt$dimss)%>% 
  RunUMAP(reduction = "harmony", dims = 1:opt$dimss, verbose = F)

p1 <- DimPlot(sce, reduction = "umap", group.by = "POOL", raster = F)
p2 <- DimPlot(sce, reduction = "umap",group.by = "l2_onek1k", label = TRUE, raster = F)
p3 <- DimPlot(sce, reduction = "umap",group.by = "cell_type_onek1k", label = TRUE, raster = F)
P.total = wrap_plots(p1,p2,p3)
ggsave(plot = P.total, filename = paste0(opt$dimss,"Umap_plot.png"),width = 48,height = 6)

for (res in c(1.2)) {
  print(res)
  sce <- FindClusters(sce,resolution = res, 
                      algorithm = 1)
}
qsave(sce, file=paste0(opt$dimss,"findcluster.qs"))

cluster_umap <- wrap_plots(
  DimPlot(sce, reduction = "umap", 
          group.by = "RNA_snn_res.1.2", 
          label = TRUE, raster = FALSE) & NoAxes(), 
  ncol = 1
)

ggsave(cluster_umap,filename = paste0(opt$dimss,"_batch_cov","Umap_plot_2.png"),width = 30, height = 10)



nohup Rscript step1_harmony_1step.R --dimss 30 > harmony_log_30.txt 2>&1 &
