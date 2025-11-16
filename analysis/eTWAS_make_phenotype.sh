



dir=("BNaive" "BMem" "CD4Naive" "CD4EM" "CD4Treg" "CD8Naive" "CD8GZMH" "CD8GZMK" "MAIT" "γδT" "NKDim" "NKBright" "cM" "ncM")

for cell in "${dir[@]}"; do
    mkdir -p "$cell"
done

cat > step1.R

args <- commandArgs(trailingOnly = TRUE)
cell <- args[1]
library(Seurat)
library(dplyr)
library(Matrix)
library(qs)
library(data.table)
library(tibble)
library(tidyr)

print(paste("Processing cell type:", cell))
cell_data_path <- paste0(cell, "/Step4.", cell, ".data.qs")
sce1 <- qread(cell_data_path)
cellnumber <- 10
donor_id_counts <- table(sce1$POOL_sampleID)
names_above_10 <- names(donor_id_counts[donor_id_counts >= cellnumber])
sce2 <- subset(sce1, cells = Cells(sce1)[sce1$POOL_sampleID %in% names_above_10])
sce2_counts_df <- as.data.frame(as.matrix(GetAssayData(sce2, slot = "counts")))
colnames(sce2_counts_df) <- sapply(strsplit(colnames(sce2_counts_df), "-"), `[`, 1)
output_count_path <- paste0("./", cell, "/", cell, "_reads_count.txt")
sce2_counts_df$rownames <- rownames(sce2_counts_df)
sce2_counts_df <- sce2_counts_df[, c("rownames", setdiff(names(sce2_counts_df), "rownames"))]
fwrite(sce2_counts_df, output_count_path, sep = "\t", quote = FALSE, col.names = TRUE)
sce2_meta <- sce2@meta.data %>%
  rownames_to_column(var = "original_rowname") %>%
  separate(original_rowname, into = c("pool", "sample", "Barcode"), sep = "@", remove = FALSE) %>%
  mutate(Barcode = sapply(strsplit(Barcode, "-"), `[`, 1)) %>%
  select(-original_rowname)
output_meta_path <- paste0("./", cell, "/", cell, "_meta_data.txt")
write.table(sce2_meta[c("pool", "sample", "Barcode", "age_onek1k", "sex_onek1k")],
            output_meta_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

print(paste("Finished processing cell type:", cell))

Rscript step1.R "BMem"
Rscript step1.R "BNaive"
Rscript step1.R "CD4EM"
Rscript step1.R "CD4Naive"
Rscript step1.R "CD4Treg"
Rscript step1.R "CD8GZMH"
Rscript step1.R "CD8GZMK"
Rscript step1.R "CD8Naive"
Rscript step1.R "MAIT"
Rscript step1.R "NKBright"
Rscript step1.R "NKDim"
Rscript step1.R "cM"
Rscript step1.R "ncM"
Rscript step1.R "γδT"



cat > step2.py
import pandas as pd
import numpy as np
import os
import sys
import qtl.norm

cell_type = sys.argv[1]

work_dir = f"{cell_type}"
os.chdir(work_dir)

meta_data = pd.read_csv(f"{cell_type}_meta_data.txt", sep="\t")
meta_data["id"] = meta_data["Barcode"]
meta_data["sample"] = meta_data["sample"]
meta_data["cell_id"] = meta_data["sample"].astype(str) + "@" + cell_type
meta_data.index = meta_data["pool"] + "@" + meta_data["sample"] + "@" + meta_data["Barcode"]

feather_file = f"{cell_type}_reads_count.feather"
txt_file = f"{cell_type}_reads_count.txt"

if os.path.exists(feather_file):
    print(f"Loading from feather file: {feather_file}")
    count_data = pd.read_feather(feather_file).set_index("index")
else:
    print(f"Reading TXT file: {txt_file}")
    count_data = pd.read_csv(
        txt_file,
        sep="\t",
        index_col=0,
        #dtype='int32',        
        low_memory=False,
        engine="c",            
        memory_map=True,      
        na_filter=False      
    )
    print("Saving as feather for future fast access...")
    count_data.reset_index().to_feather(feather_file)

count_data.index.name = None
count_data = count_data.T
count_data = count_data.loc[meta_data.index, :]
count_data["cell_id"] = list(meta_data["cell_id"])
print("Shape before aggregation:", count_data.shape)

cell_data = count_data.groupby(by="cell_id").sum().T
index = (cell_data == 0).sum(axis=1) < len(cell_data.columns) * 0.6
cell_data = cell_data.loc[index, :]
print("Shape after filtering:", cell_data.shape)

python step2.py "BMem"
python step2.py "BNaive"
python step2.py "CD4EM"
python step2.py "CD4Naive"
python step2.py "CD4Treg"
python step2.py "CD8GZMH"
python step2.py "CD8GZMK"
python step2.py "CD8Naive"
python step2.py "MAIT"
python step2.py "NKBright"
python step2.py "NKDim"
python step2.py "cM"
python step2.py "ncM"
python step2.py "γδT"



cat > step3.R

library(getopt)
library(data.table)
library(stringr)
library(dplyr)
options(future.globals.maxSize = 1000 * 1024^3)
options(stringsAsFactors = FALSE)

spec <- matrix(c(
  'pheno', 'p', 1, "character", "Path to phenotype file (e.g., int_norm.txt)",
  'ref', 'r', 1, "character", "Path to reference gene information file (e.g., gene_info.txt)",
  'cell', 'c', 1, "character", "Cell type (e.g., BMem)"
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

if (is.null(opt$pheno) | is.null(opt$ref) | is.null(opt$cell)) {
  stop("Usage: Rscript step3.R --pheno int_norm.txt --ref gene_info.txt --cell cell_type")
}

cell <- opt$cell

pheno_data <- fread(opt$pheno, header = TRUE)
ref_data <- fread(opt$ref, header = TRUE)

pheno_data <- rename(pheno_data, "geneName" = "V1")
pheno_data <- data.frame(pheno_data, check.names = FALSE)

colnames(ref_data) <- c("chr", "start", "end", "strand", "gene_ESMID", "geneName", "geneType")
ref_data <- ref_data[ref_data$chr %in% paste0("chr", 1:22),]

pheno_data <- data.frame(merge(ref_data, pheno_data), check.names = FALSE)
pheno_data <- pheno_data[, c("chr", "start", "end", "gene_ESMID", colnames(pheno_data)[8:ncol(pheno_data)])]
pheno_data$chr <- str_remove(pheno_data$chr, "chr")
colnames(pheno_data) <- str_split_fixed(colnames(pheno_data), "@", 2)[,1]
colnames(pheno_data) <- c("#chr", "start", "end", "gene_ESMID", colnames(pheno_data)[5:ncol(pheno_data)])

fwrite(pheno_data, paste0("./", cell, "/", cell, "_expr_bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
fwrite(pheno_data[, c("#chr", "start", "end", "gene_ESMID")], paste0("./", cell, "/", cell, "_expr_meta"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


dir=("BNaive" "BMem" "CD4Naive" "CD4EM" "CD4Treg" "CD8Naive" "CD8GZMH" "CD8GZMK" "MAIT" "γδT" "NKDim" "NKBright" "cM" "ncM")

for cell in "${dir[@]}"; do
    echo "Processing $cell"
    Rscript step3.R --pheno "./${cell}/${cell}.txt" --ref "tt1.txt" --cell ${cell}
    bedtools sort -header -i ./${cell}/${cell}_expr_bed
done





