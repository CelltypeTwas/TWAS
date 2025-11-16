

library(stringr)
library(dplyr)
library(data.table)
library(ggplot2)
library(locuscomparer)
library(cowplot)


args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Rscript batch_coloc_analysis.R GENE INTRON TISSUE\n",
       "Rscript batch_coloc_analysis.R AHI1 chr6:135492291:135495313 NKDim")
}

GENE <- args[1]
INTRON <- args[2]
TISSUE <- args[3]

message(GENE, " | ", INTRON, " | ", TISSUE)


output_dir <- paste0("/fig/",
                     GENE, "_", INTRON, "_", TISSUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)


chr_num <- gsub("chr", "", str_extract(INTRON, "chr[0-9]+"))

eqtl_file <- paste0("/",
                    GENE, "/eqtl_", GENE, "_", TISSUE, ".txt")
sqtl_file <- paste0("/",
                    GENE, "/sqtl_", GENE, "_", INTRON, "_", TISSUE, ".txt")
rs_file <- paste0("/merged_chr",
                  chr_num, ".txt")
vcf_file <- paste0("/filter.chr",
                   chr_num, ".dose.vcf.gz")


eqtl <- read.table(eqtl_file, header = FALSE, sep = "", stringsAsFactors = FALSE)
colnames(eqtl) <- c("gene_esmb", "variant_id", "beta_eqtl", "se_eqtl", "pval_eqtl")

sqtl <- read.table(sqtl_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

rs_map <- read.table(rs_file, header = FALSE, sep = "", stringsAsFactors = FALSE)
colnames(rs_map) <- c("pos", "rsid")

merged <- merge(eqtl, sqtl, by = "variant_id")

merged_new <- merged %>%
  select(variant_id, pval_eqtl, pval_nominal) %>%
  mutate(
    `-log10(pval_eqtl)` = -log10(pval_eqtl),
    `-log10(pval_nominal)` = -log10(pval_nominal),
    variant_id2 = sub(":[ACGT]+:[ACGT]+$", "", variant_id)
  )

merged_new$pos <- sub("^chr", "", merged_new$variant_id2)
merged_final2 <- merge(merged_new, rs_map, by = "pos")

merged_final <- merged_final2[, c("variant_id", "rsid", "pval_eqtl", "pval_nominal",
                                   "-log10(pval_eqtl)", "-log10(pval_nominal)", "variant_id2")]

merged_clean <- merged_final %>%
  mutate(
    `Variation ID` = variant_id,
    dbSNP = rsid,
    pos = as.integer(str_extract(variant_id, "(?<=:)[0-9]+")),
    epval = pval_eqtl,
    spval = pval_nominal,
    logpe = `-log10(pval_eqtl)`,
    logps = `-log10(pval_nominal)`
  ) %>%
  select(`Variation ID`, dbSNP, pos, epval, spval, logpe, logps) %>%
  as.data.table()


outfile <- paste0("merged_new_with_rs_chr", chr_num, ".txt")
write.table(merged_clean, outfile, sep = "\t", quote = FALSE, row.names = FALSE)

setDT(merged_clean)
target_id <- merged_clean[which.min(epval * spval), dbSNP]

variant_for_target <- merged_clean$`Variation ID`[merged_clean$dbSNP == target_id]
write.table(variant_for_target, "variant_for_target.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

variant_list <- merged_clean$`Variation ID`
write.table(variant_list, "variant_id.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

system(paste("plink --vcf", vcf_file,
             "--ld-snp-list variant_for_target.txt",
             "--extract variant_id.txt",
             "--r2 --ld-window 1000000 --ld-window-kb 1000 --ld-window-r2 0",
             "--out LD_result > /dev/null 2>&1"))
system("awk '{print $3, $6, $7}' LD_result.ld > LD_result_1")

LD_result_1 <- read.table("LD_result_1", header = TRUE, sep = "", 
                          stringsAsFactors = FALSE)

LD_result_1 <- LD_result_1 %>%
  mutate(
    SNP_A_fmt = str_replace(SNP_A, "^chr", ""),
    SNP_A_fmt = str_replace(SNP_A_fmt, ":[ACGT]+:[ACGT]+$", ""),
    SNP_B_fmt = str_replace(SNP_B, "^chr", ""),
    SNP_B_fmt = str_replace(SNP_B_fmt, ":[ACGT]+:[ACGT]+$", "")
  ) %>%
  left_join(rs_map, by = c("SNP_A_fmt" = "pos")) %>%
  rename(rsid_A = rsid) %>%
  left_join(rs_map, by = c("SNP_B_fmt" = "pos")) %>%
  rename(rsid_B = rsid) %>%
  mutate(
    SNP_A_final = ifelse(!is.na(rsid_A), rsid_A, SNP_A),
    SNP_B_final = ifelse(!is.na(rsid_B), rsid_B, SNP_B)
  ) %>%
  select(SNP_A_final, SNP_B_final, R2)

colnames(LD_result_1) <- c("SNP_A", "SNP_B", "R2")

add_label <- function(merged, snp) {
  merged$label <- ifelse(merged$rsid %in% snp, merged$rsid, '')
  return(merged)
}

merged <- merged_clean
colnames(merged)[which(colnames(merged) == "dbSNP")] <- "rsid"
merged <- merged %>% select(rsid, pos, epval, spval, logpe, logps)
colnames(merged) <- c("rsid", "pos", "pval1", "pval2", "logp1", "logp2")
merged$chr <- chr_num
merged <- merged[!duplicated(merged[, c("pos")]), ]

ld <- as.data.table(LD_result_1)
setnames(ld, c("SNP_A", "SNP_B", "R2"))
ld <- ld[, .(SNP_A = as.character(SNP_A),
             SNP_B = as.character(SNP_B),
             R2 = as.numeric(R2))]

snp <- merged[which.min(merged$pval1 * merged$pval2), 'rsid']
snp <- get_lead_snp(merged, snp)

color <- assign_color(merged$rsid, snp, ld)
shape <- ifelse(merged$rsid == snp, 23, 21)
names(shape) <- merged$rsid
size <- ifelse(merged$rsid == snp, 4, 3)
names(size) <- merged$rsid
merged <- add_label(merged, snp)

p1 <- make_scatterplot(merged, paste(GENE, 'cis-eQTL'), paste(INTRON, 'trans-sQTL'), 
                       color, shape, size, legend = TRUE, 
                       legend_position = c('bottomright', 'topright', 'topleft'))
ggsave("single.pdf", p1, device = "pdf", width = 4, height = 4)

p2 <- make_combined_plot(merged, 'cis-eQTL', 'trans-sQTL', ld, 
                         chr = as.numeric(chr_num), snp = snp, 
                         combine = FALSE, legend = TRUE)

ggsave("locuscompare.pdf", plot = p2$locuscompare, device = "pdf", width = 6, height = 6)
ggsave("locuszoom1.pdf", plot = p2$locuszoom1, device = "pdf", width = 6, height = 6)
ggsave("locuszoom2.pdf", plot = p2$locuszoom2, device = "pdf", width = 6, height = 6)
