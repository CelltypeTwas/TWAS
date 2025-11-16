
VERSION <- "1.0.0"
USAGE <- "
COLOC Analysis Pipeline (v%s)

Usage:
  Rscript coloc_analysis.R <file_name> <gene_name> <celltype> <chromosome>

Arguments:
  <file_name>    Trait name (e.g., BMI, MS)
  <gene_name>    Target gene (e.g., PRR14)
  <celltype>     Cell type (e.g., NKBright)
  <chromosome>   Chromosome (e.g., chr16)

Output:
  Results saved in 'result/' directory
"

cc_trait <- c("MS","SLE","T1D","ovariancancer","RA","CD","UC",
              "endometrialcarcinoma","Lungadenocarcinoma","lungcarcinoma",
              "smallcelllungcarcinoma","Squamouscelllungcarcinoma",
              "Breastcancer","Kidneycancer","cervicalcancer","GD","SS",
              "Type2diabetes","Basalcellcarcinoma","squamouscellcarcinoma",
              "PA","hyperlipidemia","asthma","Prostatecancer","gout")

cc_trait_case <- c(4888, 5201, 18942, 15588, 22350, 20873, 23252, 8758, 
                   11273, 29266, 2664, 7426, 76192, 29020, 363, 4487, 1599, 
                   84224, 20791, 7402, 5065, 17485, 56167, 122188, 37105)

cc_trait_control <- c(10395, 9066, 501638, 105724, 74823, 346719, 352256, 
                      46126, 55483, 56450, 21444, 55627, 63082, 835670, 861, 
                      629598, 658316, 583280, 286893, 286892, 21286, 331737, 
                      352255, 604640, 1448128)

quant_trait <- c("Hypertension","BAS","EOS","Hb","Ht","LYM","MCHC","MCH", 
                 "MCV", "MON","NEU", "PLT", "RBC", "WBC", "BMI","height",
                 "Metabolicsyndrome")

quant_trait_num <- c(338391, 474001, 474237, 563946, 562259, 524923, 491553, 
                     486823, 544127, 349856, 519288, 542827, 545203, 562243, 
                     523818, 360388, 1384348)

cell_type <- c("BNaive","BMem", "CD4Naive", "CD4EM", "CD4Treg", "CD8GZMH", 
               "CD8GZMK", "CD8Naive", "MAIT", "γδT", "NKDim", "NKBright", 
               "cM", "ncM")

celltype_sample <- c(873, 866, 972, 977, 748, 921, 934, 755, 296, 255, 974, 
                     289, 569, 446)


main <- function() {

  options(stringsAsFactors = FALSE)
  suppressPackageStartupMessages(library(coloc))
  

  args <- commandArgs(TRUE)
  
  if ("--help" %in% args || "-h" %in% args) {
    cat(sprintf(USAGE, VERSION))
    quit(status = 0)
  }

  if (length(args) != 4) {
    cat("Error: Incorrect number of arguments\n")
    cat(sprintf(USAGE, VERSION))
    quit(status = 1)
  }
  
  file_name <- args[1]
  gene_name <- args[2]
  celltype <- args[3]
  chromosome <- args[4]
  
  num_1 <- which(cc_trait == file_name)
  num_2 <- which(quant_trait == file_name)
  
  if (length(num_1) > 0) {
    trait_type <- "cc"
    trait_sample_size <- cc_trait_case[num_1] + cc_trait_control[num_1]
    trait_s <- cc_trait_case[num_1] / trait_sample_size
  } else if (length(num_2) > 0) {
    trait_type <- "quant"
    trait_sample_size <- quant_trait_num[num_2]
  } else {
    stop("Error: Unrecognized trait name '", file_name, "'")
  }

  num_3 <- which(cell_type == celltype)
  if (length(num_3) == 0) {
    stop("Error: Unrecognized cell type '", celltype, "'")
  }
  
  chromosome_number <- as.numeric(sub("chr", "", chromosome, ignore.case = TRUE))
  if (is.na(chromosome_number)) {
    stop("Error: Invalid chromosome format '", chromosome, "'. Expected format like 'chr16'")
  }


  base_dir <- "/coloc"
  
  file_path_gwas <- file.path(base_dir, celltype, file_name, "coloc", "gwas", 
                             paste0(gene_name, "_gwas1"))
  file_path_eqtl <- file.path(base_dir, celltype, file_name, "coloc", "eqtl", 
                             gene_name)
  file_path_maf <- file.path("/rename_maf", 
                            paste0(chromosome, "_all"))
  file_path_ref <- file.path("/output", 
                           paste0("merged_", chromosome, ".txt"))

  required_files <- c(file_path_gwas, file_path_eqtl, file_path_maf, file_path_ref)
  missing_files <- required_files[!file.exists(required_files)]
  
  if (length(missing_files) > 0) {
    stop("Error: Missing required files:\n  ", paste(missing_files, collapse = "\n  "))
  }
  
  safe_read <- function(path, ...) {
    tryCatch(
      {
        data <- read.table(path, ...)
        cat("   ✔ Loaded:", path, "\n")
        return(data)
      },
      error = function(e) {
        stop("Error loading ", path, ":\n", e$message)
      }
    )
  }

  data_gwas <- safe_read(file_path_gwas, header = FALSE, sep = " ")
  colnames(data_gwas) <- c("snp_gwas", "ref_gwas", "alt_gwas", "beta_gwas", 
                          "se_gwas", "P_gwas", "se2_gwas")
  data_gwas$snp_gwas <- gsub("_", ":", data_gwas$snp_gwas)
  data_gwas$snp_gwas <- sub(":b38", "", data_gwas$snp_gwas)


  data_eqtl <- safe_read(file_path_eqtl, header = FALSE, sep = " ")
  colnames(data_eqtl) <- c("esmble_eqtl", "snp_eqtl", "pval_eqtl", "beta_eqtl", 
                          "se_eqtl")
  data_eqtl$snp_eqtl <- gsub("chr", "", data_eqtl$snp_eqtl)

  data_maf <- safe_read(file_path_maf, header = FALSE, sep = "")
  colnames(data_maf) <- c("snp_maf", "maf")
  
  data_ref <- safe_read(file_path_ref, header = FALSE, sep = " ")
  colnames(data_ref) <- c("snp_pos", "snp_rs")
  

  merged_data <- merge(data_gwas, data_eqtl, 
                       by.x = "snp_gwas", by.y = "snp_eqtl", 
                       all = FALSE)
  cat("   ✔ GWAS + eQTL: ", nrow(merged_data), "SNPs\n")

  merged_maf <- merge(merged_data, data_maf, 
                      by.x = "snp_gwas", by.y = "snp_maf", 
                      all = FALSE)
  merged_maf$snp_ref <- sub(":[^:]+:[^:]+$", "", merged_maf$snp_gwas)
  cat("   ✔ + MAF data: ", nrow(merged_maf), "SNPs\n")

  merged_final <- merge(merged_maf, data_ref, 
                        by.x = "snp_ref", by.y = "snp_pos", 
                        all = FALSE)
  cat("   ✔ + Reference: ", nrow(merged_final), "final SNPs\n")

  if (nrow(merged_final) < 10) {
    stop("Error: Only ", nrow(merged_final), " overlapping SNPs found (minimum 10 required)")
  }


  tryCatch(
    {
      if (trait_type == "cc") {
        result <- coloc.abf(
          dataset1 = list(
            pvalues = merged_final$P_gwas,
            type = "cc",
            beta = merged_final$beta_gwas,
            varbeta = merged_final$se2_gwas,
            s = trait_s,
            N = trait_sample_size,
            snp = merged_final$snp_rs
          ),
          dataset2 = list(
            pvalues = merged_final$pval_eqtl,
            type = "quant",
            N = celltype_sample[num_3],
            snp = merged_final$snp_rs
          ),
          MAF = merged_final$maf
        )
      } else {
        result <- coloc.abf(
          dataset1 = list(
            pvalues = merged_final$P_gwas,
            type = "quant",
            beta = merged_final$beta_gwas,
            varbeta = merged_final$se2_gwas,
            N = trait_sample_size,
            snp = merged_final$snp_rs
          ),
          dataset2 = list(
            pvalues = merged_final$pval_eqtl,
            type = "quant",
            N = celltype_sample[num_3],
            snp = merged_final$snp_rs
          ),
          MAF = merged_final$maf
        )
      }
      
      cat("   ✔ Analysis completed successfully\n")
    },
    error = function(e) {
      stop("COLOC analysis failed:\n", e$message)
    }
  )
  
  if (!dir.exists("result")) {
  dir.create("result")
}

output_dir <- file.path("result", celltype, file_name)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

output_file <- file.path(output_dir, paste0("coloc_", 
                                            trait_type, "@", 
                                            celltype, "@", 
                                            file_name, "@", 
                                            gene_name, "@.txt"))
  
tryCatch(
  {
    write.table(result$summary, file = output_file, 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    cat("   ✔ Results saved to:", output_file, "\t")
  },
  error = function(e) {
    stop("Error saving results:\n", e$message)
  }
)

}

tryCatch(
  {
    main()
  },
  error = function(e) {
    cat(e$message, "\n\n")
    cat("For help, run: Rscript coloc_analysis.R --help\n")
    quit(status = 1)
  }
)