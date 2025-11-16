


cat > GWAS_CD4EM_splicing.sh


nohup bash GWAS_CD4EM_splicing.sh > GWAS_CD4EM_splicing.log 2>&1 &


cat > gwas_processing_splicing.sh
TRAIT_NAME=$1
GWAS_FILE=$2
MODEL_NAME=$3
N=$4

GWAS_DIR="GWAS_catalog"
REF_FILE="TWAS_1"
ANNOTATION_FILE="merged_all_chr.txt"
OUTPUT_DIR="${GWAS_DIR}/${TRAIT_NAME}"
MODEL_DB_PATH="${MODEL_NAME}/model/P/X5/data/output/dbs/gEUVADIS_HapMap_alpha0.5_window1e6_filtered.db"
COVARIANCE_PATH="${MODEL_NAME}/model/P/X5/data/output/allCovariances/gEUVADIS_HapMap_alpha0.5_window1e6.txt.gz"
GWAS_FOLDER="${OUTPUT_DIR}/GWAS"

mkdir -p "$OUTPUT_DIR" && cd "$OUTPUT_DIR"

echo "Running SPrediXcan analysis..."
/MetaXcan/MetaXcan-master/software/SPrediXcan.py \
    --model_db_path "$MODEL_DB_PATH" \
    --covariance "$COVARIANCE_PATH" \
    --gwas_folder "$GWAS_FOLDER" \
    --gwas_file_pattern ".*gz" \
    --beta_column BETA \
    --pvalue_column P \
    --snp_column SNP  \
    --non_effect_allele_column A2 \
    --effect_allele_column A1 \
    --keep_non_rsid \
    --output_file "${TRAIT_NAME}_${MODEL_NAME}_C6_splicing.csv"

echo "All processing completed for ${TRAIT_NAME} with model ${MODEL_NAME}."




bash gwas_processing_splicing.sh MS AndlauerTF_SciAdv_2016_MS.txt.gz CD4EM 15283
bash gwas_processing_splicing.sh SLE BenthamJ_NatGenet_2015_SLE.txt.gz CD4EM 14267 
bash gwas_processing_splicing.sh Hypertension CareyCE_NatHumBehav_2024_Hypertension.txt.gz CD4EM 338391
bash gwas_processing_splicing.sh BAS ChenMH_CELL_2020_BAS.txt.gz CD4EM 474001
bash gwas_processing_splicing.sh EOS ChenMH_CELL_2020_EOS.txt.gz CD4EM 474237
bash gwas_processing_splicing.sh Hb ChenMH_CELL_2020_Hb.txt.gz CD4EM 563946
bash gwas_processing_splicing.sh Ht ChenMH_CELL_2020_Ht.txt.gz CD4EM 562259
bash gwas_processing_splicing.sh LYM ChenMH_CELL_2020_LYM.txt.gz CD4EM 524923
bash gwas_processing_splicing.sh MCHC ChenMH_CELL_2020_MCHC.txt.gz CD4EM 491553
bash gwas_processing_splicing.sh MCH ChenMH_CELL_2020_MCH.txt.gz CD4EM 486823
bash gwas_processing_splicing.sh MCV ChenMH_CELL_2020_MCV.txt.gz CD4EM 544127
bash gwas_processing_splicing.sh MON ChenMH_CELL_2020_MON.txt.gz CD4EM 349856
bash gwas_processing_splicing.sh NEU ChenMH_CELL_2020_NEU.txt.gz CD4EM 519288
bash gwas_processing_splicing.sh PLT ChenMH_CELL_2020_PLT.txt.gz CD4EM 542827
bash gwas_processing_splicing.sh RBC ChenMH_CELL_2020_RBC.txt.gz CD4EM 545203
bash gwas_processing_splicing.sh WBC ChenMH_CELL_2020_WBC.txt.gz CD4EM 562243
bash gwas_processing_splicing.sh T1D ChiouJ_Nature_2021_T1D.txt.gz CD4EM 520580
bash gwas_processing_splicing.sh ovariancancer DarengEO_AmJHumGenet_2024_ovariancancer.txt.gz CD4EM 121312
bash gwas_processing_splicing.sh RA Ishigaki_NatGenet_2022_RA.txt.gz CD4EM 97173
bash gwas_processing_splicing.sh CD LiuZ_NatGenet_2023_CD.txt.gz CD4EM 367592
bash gwas_processing_splicing.sh UC LiuZ_NatGenet_2023_UC.txt.gz CD4EM 375508
bash gwas_processing_splicing.sh endometrialcarcinoma MaraTA_NatCommun_2018_endometrialcarcinoma.txt.gz CD4EM 54884
bash gwas_processing_splicing.sh Lungadenocarcinoma McKayJD_NatGenet_2017_Lungadenocarcinoma.txt.gz CD4EM 66756
bash gwas_processing_splicing.sh lungcarcinoma McKayJD_NatGenet_2017_lungcarcinoma.txt.gz CD4EM 85716
bash gwas_processing_splicing.sh smallcelllungcarcinoma McKayJD_NatGenet_2017_smallcelllungcarcinoma.txt.gz CD4EM 24108
bash gwas_processing_splicing.sh Squamouscelllungcarcinoma McKayJD_NatGenet_2017_Squamouscelllungcarcinoma.txt.gz CD4EM 63053
bash gwas_processing_splicing.sh Breastcancer MichailidouK_Nature_2017_Breastcancer.txt.gz CD4EM 139274
bash gwas_processing_splicing.sh Metabolicsyndrome ParkS_NatGenet_2024_Metabolicsyndrome.txt.gz CD4EM 1384348
bash gwas_processing_splicing.sh Kidneycancer PurdueMP_NatGenet_2024_Kidneycancer.txt.gz CD4EM 864690
bash gwas_processing_splicing.sh cervicalcancer RamachandranD_HumMolGenet_2022_cervicalcancer.txt.gz CD4EM 1224
bash gwas_processing_splicing.sh BMI SakaueS_NatGenet_2021_BMI.txt.gz CD4EM 523818
bash gwas_processing_splicing.sh GD SakaueS_NatGenet_2021_GD.txt.gz CD4EM 634085
bash gwas_processing_splicing.sh height SakaueS_NatGenet_2021_height.txt.gz CD4EM 525444
bash gwas_processing_splicing.sh SS SakaueS_NatGenet_2021_SS.txt.gz CD4EM 659915
bash gwas_processing_splicing.sh Type2diabetes SakaueS_NatGenet_2021_Type2diabetes.txt.gz CD4EM 667504
bash gwas_processing_splicing.sh Basalcellcarcinoma SeviiriM_NatCommun_2022_Basalcellcarcinoma.txt.gz CD4EM 307684
bash gwas_processing_splicing.sh squamouscellcarcinoma SeviiriM_NatCommun_2022_squamouscellcarcinoma.txt.gz CD4EM 294294
bash gwas_processing_splicing.sh PA SoomroM_ArthritisRheumatol_2022_PA.txt.gz CD4EM 26351
bash gwas_processing_splicing.sh hyperlipidemia TrinderM_Atherosclerosis_2021_hyperlipidemia.txt.gz CD4EM 349222
bash gwas_processing_splicing.sh asthma ValetteK_CommunBiol_2021_asthma.txt.gz CD4EM 408442
bash gwas_processing_splicing.sh Prostatecancer WangA_NatGenet_2023_Prostatecancer.txt.gz CD4EM 726828
bash gwas_processing_splicing.sh gout ZhouW_CellGenom_2022_gout.txt.gz CD4EM 1485233
