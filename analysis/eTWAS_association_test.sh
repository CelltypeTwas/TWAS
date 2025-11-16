


cat > gwas_processing.sh

TRAIT_NAME=$1
GWAS_FILE=$2
MODEL_NAME=$3
N=$4

GWAS_DIR="/GWAS_catalog"
REF_FILE="/GWAS_catalog/TWAS_1"
ANNOTATION_FILE="merged_all_chr.txt"
OUTPUT_DIR="${GWAS_DIR}/${TRAIT_NAME}"
MODEL_DB_PATH="/${MODEL_NAME}/model/P/X5/data/output/dbs/gEUVADIS_HapMap_alpha0.5_window1e6_filtered.db"
COVARIANCE_PATH="/${MODEL_NAME}/model/P/X5/data/output/allCovariances/gEUVADIS_HapMap_alpha0.5_window1e6.txt.gz"
GWAS_FOLDER="${OUTPUT_DIR}/GWAS"

mkdir -p "$OUTPUT_DIR" && cd "$OUTPUT_DIR"

echo "Preprocessing GWAS File..."
zcat "/GWAS/summary/${GWAS_FILE}" | \
awk -v ref="$REF_FILE" '
BEGIN {
    while ((getline < ref) > 0) {
        twas[$1] = 1
    }
}
NR == 1 { next }
{
    split($2, a, ":")
    if (length(a) != 4) next

    key1 = a[1] ":" a[2] ":" a[3] ":" a[4]
    key2 = a[1] ":" a[2] ":" a[4] ":" a[3]

    if (key1 in twas) {
        print $1, $2, $3, $4, $5
    } else if (key2 in twas) {
        neg_beta = ($3 == "NA") ? "NA" : -1 * $3
        newkey = a[1] ":" a[2] ":" a[4] ":" a[3]
        print $1, newkey, neg_beta, $4, $5
    }
}' | tr ' ' '\t' > "${TRAIT_NAME}_1"

cat "${TRAIT_NAME}_1" | awk '{print $2}' > "${TRAIT_NAME}_GWAS_3"
awk '{ split($2, a, ":"); $1 = a[1]":"a[2]; print }' OFS='\t' "${TRAIT_NAME}_1" > "${TRAIT_NAME}_1_modified.txt"

echo "Annotating GWAS positions..."
awk 'FNR==NR {map[$1]=$2; next} ($1 in map) {print $0, map[$1]}' \
    "$ANNOTATION_FILE" \
    "${TRAIT_NAME}_1_modified.txt" | tr ' ' '\t' > "${TRAIT_NAME}_1_annotated.txt"

echo "Building GWAS output..."
awk '{print $2}' "${TRAIT_NAME}_1_annotated.txt" | awk -F ':' '{print $1,$2,$3,$4}' > "${TRAIT_NAME}_part1"
awk '{print $2}' "${TRAIT_NAME}_1_annotated.txt" | awk -F ':' '{print $1"_"$2"_"$3"_"$4"_b38"}' > "${TRAIT_NAME}_SNP"

(echo -e "CHR\tSNP\tBP\tA1\tA2\tN\tBETA\tSE\tP"; \
paste "${TRAIT_NAME}_part1" "${TRAIT_NAME}_1_annotated.txt" "${TRAIT_NAME}_SNP" | \
awk -v N="$N" '{print $1, $NF, $2, $4, $3, N, $7, $8, $9}') | tr ' ' '\t' > "${TRAIT_NAME}_final_output.txt"

echo "Compressing and moving GWAS results..."
mkdir -p "$GWAS_FOLDER"
gzip -c "${TRAIT_NAME}_final_output.txt" > "${GWAS_FOLDER}/${TRAIT_NAME}_MM.gz"

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
    --output_file "${TRAIT_NAME}_${MODEL_NAME}_C6.csv"

echo "All processing completed for ${TRAIT_NAME} with model ${MODEL_NAME}."


bash gwas_processing.sh MS AndlauerTF_SciAdv_2016_MS.txt.gz CD4EM 15283
bash gwas_processing.sh SLE BenthamJ_NatGenet_2015_SLE.txt.gz CD4EM 14267 
bash gwas_processing.sh Hypertension CareyCE_NatHumBehav_2024_Hypertension.txt.gz CD4EM 338391
bash gwas_processing.sh BAS ChenMH_CELL_2020_BAS.txt.gz CD4EM 474001
bash gwas_processing.sh EOS ChenMH_CELL_2020_EOS.txt.gz CD4EM 474237
bash gwas_processing.sh Hb ChenMH_CELL_2020_Hb.txt.gz CD4EM 563946
bash gwas_processing.sh Ht ChenMH_CELL_2020_Ht.txt.gz CD4EM 562259
bash gwas_processing.sh LYM ChenMH_CELL_2020_LYM.txt.gz CD4EM 524923
bash gwas_processing.sh MCHC ChenMH_CELL_2020_MCHC.txt.gz CD4EM 491553
bash gwas_processing.sh MCH ChenMH_CELL_2020_MCH.txt.gz CD4EM 486823
bash gwas_processing.sh MCV ChenMH_CELL_2020_MCV.txt.gz CD4EM 544127
bash gwas_processing.sh MON ChenMH_CELL_2020_MON.txt.gz CD4EM 349856
bash gwas_processing.sh NEU ChenMH_CELL_2020_NEU.txt.gz CD4EM 519288
bash gwas_processing.sh PLT ChenMH_CELL_2020_PLT.txt.gz CD4EM 542827
bash gwas_processing.sh RBC ChenMH_CELL_2020_RBC.txt.gz CD4EM 545203
bash gwas_processing.sh WBC ChenMH_CELL_2020_WBC.txt.gz CD4EM 562243
bash gwas_processing.sh T1D ChiouJ_Nature_2021_T1D.txt.gz CD4EM 520580
bash gwas_processing.sh ovariancancer DarengEO_AmJHumGenet_2024_ovariancancer.txt.gz CD4EM 121312
bash gwas_processing.sh RA Ishigaki_NatGenet_2022_RA.txt.gz CD4EM 97173
bash gwas_processing.sh CD LiuZ_NatGenet_2023_CD.txt.gz CD4EM 367592
bash gwas_processing.sh UC LiuZ_NatGenet_2023_UC.txt.gz CD4EM 375508
bash gwas_processing.sh endometrialcarcinoma MaraTA_NatCommun_2018_endometrialcarcinoma.txt.gz CD4EM 54884
bash gwas_processing.sh Lungadenocarcinoma McKayJD_NatGenet_2017_Lungadenocarcinoma.txt.gz CD4EM 66756
bash gwas_processing.sh lungcarcinoma McKayJD_NatGenet_2017_lungcarcinoma.txt.gz CD4EM 85716
bash gwas_processing.sh smallcelllungcarcinoma McKayJD_NatGenet_2017_smallcelllungcarcinoma.txt.gz CD4EM 24108
bash gwas_processing.sh Squamouscelllungcarcinoma McKayJD_NatGenet_2017_Squamouscelllungcarcinoma.txt.gz CD4EM 63053
bash gwas_processing.sh Breastcancer MichailidouK_Nature_2017_Breastcancer.txt.gz CD4EM 139274
bash gwas_processing.sh Metabolicsyndrome ParkS_NatGenet_2024_Metabolicsyndrome.txt.gz CD4EM 1384348
bash gwas_processing.sh Kidneycancer PurdueMP_NatGenet_2024_Kidneycancer.txt.gz CD4EM 864690
bash gwas_processing.sh cervicalcancer RamachandranD_HumMolGenet_2022_cervicalcancer.txt.gz CD4EM 1224
bash gwas_processing.sh BMI SakaueS_NatGenet_2021_BMI.txt.gz CD4EM 523818
bash gwas_processing.sh GD SakaueS_NatGenet_2021_GD.txt.gz CD4EM 634085
bash gwas_processing.sh height SakaueS_NatGenet_2021_height.txt.gz CD4EM 525444
bash gwas_processing.sh SS SakaueS_NatGenet_2021_SS.txt.gz CD4EM 659915
bash gwas_processing.sh Type2diabetes SakaueS_NatGenet_2021_Type2diabetes.txt.gz CD4EM 667504
bash gwas_processing.sh Basalcellcarcinoma SeviiriM_NatCommun_2022_Basalcellcarcinoma.txt.gz CD4EM 307684
bash gwas_processing.sh squamouscellcarcinoma SeviiriM_NatCommun_2022_squamouscellcarcinoma.txt.gz CD4EM 294294
bash gwas_processing.sh PA SoomroM_ArthritisRheumatol_2022_PA.txt.gz CD4EM 26351
bash gwas_processing.sh hyperlipidemia TrinderM_Atherosclerosis_2021_hyperlipidemia.txt.gz CD4EM 349222
bash gwas_processing.sh asthma ValetteK_CommunBiol_2021_asthma.txt.gz CD4EM 408442
bash gwas_processing.sh Prostatecancer WangA_NatGenet_2023_Prostatecancer.txt.gz CD4EM 726828
bash gwas_processing.sh gout ZhouW_CellGenom_2022_gout.txt.gz CD4EM 1485233

