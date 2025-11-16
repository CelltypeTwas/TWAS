



trait_file="use.trait.txt"
cells=("BNaive" "BMem" "CD4Naive" "CD4EM" "CD4Treg" "CD8Naive" "CD8GZMH" "CD8GZMK" "MAIT" "γδT" "NKDim" "NKBright" "cM" "ncM")

while read trait; do
  for cell in "${cells[@]}"; do
    f1="/GWAS_catalog/${trait}/${cell}/${trait}_${cell}_C6.csv.merge.sort.sig_ok"
    f2="/GWAS_catalog/${trait}/${cell}_splicing/${trait}_${cell}_C6_splicing_ok_bef.txt"
    if [[ -f $f1 && -f $f2 ]]; then
      echo "==> ${trait}_${cell}"
      awk 'NR==FNR{a[$2]; next} $2 in a{print $2}' "$f2" "$f1"
    fi
  done
done < "$trait_file" > all1.txt




INPUT_FILE="/all1.txt"
REF_FILE="/gene_type_id_region.txt"
OUTPUT_FILE="/all1_with_info.txt"


awk -v ref_file="${REF_FILE}" '
BEGIN {
    while ((getline < ref_file) > 0) {
        gene_name = $NF
        gene_info[gene_name] = $0
    }
    close(ref_file)
}
{
    if ($0 ~ /^==>/) {
        print $0
    } else {

        gene_name = $1

        if (gene_name in gene_info) {
            print gene_info[gene_name]
        } else {

            print "# Gene not found: " $0
        }
    }
}
' ${INPUT_FILE} > ${OUTPUT_FILE}