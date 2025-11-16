


cat run_pipeline_pred_r_2.sh


TYPE=$1
BASEDIR=$2

mkdir -p ${BASEDIR}/${TYPE} && cd ${BASEDIR}/${TYPE}
mkdir -p model && cd model
cp ../scale_data_3_${TYPE}.txt ./scale_data_3.txt
cp ../gencode.v12.genes_${TYPE}.gtf  ./gencode.v12.genes.gtf
cp  ../geuvadis.annot_${TYPE}.txt ./geuvadis.annot.txt
cp  ../geuvadis_${TYPE}.snps.txt ./geuvadis.snps.txt
mv scale_data_3.txt geuvadis.expr.txt

mkdir P && cd P
for i in $(seq 1 9); do
mkdir X${i}
cd ./X${i}
mkdir data
mkdir data/input
mkdir data/input/annotations
mkdir data/input/annotations/gene_annotation
mkdir data/input/annotations/snp_annotation
mkdir data/input/expression_phenotypes
mkdir data/input/genotypes

mkdir data/intermediate
mkdir -p data/intermediate/annotations/gene_annotation
mkdir -p data/intermediate/annotations/snp_annotation
mkdir -p data/intermediate/expression_phenotypes
mkdir -p data/intermediate/genotypes
mkdir -p data/intermediate/model_by_chr

mkdir data/output
mkdir -p data/output/allBetas
mkdir -p data/output/allCovariances
mkdir -p data/output/allLogs
mkdir -p data/output/allMetaData
mkdir -p data/output/allResults
mkdir -p data/output/dbs

mkdir joblogs
cp -r /joblogs/example ./joblogs
cp -r /PredictDBPipeline/scripts ./


rm -r ./joblogs/example/predictdb_example
rm    ./joblogs/example/predictdb_example.tar.gz
rm    ./joblogs/example/sh_list/log/*

cp ../../geuvadis.annot.txt       ./data/input/annotations/snp_annotation
cp ../../gencode.v12.genes.gtf    ./data/input/annotations/gene_annotation
cp ../../geuvadis.snps.txt        ./data/input/genotypes
cp ../../geuvadis.expr.txt        ./data/input/expression_phenotypes
cd ../
done

cd ..
cd ./P/X1/joblogs/example
python preprocess_example.py
cd ../../../

for i in $(seq 2 9); do
cp ./X1/data/intermediate/annotations/gene_annotation/* ./X${i}/data/intermediate/annotations/gene_annotation
cp ./X1/data/intermediate/annotations/snp_annotation/* ./X${i}/data/intermediate/annotations/snp_annotation
cp ./X1/data/intermediate/expression_phenotypes/*  ./X${i}/data/intermediate/expression_phenotypes
cp ./X1/data/intermediate/genotypes/*  ./X${i}/data/intermediate/genotypes
done

cd ${BASEDIR}/${TYPE}/model/P/X5/joblogs/example


for i in $(seq 1 22)
do
cd ${BASEDIR}/${TYPE}/model/P/X5/joblogs/example
Rscript splicing.R ../../data/intermediate/expression_phenotypes/geuvadis.expr.RDS ../../data/intermediate/genotypes/geuvadis.snps.chr${i}.txt ../../data/intermediate/annotations/gene_annotation/gencode.v12.genes.parsed.RDS ../../data/intermediate/annotations/snp_annotation/geuvadis.annot.chr${i}.RDS 10 0.5 ../../data/intermediate/model_by_chr/ ${i} HapMap 1e6 gEUVADIS
done



bash run_pipeline_pred_r_2.sh NKDim 2 > NKDim_run_2_splicing.log 2>&1