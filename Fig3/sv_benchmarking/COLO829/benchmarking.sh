cd /staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/COLO829
source /staging/leuven/stg_00096/home/rforsyth/miniconda3/etc/profile.d/conda.sh
conda activate minda
MINDA_PATH=/staging/leuven/stg_00096/home/rforsyth/software/minda/minda.py
RAW_COLO_ONT=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/data/lr_somatic/out/COLO829-ONT/variants/severus/somatic_SVs/severus_somatic.vcf.gz
RAW_COLO_ONT_TO=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/data/lr_somatic/out/COLO829-ONT-TO/variants/severus/somatic_SVs/severus_somatic.vcf.gz
RAW_COLO_PB=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/data/lr_somatic/out/COLO829-PB/variants/severus/somatic_SVs/severus_somatic.vcf.gz
RAW_COLO_PB_TO=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/data/lr_somatic/out/COLO829-PB-TO/variants/severus/somatic_SVs/severus_somatic.vcf.gz

TRUTH_SET=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/severus_files/variant_calls_and_benchmarks/truthsets_beds/truthset_somaticSVs_COLO829_hg38lifted_chr_alt_vaf.vcf.gz

### ADD VAF TO RAW FILES
VAF_SCRIPT_PATH=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/add_vaf.py

VAF_SEVERUS_COMP=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/severus_files/variant_calls_and_benchmarks/COLO829/vaf_vcfs/severus_pb.vcf.gz
VAF_COLO_ONT=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/COLO829/ONT_vaf.vcf.gz
VAF_COLO_ONT_TO=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/COLO829/ONT_TO_vaf.vcf.gz
VAF_COLO_PB=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/COLO829/PB_vaf.vcf.gz
VAF_COLO_PB_TO=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/COLO829/PB_TO_vaf.vcf.gz

python3 $VAF_SCRIPT_PATH $RAW_COLO_ONT $VAF_COLO_ONT
python3 $VAF_SCRIPT_PATH $RAW_COLO_ONT_TO $VAF_COLO_ONT_TO
python3 $VAF_SCRIPT_PATH $RAW_COLO_PB $VAF_COLO_PB
python3 $VAF_SCRIPT_PATH $RAW_COLO_PB_TO $VAF_COLO_PB_TO

### RUN BENCHMARK
$MINDA_PATH truthset --base $TRUTH_SET --vcfs $VAF_SEVERUS_COMP $VAF_COLO_ONT $VAF_COLO_ONT_TO $VAF_COLO_PB $VAF_COLO_PB_TO --min_size 50 --tolerance 500 --out_dir minda_out --vaf 0.1

### GET TSV FILES
conda activate vcf_tsv
vcf2tsvpy --input_vcf $VAF_COLO_ONT --out_tsv COLO_ONT_vaf.tsv
sed -i '1d' COLO_ONT_vaf.tsv

vcf2tsvpy --input_vcf $VAF_COLO_ONT_TO --out_tsv COLO_ONT_TO_vaf.tsv
sed -i '1d' COLO_ONT_TO_vaf.tsv

vcf2tsvpy --input_vcf $VAF_COLO_PB --out_tsv COLO_PB_vaf.tsv
sed -i '1d' COLO_PB_vaf.tsv

vcf2tsvpy --input_vcf $VAF_COLO_PB_TO --out_tsv COLO_PB_TO_vaf.tsv
sed -i '1d' COLO_PB_TO_vaf.tsv


VAF_COLO_ONT2=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/COLO829/ONT_vaf2.vcf.gz
VAF_COLO_ONT_TO2=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/COLO829/ONT_TO_vaf2.vcf.gz
VAF_COLO_PB2=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/COLO829/PB_vaf2.vcf.gz
VAF_COLO_PB_TO2=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/COLO829/PB_TO_vaf2.vcf.gz

bcftools filter -i "INFO/VAF<.9" $VAF_COLO_ONT -Oz -o $VAF_COLO_ONT2
bcftools filter -i "INFO/VAF<.9" $VAF_COLO_ONT_TO -Oz -o $VAF_COLO_ONT_TO2
bcftools filter -i "INFO/VAF<.9" $VAF_COLO_PB -Oz -o $VAF_COLO_PB2
bcftools filter -i "INFO/VAF<.9" $VAF_COLO_PB_TO -Oz -o $VAF_COLO_PB_TO2

conda activate minda

$MINDA_PATH truthset --base $TRUTH_SET --vcfs $VAF_COLO_ONT2 $VAF_COLO_ONT_TO2 $VAF_COLO_PB2 $VAF_COLO_PB_TO2 --min_size 50 --tolerance 500 --out_dir minda_out2 --vaf 0.1
