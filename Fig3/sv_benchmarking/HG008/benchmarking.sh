cd /staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/HG008
source /staging/leuven/stg_00096/home/rforsyth/miniconda3/etc/profile.d/conda.sh
conda activate minda
MINDA_PATH=/staging/leuven/stg_00096/home/rforsyth/software/minda/minda.py
RAW_HG008_ONT=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/data/lr_somatic/out/HG008-ONT/variants/severus/somatic_SVs/severus_somatic.vcf.gz
RAW_HG008_ONT_TO=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/data/lr_somatic/out/HG008-ONT-TO/variants/severus/somatic_SVs/severus_somatic.vcf.gz
RAW_HG008_PB=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/data/lr_somatic/out/HG008-PB/variants/severus/somatic_SVs/severus_somatic.vcf.gz
RAW_HG008_PB_TO=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/data/lr_somatic/out/HG008-PB-TO/variants/severus/somatic_SVs/severus_somatic.vcf.gz

TRUTH_SET_RAW=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/GRCh38_HG008-T-V0.4_somatic-stvar_PASS.draftbenchmark.vcf.gz

VAF_HG008_ONT=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/HG008/ONT_vaf.vcf.gz
VAF_HG008_ONT_TO=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/HG008/ONT_TO_vaf.vcf.gz
VAF_HG008_PB=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/HG008/PB_vaf.vcf.gz
VAF_HG008_PB_TO=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/HG008/PB_TO_vaf.vcf.gz
BED=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/GRCh38_HG008-T-V0.4_somatic-stvar.draftbenchmark.bed

#bcftools view -R $BED $RAW_HG008_ONT_TO | bcftools filter -i "VAF>.1" -Oz -o $VAF_HG008_ONT_TO
#bcftools view -R $BED $RAW_HG008_ONT | bcftools filter -i "VAF>.1" -Oz -o $VAF_HG008_ONT
#bcftools view -R $BED $RAW_HG008_PB | bcftools filter -i "VAF>.1" -Oz -o $VAF_HG008_PB
#bcftools view -R $BED $RAW_HG008_PB_TO | bcftools filter -i "VAF>.1" -Oz -o $VAF_HG008_PB_TO

bcftools filter -i "VAF>.1" $RAW_HG008_ONT_TO -Oz -o $VAF_HG008_ONT_TO
bcftools filter -i "VAF>.1" $RAW_HG008_ONT -Oz -o $VAF_HG008_ONT
bcftools filter -i "VAF>.1" $RAW_HG008_PB -Oz -o $VAF_HG008_PB
bcftools filter -i "VAF>.1" $RAW_HG008_PB_TO -Oz -o $VAF_HG008_PB_TO

tabix -p vcf $VAF_HG008_ONT_TO
tabix -p vcf $VAF_HG008_ONT
tabix -p vcf $VAF_HG008_PB
tabix -p vcf $VAF_HG008_PB_TO

REGION_HG008_ONT=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/HG008/ONT_vaf_region.vcf.gz
REGION_HG008_ONT_TO=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/HG008/ONT_TO_vaf_region.vcf.gz
REGION_HG008_PB=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/HG008/PB_vaf_region.vcf.gz
REGION_HG008_PB_TO=/staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/HG008/PB_TO_vaf_region.vcf.gz

bcftools view -R $BED $VAF_HG008_ONT_TO -Oz -o $REGION_HG008_ONT_TO
bcftools view -R $BED $VAF_HG008_ONT -Oz -o $REGION_HG008_ONT
bcftools view -R $BED $VAF_HG008_PB -Oz -o $REGION_HG008_PB
bcftools view -R $BED $VAF_HG008_PB_TO -Oz -o $REGION_HG008_PB_TO


TRUTH_SET=filtered_truth.vcf.gz

bcftools view -R $BED $TRUTH_SET_RAW -Oz -o $TRUTH_SET

### RUN BENCHMARK

$MINDA_PATH truthset --base $TRUTH_SET_RAW --vcfs $VAF_HG008_ONT $VAF_HG008_ONT_TO $VAF_HG008_PB $VAF_HG008_PB_TO --min_size 50 --tolerance 500 --out_dir minda_out

$MINDA_PATH truthset --base $TRUTH_SET --vcfs $REGION_HG008_ONT $REGION_HG008_ONT_TO $REGION_HG008_PB $REGION_HG008_PB_TO --min_size 50 --tolerance 500 --out_dir minda_out_2

### GET TSV FILES
conda activate vcf_tsv
vcf2tsvpy --input_vcf $VAF_HG008_ONT --out_tsv HG008_ONT_vaf.tsv
sed -i '1d' HG008_ONT_vaf.tsv

vcf2tsvpy --input_vcf $VAF_COLO_ONT_TO --out_tsv HG008_ONT_TO_vaf.tsv
sed -i '1d' HG008_ONT_TO_vaf.tsv

vcf2tsvpy --input_vcf $VAF_HG008_PB --out_tsv HG008_PB_vaf.tsv
sed -i '1d' HG008_PB_vaf.tsv

vcf2tsvpy --input_vcf $VAF_HG008_PB_TO --out_tsv HG008_PB_TO_vaf.tsv
sed -i '1d' HG008_PB_TO_vaf.tsv
