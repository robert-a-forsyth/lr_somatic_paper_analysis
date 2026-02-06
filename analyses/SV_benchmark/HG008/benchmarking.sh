cd /staging/leuven/stg_00096/home/rforsyth/LRSomatic/CCS15_LRSomatic_Testing/gr38_CCS15/sv_benchmarking/HG008
conda activate minda

RAW_HG008_ONT=lr_somatic_paper_analysis/analyses/SV_benchmark/data/HG008-ONT/severus_somatic.vcf.gz
RAW_HG008_ONT_TO=lr_somatic_paper_analysis/analyses/SV_benchmark/data/HG008-ONT-TO/severus_somatic.vcf.gz
RAW_HG008_PB=lr_somatic_paper_analysis/analyses/SV_benchmark/data/HG008-PB/severus_somatic.vcf.gz
RAW_HG008_PB_TO=lr_somatic_paper_analysis/analyses/SV_benchmark/data/HG008-PB-TO/severus_somatic.vcf.gz

TRUTH_SET_RAW=lr_somatic_paper_analysis/analyses/SV_benchmark/HG008/GRCh38_HG008-T-V0.4_somatic-stvar_PASS.draftbenchmark.vcf.gz

VAF_HG008_ONT=lr_somatic_paper_analysis/analyses/SV_benchmark/HG008/vaf_vcfs/ONT_vaf.vcf.gz
VAF_HG008_ONT_TO=lr_somatic_paper_analysis/analyses/SV_benchmark/HG008/vaf_vcfs/ONT_TO_vaf.vcf.gz
VAF_HG008_PB=lr_somatic_paper_analysis/analyses/SV_benchmark/HG008/vaf_vcfs/PB_vaf.vcf.gz
VAF_HG008_PB_TO=lr_somatic_paper_analysis/analyses/SV_benchmark/HG008/vaf_vcfs/PB_TO_vaf.vcf.gz
BED=lr_somatic_paper_analysis/analyses/SV_benchmark/HG008/GRCh38_HG008-T-V0.4_somatic-stvar.draftbenchmark.bed

bcftools filter -i "VAF>.1" $RAW_HG008_ONT_TO -Oz -o $VAF_HG008_ONT_TO
bcftools filter -i "VAF>.1" $RAW_HG008_ONT -Oz -o $VAF_HG008_ONT
bcftools filter -i "VAF>.1" $RAW_HG008_PB -Oz -o $VAF_HG008_PB
bcftools filter -i "VAF>.1" $RAW_HG008_PB_TO -Oz -o $VAF_HG008_PB_TO

tabix -p vcf $VAF_HG008_ONT_TO
tabix -p vcf $VAF_HG008_ONT
tabix -p vcf $VAF_HG008_PB
tabix -p vcf $VAF_HG008_PB_TO

REGION_HG008_ONT=lr_somatic_paper_analysis/analyses/SV_benchmark/HG008/vaf_vcfs/ONT_vaf_region.vcf.gz
REGION_HG008_ONT_TO=lr_somatic_paper_analysis/analyses/SV_benchmark/HG008/vaf_vcfs/ONT_TO_vaf_region.vcf.gz
REGION_HG008_PB=lr_somatic_paper_analysis/analyses/SV_benchmark/HG008/vaf_vcfs/PB_vaf_region.vcf.gz
REGION_HG008_PB_TO=lr_somatic_paper_analysis/analyses/SV_benchmark/HG008/vaf_vcfs/PB_TO_vaf_region.vcf.gz

bcftools view -R $BED $VAF_HG008_ONT_TO -Oz -o $REGION_HG008_ONT_TO
bcftools view -R $BED $VAF_HG008_ONT -Oz -o $REGION_HG008_ONT
bcftools view -R $BED $VAF_HG008_PB -Oz -o $REGION_HG008_PB
bcftools view -R $BED $VAF_HG008_PB_TO -Oz -o $REGION_HG008_PB_TO


TRUTH_SET=filtered_truth.vcf.gz

bcftools view -R $BED $TRUTH_SET_RAW -Oz -o $TRUTH_SET

### RUN BENCHMARK

$MINDA_PATH truthset --base $TRUTH_SET_RAW --vcfs $VAF_HG008_ONT $VAF_HG008_ONT_TO $VAF_HG008_PB $VAF_HG008_PB_TO --min_size 50 --tolerance 500 --out_dir minda_out

$MINDA_PATH truthset --base $TRUTH_SET --vcfs $REGION_HG008_ONT $REGION_HG008_ONT_TO $REGION_HG008_PB $REGION_HG008_PB_TO --min_size 50 --tolerance 500 --out_dir minda_out_2
